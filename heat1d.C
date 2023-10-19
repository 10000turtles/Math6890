
// Solvetheheatequationinone-dimension
//

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <omp.h>

// defineanewtype"Real"whichisequivalenttoa"double"
typedef double Real;

#include <string>

using std::max;
using std::string;
using namespace std;

#include <ctime>
//--------------------------------------------
// Returnthecurrentwall-clocktimeinseconds
//--------------------------------------------
inline double getCPU()
{
#if defined(_OPENMP)
    return omp_get_wtime();
#else
    return (1.0 * std::clock()) / CLOCKS_PER_SEC;
#endif
}

//--------------------------------------------------------------------------------------
// Functiontosaveavectortoamatlabfile.
// matlabFile (input):savevectortothisfile
// u_p (input):arrayofvectorvalues
// name (input):nameforarray
// (nd1a:nd1b)(input):arraydimensions
//--------------------------------------------------------------------------------------
int writeMatlabVector(FILE *matlabFile, Real *u_p, const char *name, int nd1a, int nd1b)
{
#define u(i) u_p[i - nd1a]

    const int numPerLine = 8; // numberofentriesperline
    // Savethevectoras:
    //  name=[numnumnumnumnum...
    //  numnumnumnumnum];
    fprintf(matlabFile, "%s=[", name);
    for (int i = nd1a; i <= nd1b; i++)
    {
        fprintf(matlabFile, "%20.15e ", u(i));
        if ((i - nd1a) % numPerLine == numPerLine - 1)
            fprintf(matlabFile, "...\n"); // continuationline
    }
    fprintf(matlabFile, "];\n");

    return 0;
}

int main(int argc, char *argv[])
{
    printf("Usage: heat1d [Nx] [num_threads] [matlabFileName.m]\n"
           " Nx = number of grid cells.\n"
           " matlabFileName.m : save results to this file.\n");

#define TRIG_DD 1
#define TRIG_NN 2
#define POLY_DD 3
#define POLY_NN 4

//=====Choosethesolutionhereorcompilewith-DSOLUTION=[1|2|3|4]=====
#ifndef SOLUTION
#define SOLUTION TRIG_DD
// #define SOLUTION TRIG_NN
//  #define SOLUTION POLY_DD
//  #define SOLUTION POLY_NN
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
    const Real pi = M_PI;

    int debug = 0; // setto1fordebuginfo
    Real xa = 0., xb = 1.;
    Real kappa = .1;
    Real tFinal = .01;
    Real cfl = .9; // time-stepsafetyfactor

    int Nx = 10; // default

    string matlabFileName = "heat1d" + to_string(Nx) + ".m";

    int num_th = 1;

    if (argc >= 2) // read any command line arguments
    {
        Nx = atoi(argv[1]);
        matlabFileName = "heat1d" + to_string(Nx) + ".m";
        printf("Setting Nx=%d\n", Nx);
        if (argc >= 4)
        {
            matlabFileName = argv[3];
            printf("Setting matlabFileName=[%s]\n", matlabFileName.c_str());
        }
        num_th = atoi(argv[2]);
    }

    //=============Gridandindexing==============
    // xa xb
    // G---X---+---+---+---+--...---+---X---G
    // 0 1 2 Nx
    // n1a n1b
    // nd1a nd1b
    // Cindex: 0 1 2 3...

    Real dx = (xb - xa) / Nx;
    const int numGhost = 1;
    const int n1a = 0;
    const int n1b = Nx;
    const int nd1a = n1a - numGhost;
    const int nd1b = n1b + numGhost;
    const int nd1 = nd1b - nd1a + 1; // totalnumberofgridpoints;

    // Createanarrayofgridpoints:
    Real *x_p = new Real[nd1];
#define x(i) x_p[i - nd1a]

    for (int i = nd1a; i <= nd1b; i++)
        x(i) = xa + (i - n1a) * dx;

    if (debug > 1)
    {
        for (int i = nd1a; i <= nd1b; i++)
            printf("x(%2d)=%12.4e\n", i, x(i));
    }

    const int dirichlet = 1, neumann = 2;
    const int numberOfDimensions = 1;
    int *boundaryCondition_p = new int[2 * numberOfDimensions];
#define boundaryCondition(side, axis) boundaryCondition_p[(side) + 2 * (axis)]

    const Real kx = 3.;
    const Real kxPi = kx * pi;
    const Real kappaPiSq = kappa * kxPi * kxPi;

#if SOLUTION == TRIG_DD
    // TruesolutionfordirichletBC’s
    boundaryCondition(0, 0) = dirichlet;
    boundaryCondition(1, 0) = dirichlet;

    const char solutionName[] = "trueDD";

#define UTRUE(x, t) sin(kxPi *(x)) * exp(-kappaPiSq *(t))
#define UTRUEX(x, t) kxPi *cos(kxPi *(x)) * exp(-kappaPiSq *(t))
#define FORCE(x, t) (0.)

#elif SOLUTION == TRIG_NN

    // TruesolutionforNeumannBC’s
    boundaryCondition(0, 0) = neumann;
    boundaryCondition(1, 0) = neumann;
    const char solutionName[] = "trueNN";

#define UTRUE(x, t) cos(kxPi *(x)) * exp(-kappaPiSq *(t))
#define UTRUEX(x, t) -kxPi *sin(kxPi *(x)) * exp(-kappaPiSq *(t))
#define FORCE(x, t) (0.)

#elif (SOLUTION == POLY_DD) || (SOLUTION == POLY_NN)

// polynomialmanufacturedsolution
#if SOLUTION == POLY_DD
    const char solutionName[] = "polyDD";
    boundaryCondition(0, 0) = dirichlet;
    boundaryCondition(1, 0) = dirichlet;
#else
    const char solutionName[] = "polyNN";
    boundaryCondition(0, 0) = neumann;
    boundaryCondition(1, 0) = neumann;
#endif

    const Real b0 = 1., b1 = .5, b2 = .25;
    const Real a0 = 1., a1 = .3;
#define UTRUE(x, t) (b0 + (x) * (b1 + (x) * b2)) * (a0 + (t) * (a1))
#define UTRUEX(x, t) (b1 + 2. * (x) * b2) * (a0 + (t) * (a1))
#define UTRUET(x, t) (b0 + (x) * (b1 + (x) * b2)) * (a1)
#define UTRUEXX(x, t) (2. * b2) * (a0 + (t) * (a1))

    // force=u_t-kappa*u.xx
#define FORCE(x, t) (UTRUET(x, t) - kappa * UTRUEXX(x, t))

#else
    printf("ERROR: unknown solution");
    abort();
#endif

    Real *u_p[2]; // twoarrayswillbeusedforcurrentandnewtimes
    u_p[0] = new Real[nd1];
    u_p[1] = new Real[nd1];

// Macrostodefinefortranlikearrays
#define uc(i) u_p[cur][i - nd1a]
#define un(i) u_p[next][i - nd1a]
#define w(i) w_temp[i - nd1a]
    // initialconditions
    Real t = 0.;
    int cur = 0; //"current"solution,indexintou_p[]
    for (int i = nd1a; i <= nd1b; i++)
        uc(i) = UTRUE(x(i), t);

    if (debug > 0)
    {
        printf("After initial conditions\nu=[");
        for (int i = nd1a; i <= nd1b; i++)
            printf("%10.4e,", uc(i));
        printf("]\n");
    }

    // Time-steprestrictionis kappa*dt/dx^2<.5
    const Real dx2 = dx * dx;
    Real dt = cfl * .5 * dx2 / kappa; // dt,adjustedbelow
    const int numSteps = ceil(tFinal / dt);
    dt = tFinal / numSteps; // adjustdttoreachthefinaltime
    const Real rx = kappa * dt / dx2;

    printf("-------------------Solve the heat equation in 1D solution=%s---------------------\n",
           solutionName);
    printf(" numGhost=%d, n1a=%d, n1b=%d, nd1a=%d, nd1b=%d\n", numGhost, n1a, n1b, nd1a, nd1b);
    printf(" numSteps=%d, Nx=%d, kappa=%g, tFinal=%g, boundaryCondition[0]=%d, boundaryCondition[1]=%d\n",
           numSteps, Nx, kappa, tFinal, boundaryCondition(0, 0), boundaryCondition(1, 0));

    //----------TIME-STEPPINGLOOP--------
    Real cpu0 = getCPU();
    Real maxErr = 0.;
    int i;
    // #pragma omp parallel default(shared) num_threads(num_th)
    for (int n = 0; n < numSteps; n++)
    {
        t = n * dt; // currenttime

        const int cur = n % 2;        // currenttimelevel
        const int next = (n + 1) % 2; // nexttimelevel
        Real *w_temp = new Real[nd1];
//---updatetheinteriorpoints---
#pragma omp parallel default(shared) num_threads(num_th)
        {
#pragma omp for private(i)
            for (i = n1a; i <= n1b; i++)
            {
                un(i) = uc(i) + rx * (uc(i + 1) - 2. * uc(i) + uc(i - 1)) + dt * FORCE(x(i), t);
                // w(i) = uc(i) + rx / 2 * (uc(i + 1) - 2. * uc(i) + uc(i - 1)) + dt / 2 * FORCE(x(i), t);
            }
        }
        //----boundaryconditions---
        for (int side = 0; side <= 1; side++)
        {
            const int i = side == 0 ? n1a : n1b; // boundaryindex
            const int is = 1 - 2 * side;         // is=1onleft,-1onright
            if (boundaryCondition(side, 0) == dirichlet)
            {
                w(i) = UTRUE(x(i), t + dt);
                w(i - is) = 3. * w(i) - 3. * w(i + is) + w(i + 2 * is); // extrapolateghost
            }
            else
            {
                // NeumannBC
                w(i - is) = w(i + is) - 2. * is * dx * UTRUEX(x(i), t + dt);
            }
        }

        if (debug > 1)
        {
            printf("step %d: After update interior and real BCs\nu=[", n + 1);
            for (int i = nd1a; i <= nd1b; i++)
                printf("%12.4e, ", un(i));
            printf("]\n");
        }

        if (debug > 0)
        {
            // computetheerror
            maxErr = 0.;
            // #pragma omp parallel  default(shared) num_threads(num_th)  private(i)
            // #pragma omp parallel for private(i) reduction(max : maxErr)
#pragma omp parallel default(shared) num_threads(num_th)
            {
#pragma omp for private(i)
                for (int i = nd1a; i <= nd1b; i++)
                {
                    Real err = fabs(un(i) - UTRUE(x(i), t + dt));
                    maxErr = max(maxErr, err);
                }
            }
            printf("step=%d, t=%9.3e, maxErr=%9.2e\n", n + 1, t + dt, maxErr);
        }

    } // endtime-steppingloop

    Real cpuTimeStep = getCPU() - cpu0;

    //----checktheerror----
    t += dt; // tFinal;
    if (fabs(t - tFinal) > 1e-3 * dt / tFinal)
    {
        printf("ERROR: AFTER TIME_STEPPING: t=%16.8e IS NOT EQUAL to tFinal=%16.8e\n", t, tFinal);
    }

    Real *error_p = new Real[nd1];
#define error(i) error_p[i - nd1a]

    cur = numSteps % 2;
    maxErr = 0.;
    for (int i = nd1a; i <= nd1b; i++)
    {
        error(i) = uc(i) - UTRUE(x(i), t);
        maxErr = max(maxErr, abs(error(i)));
    }

    printf("numSteps=%4d, Nx=%3d, maxErr=%9.2e, cpu=%9.2e(s)\n", numSteps, Nx, maxErr, cpuTimeStep);

    //---Writeafileforplottinginmatlab--
    FILE *matlabFile = fopen(matlabFileName.c_str(), "w");
    fprintf(matlabFile, "%%Fileswrittenbyheat1d.C\n");
    fprintf(matlabFile, "xa=%g;xb=%g;kappa=%g;t=%g;maxErr%d=%10.3e;cpuTimeStep=%10.3e;\n", xa, xb, kappa, tFinal, Nx, maxErr, cpuTimeStep);
    fprintf(matlabFile, "Nx=%d;dx%d=%14.6e;numGhost=%d;n1a=%d;n1b=%d;nd1a=%d;nd1b=%d;\n", Nx, Nx, dx, numGhost, n1a, n1b, nd1a, nd1b);
    fprintf(matlabFile, "solutionName=\'%s\';\n", solutionName);

    string s1 = "x" + to_string(Nx);
    string s2 = "u" + to_string(Nx);
    string s3 = "err" + to_string(Nx);

    writeMatlabVector(matlabFile, x_p, s1.c_str(), nd1a, nd1b);
    writeMatlabVector(matlabFile, u_p[cur], s2.c_str(), nd1a, nd1b);
    writeMatlabVector(matlabFile, error_p, s3.c_str(), nd1a, nd1b);

    fclose(matlabFile);
    printf("Wrote file %s\n\n", matlabFileName.c_str());

    delete[] u_p[0];
    delete[] u_p[1];
    delete[] x_p;
    delete[] error_p;
    delete[] boundaryCondition_p;

    return 0;
}