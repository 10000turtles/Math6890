//
// Solvetheheatequationinone-dimension
//

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

// define anewtype"Real"whichisequivalenttoa"double "
typedef double Real;

#include <string>
using std::max;
using std::string;

#include <ctime>
//--------------------------------------------
// Returnthecurrentwall-clocktimeinseconds
//--------------------------------------------
inline double getCPU()
{
    return (1.0 * std::clock()) / CLOCKS_PER_SEC;
}

//--------------------------------------------------------------------------------------
// Functiontosaveavectortoamatlabfile.
// matlabFile(input):savevectortothisfile
// u_p(input):arrayofvectorvalues
// name(input):nameforarray
//(nd1a:nd1b)(input):arraydimensions
//--------------------------------------------------------------------------------------
int writeMatlabVector(FILE *matlabFile, Real *u_p, const char *name, int nd1a, int nd1b)
{
#define u(i) u_p[i - nd1a]

    const int numPerLine = 8; // numberofentriesperline
    // Savethevectoras:
    // name=[numnumnumnumnum...
    // numnumnumnumnum];
    fprintf(matlabFile, "%s=[", name);
    for (int i = nd1a; i <= nd1b; i++)
    {
        fprintf(matlabFile, "%20.15e", u(i));
        if ((i - nd1a) % numPerLine == numPerLine - 1)
            fprintf(matlabFile, "...\n"); // continuationline
    }
    fprintf(matlabFile, "];\n");

    return 0;
}

int main(int argc, char *argv[])
{
    printf("Usage:heat1d[Nx][matlabFileName.m]\n"
           "Nx=numberofgridcells.\n"
           "matlabFileName.m:saveresultstothisfile.\n");

#define TRIG_DD 1
#define TRIG_NN 2
#define POLY_DD 3
#define POLY_NN 4
    //=====Choosethesolutionhereorcompilewith-DSOLUTION=[1|2|3|4]=====

#define SOLUTION POLY_DD

    const Real pi = M_PI;

    int debug = 0; // setto1fordebuginfo
    Real xa = 0., xb = 1.;
    Real kappa = .1;
    Real tFinal = .2;
    Real cfl = .9; // time-stepsafetyfactor

    int Nx = 10; // default
    string matlabFileName = "heat1d.m";

    if (argc >= 2) // readanycommandlinearguments
    {
        Nx = atoi(argv[1]);
        tFinal = atof(argv[2]);
        printf("SettingNx=%d\n", Nx);
        if (argc >= 3)
        {

            printf("SettingmatlabFileName=[%s]\n", matlabFileName.c_str());
        }
    }

    //=============Gridandindexing==============
    // xaxb
    // G---X---+---+---+---+--...---+---X---G
    // Nx
    // n1an1b
    // nd1and1b
    // Cindex:3...

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

    const char solutionName[] = "polyDD";
    boundaryCondition(0, 0) = dirichlet;
    boundaryCondition(1, 0) = dirichlet;

    Real b0 = 1., b1 = .5, b2 = .25;
    Real a0 = 1., a1 = .3;

#define UTRUE(x, t) (b0 + (x) * (b1 + (x) * b2)) * (a0 + (t) * (a1))
#define UTRUEX(x, t) (b1 + 2. * (x) * b2) * (a0 + (t) * (a1))
#define UTRUET(x, t) (b0 + (x) * (b1 + (x) * b2)) * (a1)
#define UTRUEXX(x, t) (2. * b2) * (a0 + (t) * (a1))

    // force=u_t-kappa*u.xx
#define FORCE(x, t) (UTRUET(x, t) - kappa * UTRUEXX(x, t))

    Real *u_p[2]; // twoarrayswillbeusedforcurrentandnewtimes
    u_p[0] = new Real[nd1];
    u_p[1] = new Real[nd1];

// Macrostodefine fortranlikearrays
#define uc(i) u_p[cur][i - nd1a]
#define un(i) u_p[next][i - nd1a]

    // initialconditions
    Real t = 0.;
    int cur = 0; //"current"solution,indexintou_p[]
    for (int i = nd1a; i <= nd1b; i++)
        uc(i) = UTRUE(x(i), t);

    if (debug > 0)
    {
        printf("Afterinitialconditions\nu=[");
        for (int i = nd1a; i <= nd1b; i++)
            printf("%10.4e,", uc(i));
        printf("]\n");
    }

    // Time-steprestrictioniskappa*dt/dx^2<.5
    const Real dx2 = dx * dx;
    Real dt = cfl * .5 * dx2 / kappa; // dt,adjustedbelow
    const int numSteps = ceil(tFinal / dt);
    dt = tFinal / numSteps; // adjustdttoreachthefinaltime
    const Real rx = kappa * dt / dx2;

    printf("-------------------Solvetheheatequationin1Dsolution=%s---------------------\n",
           solutionName);
    printf("numGhost=%d,n1a=%d,n1b=%d,nd1a=%d,nd1b=%d\n", numGhost, n1a, n1b, nd1a, nd1b);
    printf("numSteps=%d,Nx=%d,kappa=%g,tFinal=%g,boundaryCondition(0,0)=%d,boundaryCondition(1,0)=%d\n",
           numSteps, Nx, kappa, tFinal, boundaryCondition(0, 0), boundaryCondition(1, 0));

    //----------TIME-STEPPINGLOOP--------
    Real cpu0 = getCPU();
    for (int n = 0; n < numSteps; n++)
    {
        t = n * dt; // currenttime

        const int cur = n % 2;        // currenttimelevel
        const int next = (n + 1) % 2; // nexttimelevel

        //---updatetheinteriorpoints---
        for (int i = n1a; i <= n1b; i++)
        {
            un(i) = uc(i) + rx * (uc(i + 1) - 2. * uc(i) + uc(i - 1)) + dt * FORCE(x(i), t);
        }

        //----boundaryconditions---
        for (int side = 0; side <= 1; side++)
        {
            const int i = side == 0 ? n1a : n1b; // boundaryindex
            const int is = 1 - 2 * side;         // is=1onleft,-1onright
            if (boundaryCondition(side, 0) == dirichlet)
            {
                un(i) = UTRUE(x(i), t + dt);
                un(i - is) = 3. * un(i) - 3. * un(i + is) + un(i + 2 * is); // extrapolateghost
            }
            else
            {
                // NeumannBC
                un(i - is) = un(i + is) - 2. * is * dx * UTRUEX(x(i), t + dt);
            }
        }

        if (debug > 1)
        {
            printf("step%d:AfterupdateinteriorandrealBCs\nu=[", n + 1);
            for (int i = nd1a; i <= nd1b; i++)
                printf("%12.4e,", un(i));
            printf("]\n");
        }

        if (debug > 0)
        {
            // computetheerror
            Real maxErr = 0.;
            for (int i = nd1a; i <= nd1b; i++)
            {
                Real err = fabs(un(i) - UTRUE(x(i), t + dt));
                maxErr = max(maxErr, err);
            }
            printf("step=%d,t=%9.3e,maxErr=%9.2e\n", n + 1, t + dt, maxErr);
        }

    } // endtime-steppingloop

    Real cpuTimeStep = getCPU() - cpu0;

    //----checktheerror----
    t += dt; // tFinal;
    if (fabs(t - tFinal) > 1e-3 * dt / tFinal)
    {
        printf("ERROR:AFTERTIME_STEPPING:t=%16.8eISNOTEQUALtotFinal=%16.8e\n", t, tFinal);
    }

    Real *error_p = new Real[nd1];
#define error(i) error_p[i - nd1a]

    cur = numSteps % 2;
    Real maxErr = 0.;
    for (int i = nd1a; i <= nd1b; i++)
    {
        error(i) = uc(i) - UTRUE(x(i), t);
        maxErr = max(maxErr, abs(error(i)));
    }

    printf("numSteps=%4d,Nx=%3d,maxErr=%9.2e,cpu=%9.2e(s)\n", numSteps, Nx, maxErr, cpuTimeStep);

    //---Writeafileforplottinginmatlab--
    FILE *matlabFile = fopen(matlabFileName.c_str(), "w");
    fprintf(matlabFile, "%%Filewrittenbyheat1d.C\n");
    fprintf(matlabFile, "xa=%g;xb=%g;kappa=%g;t=%g;maxErr=%10.3e;cpuTimeStep=%10.3e;\n", xa, xb, kappa, tFinal, maxErr, cpuTimeStep);
    fprintf(matlabFile, "Nx=%d;dx=%14.6e;numGhost=%d;n1a=%d;n1b=%d;nd1a=%d;nd1b=%d;\n", Nx, dx, numGhost, n1a, n1b, nd1a, nd1b);
    fprintf(matlabFile, "solutionName=\'%s\';\n", solutionName);

    writeMatlabVector(matlabFile, x_p, "x", nd1a, nd1b);
    writeMatlabVector(matlabFile, u_p[cur], "u", nd1a, nd1b);
    writeMatlabVector(matlabFile, error_p, "err", nd1a, nd1b);

    fclose(matlabFile);
    printf("Wrotefile%s\n\n", matlabFileName.c_str());

    delete[] u_p[0];
    delete[] u_p[1];
    delete[] x_p;
    delete[] error_p;
    delete[] boundaryCondition_p;

    return 0;
}