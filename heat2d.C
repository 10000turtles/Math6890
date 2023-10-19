//=====================================================================
//
// HEATEQUATIONINTWODIMENSIONS
// SolvewithA++arrays
//
//=====================================================================
#include "A++.h"

// define sometypes
typedef double Real;
typedef doubleSerialArray RealArray;
typedef intSerialArray IntegerArray;

#include <float.h>
#include <limits.h>
#include <ctime>
#include <omp.h>
#define REAL_EPSILON DBL_EPSILON
#define REAL_MIN DBL_MIN

// getCPU():Returnthecurrentwall-clocktimeinseconds
inline double getCPU()
{
#if defined(_OPENMP)
    return omp_get_wtime();
#else
    return (1.0 * std::clock()) / CLOCKS_PER_SEC;
#endif
}

// include commandstpparsecommandlinearguments
#include "parseCommand.h"

// functiontowriteanarraytoamatlabreabablefile:
#include "writeMatlabArray.h"

//---Declarethefortranroutineasa"C"function
// Somecompilersaddanunder-scoretothenameof"C"andFortranroutines
// Notealsothatthefortrannameisalllowercase
#define heat2dUpdate heat2dupdate_
extern "C"
{
    void heat2dUpdate(const int &n1a, const int &n1b, const int &n2a, const int &n2b,
                      const int &nd1a, const int &nd1b, const int &nd2a, const int &nd2b,
                      Real &un, const Real &u, const Real &rx, const Real &ry, const int &thr);
}

int main(int argc, char *argv[])
{

    const Real pi = 4. * atan2(1., 1.);

    ios::sync_with_stdio();    // SynchronizeC++andCI/Osubsystems
    Index::setBoundsCheck(on); // TurnonA++arrayboundschecking

    printf("Usage:heat2d-nx=<i>-option=[0|1|2|3]-tFinal=<f>-debug=<i>-saveMatlab=[0|1|2]matlabFile=<s>\n"
           "option:0=scalarIndexing,1=arrayIndexing,2=cIndexing,3=fortranRoutine\n");

    enum BoundaryConditionsEnum
    {
        periodic = -1,
        dirichlet = 1,
        neumann = 2
    };

    enum OptionsEnum
    {
        scalarIndexing = 0,
        arrayIndexing = 1,
        cIndexing = 2,
        fortranRoutine = 3
    };

    int option = scalarIndexing;

    const int numberOfDimensions = 2;

    int debug = 0;
    Real kappa = .1;
    Real xa = 0., xb = 1.; // domainis[xa,xb]X[ya,yb]
    Real ya = 0., yb = 1.;

    Real tFinal = .5;
    int nx = 100, ny = nx;
    int thread = 1;

    int saveMatlab = 0; // 1=saveamatlabfile,2=savesolutiontoo
    string matlabFileName = "heat2d.m";

    string line;
    for (int i = 1; i < argc; i++)
    {
        line = argv[i];
        // printf("Input:argv[%d]=[%s]\n",i,line.c_str());
        if (parseCommand(line, "-nx=", nx))
        {
            ny = nx;
        }
        else if (parseCommand(line, "-debug=", debug))
        {
        }
        else if (parseCommand(line, "-option=", option))
        {
        }
        else if (parseCommand(line, "-tFinal=", tFinal))
        {
        }
        else if (parseCommand(line, "-saveMatlab=", saveMatlab))
        {
        }
        else if (parseCommand(line, "-matlabFileName=", matlabFileName))
        {
        }
        else if (parseCommand(line, "-threads=", thread))
        {
        }
    }

    const int numGhost = 1;
    const int n1a = 0;
    const int n1b = n1a + nx;
    const int nd1a = n1a - numGhost;
    const int nd1b = n1b + numGhost;
    const int nd1 = nd1b - nd1a + 1;

    const int n2a = 0;
    const int n2b = n2a + ny;
    const int nd2a = n2a - numGhost;
    const int nd2b = n2b + numGhost;
    const int nd2 = nd2b - nd2a + 1;

    IntegerArray gridIndexRange(2, numberOfDimensions);
    IntegerArray dimension(2, numberOfDimensions);
    IntegerArray boundaryCondition(2, numberOfDimensions);

    gridIndexRange(0, 0) = n1a;
    gridIndexRange(1, 0) = n1b;
    gridIndexRange(0, 1) = n2a;
    gridIndexRange(1, 1) = n2b;

    dimension(0, 0) = nd1a;
    dimension(1, 0) = nd1b;
    dimension(0, 1) = nd2a;
    dimension(1, 1) = nd2b;

    boundaryCondition(0, 0) = dirichlet; // left
    boundaryCondition(1, 0) = dirichlet; // right

    boundaryCondition(0, 1) = dirichlet; // bottom
    boundaryCondition(1, 1) = dirichlet; // top

    // Gridpoints
    Range Rx(nd1a, nd1b), Ry(nd2a, nd2b);
    RealArray x(Rx, Ry, 2);

    Real dx[2];
    dx[0] = (xb - xa) / nx;
    dx[1] = (yb - ya) / ny;

    int i1, i2;
#pragma omp parallel default(shared) num_threads(thread)
    {
#pragma omp for private(i2, i1)
        for (i2 = nd2a; i2 <= nd2b; i2++)
            for (i1 = nd1a; i1 <= nd1b; i1++)
            {
                x(i1, i2, 0) = xa + (i1 - n1a) * dx[0];
                x(i1, i2, 1) = ya + (i2 - n2a) * dx[1];
            }
    }
    const Real kx = 2., ky = 3;
    const Real kxp = kx * pi;
    const Real kyp = ky * pi;
#define UTRUE(x, y, t) sin(kxp *(x)) * sin(kyp *(y)) * exp(-kappa *(kxp * kxp + kyp * kyp) * (t))

    string optionName = option == scalarIndexing ? "scalarIndexing" : option == arrayIndexing ? "arrayIndexing"
                                                                  : option == cIndexing       ? "cIndexing"
                                                                  : option == fortranRoutine  ? "fortranRoutine"
                                                                                              : "unknown";

    printf("-----SolvetheHeatEquationintwodimensions------\n");
    printf("option=%d:%s\n", option, optionName.c_str());
    printf("saveMatlab=%d,matlabFileName=%s\n", saveMatlab, matlabFileName.c_str());
    printf("kappa=%.3g,nx=%d,ny=%d,tFinal=%6.2f,kx=%g,ky=%g\n", kappa, nx, ny, tFinal, kx, ky);

    // westoretwotimelevels
    RealArray ua[2];
    ua[0].redim(Rx, Ry);
    ua[0] = 0.;
    ua[1].redim(Rx, Ry);
    ua[1] = 0.;

    // initialconditions
    RealArray &u0 = ua[0];
    Real t = 0.;
#pragma omp parallel default(shared) num_threads(thread)
    {
#pragma omp for private(i2, i1)
        for (i2 = nd2a; i2 <= nd2b; i2++)
            for (i1 = nd1a; i1 <= nd1b; i1++)
            {
                u0(i1, i2) = UTRUE(x(i1, i2, 0), x(i1, i2, 1), t);
            }
    }
    // Time-steprestriction:
    // ForwardEuler:kappa*dt*(1/dx^2+1/dy^2)<cfl*.5
    Real cfl = .9;
    Real dt = cfl * (.5 / kappa) / (1. / (dx[0] * dx[0]) + 1. / (dx[1] * dx[1]));

    int numSteps = ceil(tFinal / dt);
    dt = tFinal / numSteps; // adjustdttoreachthefinaltime

    Real rx = kappa * (dt / (dx[0] * dx[0]));
    Real ry = kappa * (dt / (dx[1] * dx[1]));
    int cur = 0;
    int i, n;
    Index I1 = Range(n1a, n1b);
    Index I2 = Range(n1a, n1b);

    Real cpu1 = getCPU();
    for (n = 0; n < numSteps; n++)
    {
        t = n * dt; // curtime

        int next = (cur + 1) % 2;
        doubleArray &u = ua[cur];
        doubleArray &un = ua[next];

        if (option == scalarIndexing)
        {
#pragma omp parallel default(shared) num_threads(thread)
            {
#pragma omp for private(i2, i1)
                for (i2 = n2a; i2 <= n2b; i2++)
                    for (i1 = n1a; i1 <= n1b; i1++)
                    {
                        un(i1, i2) = u(i1, i2) + rx * (u(i1 + 1, i2) - 2. * u(i1, i2) + u(i1 - 1, i2)) + ry * (u(i1, i2 + 1) - 2. * u(i1, i2) + u(i1, i2 - 1));
                    }
            }
        }
        else if (option == arrayIndexing)
        {
            un(I1, I2) = u(I1, I2) + rx * (u(I1 + 1, I2) - 2. * u(I1, I2) + u(I1 - 1, I2)) + ry * (u(I1, I2 + 1) - 2. * u(I1, I2) + u(I1, I2 - 1));
        }
        else if (option == cIndexing)
        {
            // IndexasCarrays
            const double *u_p = u.getDataPointer();
            double *un_p = un.getDataPointer();
#define U(i1, i2) u_p[(i1 - nd1a) + nd1 * (i2 - nd2a)]
#define UN(i1, i2) un_p[(i1 - nd1a) + nd1 * (i2 - nd2a)]
#pragma omp parallel default(shared) num_threads(thread)
            {
#pragma omp for private(i2, i1)
                for (i2 = n2a; i2 <= n2b; i2++)
                    for (i1 = n1a; i1 <= n1b; i1++)
                    {
                        UN(i1, i2) = U(i1, i2) + rx * (U(i1 + 1, i2) - 2. * U(i1, i2) + U(i1 - 1, i2)) + ry * (U(i1, i2 + 1) - 2. * U(i1, i2) + U(i1, i2 - 1));
                    }
            }
        }
        else if (option == fortranRoutine)
        {
            // callafortranroutine
            // Note:passfirstelementofunanduarraystoFortran(callbyreference)
            if (true)
            {
                const double *u_p = u.getDataPointer();
                double *un_p = un.getDataPointer();
                heat2dUpdate(n1a, n1b, n2a, n2b,
                             u.getBase(0), u.getBound(0), u.getBase(1), u.getBound(1), // passarraydimensions
                             *un_p, *u_p, rx, ry, thread);
            }
            else
            {
                // thiswillalsowork--passfirstelementofunandu
                heat2dUpdate(n1a, n1b, n2a, n2b,
                             u.getBase(0), u.getBound(0), u.getBase(1), u.getBound(1), // passarraydimensions
                             un(nd1a, nd2a), u(nd1a, nd2a), rx, ry, thread);
            }
        }
        else
        {
            printf("ERROR:unknownoption=%d\n", option);
            abort();
        }

        //---boundaryconditions--
        for (int axis = 0; axis < numberOfDimensions; axis++)
            for (int side = 0; side <= 1; side++)
            {
                int is = 1 - 2 * side; // is=1onleft,is=-1onright
                if (boundaryCondition(side, axis) == dirichlet)
                {
                    if (axis == 0)
                    { // leftorrightside
                        i1 = gridIndexRange(side, axis);
                        int i1g = i1 - is; // indexofghostpoint
#pragma omp parallel default(shared) num_threads(thread)
                        {
#pragma omp for private(i2)
                            for (i2 = nd2a; i2 <= nd2b; i2++)
                            {
                                un(i1, i2) = UTRUE(x(i1, i2, 0), x(i1, i2, 1), t + dt);
                                un(i1g, i2) = 3. * un(i1, i2) - 3. * u(i1 + is, i2) + un(i1 + 2 * is, i2); // extrapghost
                            }
                        }
                    }
                    else
                    { // bottomortop
                        i2 = gridIndexRange(side, axis);
                        int i2g = i2 - is; // indexofghostpoint
#pragma omp parallel default(shared) num_threads(thread)
                        {
#pragma omp for private(i1)
                            for (i1 = nd1a; i1 <= nd1b; i1++)
                            {
                                un(i1, i2) = UTRUE(x(i1, i2, 0), x(i1, i2, 1), t + dt);
                                un(i1, i2g) = 3. * un(i1, i2) - 3. * u(i1, i2 + is) + un(i1, i2 + 2 * is); // extrapghost
                            }
                        }
                    }
                }
                else
                {
                    printf("ERROR:unknownboundaryCondition=%d\n", boundaryCondition(side, axis));
                    abort();
                }

            } // endforaxis

        cur = next;
    }
    Real cpuTimeStep = getCPU() - cpu1;

    //---computeerrors--
    RealArray &uc = ua[cur];
    RealArray err(Rx, Ry);

    Real maxErr = 0., maxNorm = 0.;
#pragma omp parallel default(shared) num_threads(thread)
    {
#pragma omp for private(i2, i1) reduction(max : maxErr) reduction(max : maxNorm)
        for (i2 = n2a; i2 <= n2b; i2++)
            for (i1 = n1a; i1 <= n1b; i1++)
            {
                err(i1, i2) = fabs(uc(i1, i2) - UTRUE(x(i1, i2, 0), x(i1, i2, 1), tFinal));
                maxErr = max(err(i1, i2), maxErr);
                maxNorm = max(uc(i1, i2), maxNorm);
            }
    }
    maxErr /= max(maxNorm, REAL_MIN); // relativeerror

    printf("option=%s:numSteps=%dnx=%dmaxNorm=%8.2emaxRelErr=%8.2ecpuTimeStep=%9.2e(s)\n",
           optionName.c_str(), numSteps, nx, maxNorm, maxErr, cpuTimeStep);

    if (nx <= 10)
    {
        uc.display("ua[cur]");
        err.display("err");
    }

    //---OPTIONALLYwriteamatlabfileforplottinginmatlab--
    if (saveMatlab)
    {
        FILE *matlabFile = fopen(matlabFileName.c_str(), "w");
        fprintf(matlabFile, "%%Filewrittenbyheat2d.C\n");
        fprintf(matlabFile, "xa=%g;xb=%g;ya=%g;yb=%g;kappa=%g;t=%g;maxErr=%10.3e;cpuTimeStep=%10.3e;\n",
                xa, xb, ya, yb, kappa, tFinal, maxErr, cpuTimeStep);

        fprintf(matlabFile, "n1a=%d;n1b=%d;nd1a=%d;nd1b=%d;\n", n1a, n1b, nd1a, nd1b);
        fprintf(matlabFile, "n2a=%d;n2b=%d;nd2a=%d;nd2b=%d;\n", n2a, n2b, nd2a, nd2b);
        fprintf(matlabFile, "dx(1)=%14.6e;dx(2)=%14.6e;numGhost=%d;\n", dx[0], dx[1], numGhost);
        fprintf(matlabFile, "option=%d;optionName=\'%s\';\n", option, optionName.c_str());

        if (saveMatlab > 1)
        {
            writeMatlabArray(matlabFile, x, "x", 2, dimension);
            writeMatlabArray(matlabFile, ua[cur], "u", 1, dimension);
            writeMatlabArray(matlabFile, err, "err", 1, dimension);
        }
        fclose(matlabFile);
        printf("Wrotefile[%s]\n", matlabFileName.c_str());
    }
    return 0;
}