//=============================================
// Heat equation in 2D; solved with A++ arrays
//=============================================

/// home/henshw/software/mpich-3.3.1-install/bin/mpiexec -n 2 /home/mucelj/Math6890/HW8/heat1d 200 output/temp

#include "A++.h"
#include <mpi.h>

// typdefing & defining marcos
#define REAL_EPSILON DBL_EPSILON
#define REAL_MIN DBL_MIN

typedef double Real;
typedef doubleSerialArray RealArray;
typedef intSerialArray IntegerArray;

// include necessary libraries
#include "getCPU.h"           //return clock-wall time in secs
#include "parseCommand.h"     //pasrse command line arguments
#include "writeMatlabArray.h" //write array to matlab file
#include <float.h>
#include <limits.h>
#include "getLocalIndexBounds.h"

#define heat2dUpdate heat2dupdate_ // declare fortan rountine as "C" function
extern "C"
{
    void heat2dUpdate(const int &n1a, const int &n1b, const int &n2a, const int &n2b,
                      const int &nd1a, const int &nd1b, const int &nd2a, const int &nd2b,
                      Real &un, const Real &u, const Real &rx, const Real &ry);
}

//================== MAIN FUNCTION =====================
int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);

    int num_processors;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

    string debug_file_name;
    debug_file_name = "debug/" + to_string(rank) + ".txt";
    FILE *debug_file = fopen(debug_file_name.c_str(), "w");

    const Real pi = 4. * atan2(1., 1.);

    ios::sync_with_stdio();    // sync C++ and C I/O subsystems
    Index::setBoundsCheck(on); // turn on A++ array bound checking

    fprintf(debug_file, "Usage: heat2d -nx = <i> -option=[0|1|2|3] -tFinal=<f> -debug = <i> -saveMatlab=[0|1|2] -matlabFile=<s>\n"
                        "   option : 0 = scalarIndexing, 1 = arrayIndexing, 2=cIndexing, 3=fortranRoutine\n");

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

    // define some variables
    int debug = 0;
    Real kappa = .1;
    Real xa = 0., xb = 1.; // domain is [xa, ab] x [ya, yb]
    Real ya = 0., yb = 1.;

    Real tFinal = .5;
    int nx = 100, ny = nx;

    int saveMatlab = 0; // 1 = save a matlab file, 2 = save solution too
    string matlabFileName = "heat2d.m";

    string line;
    // parse command line arguments
    for (int i = 1; i < argc; i++)
    {
        line = argv[i];
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
    }

    int nx_l;
    int n1a_l;
    int n1b_l;

    getLocalIndexBounds(rank, num_processors, nx, nx_l, n1a_l, n1b_l);
    cout << n1a_l << " " << n1b_l << endl;

    // define some more variables; for grid + indexing
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

    // fill in dimensions & boundary conditions
    IntegerArray gridIndexRange(2, numberOfDimensions);
    IntegerArray dimension(2, numberOfDimensions);
    IntegerArray boundaryCondition(2, numberOfDimensions);

    gridIndexRange(0, 0) = n1a;
    gridIndexRange(1, 0) = n1b; // left, right
    gridIndexRange(0, 1) = n2a;
    gridIndexRange(1, 1) = n2b; // down, up

    dimension(0, 0) = nd1a;
    dimension(1, 0) = nd1b;
    dimension(0, 1) = nd2a;
    dimension(1, 1) = nd2b;

    boundaryCondition(0, 0) = dirichlet; // left
    boundaryCondition(1, 0) = dirichlet; // right
    boundaryCondition(0, 1) = dirichlet; // bottom
    boundaryCondition(1, 1) = dirichlet; // top

    // grid points
    Range Rx(nd1a, nd1b), Ry(nd2a, nd2b);
    RealArray x(Rx, Ry, 2);

    Real dx[2];
    dx[0] = (xb - xa) / nx;
    dx[1] = (yb - ya) / ny; // dy

    int i1, i2;
    for (i2 = nd2a; i2 <= nd2b; i2++)
    {
        for (i1 = nd1a; i1 <= nd1b; i1++)
        {
            x(i1, i2, 0) = xa + (i1 - n1a) * dx[0];
            x(i1, i2, 1) = ya + (i2 - n2a) * dx[1];
        }
    }

    const Real kx = 2., ky = 3;
    const Real kxp = kx * pi;
    const Real kyp = ky * pi;
    // #define UTRUE(x, y, t) sin(kxp *(x)) * sin(kyp *(y)) * exp(-kappa *(kxp * kxp + kyp * kyp) * (t))
    static const Real c0 = .2, c1 = .1, c2 = .3;
    static const Real b0 = 1., b1 = .5, b2 = .25;
    static const Real a0 = 1., a1 = .3, a2 = 0.;
#define UTRUE(x, y, t) (b0 + (x) * (b1 + (x)*b2)) * (c0 + (y) * (c1 + (y)*c2)) * (a0 + (t) * (a1 + (t)*a2))
#define UTRUET(x, y, t) (b0 + (x) * (b1 + (x)*b2)) * (c0 + (y) * (c1 + (y)*c2)) * (a1 + 2. * (t)*a2)
#define UTRUEXX(x, y, t) (2. * b2) * (c0 + (y) * (c1 + (y)*c2)) * (a0 + (t) * (a1 + (t)*a2))
#define UTRUEYY(x, y, t) (b0 + (x) * (b1 + (x)*b2)) * (2. * c2) * (a0 + (t) * (a1 + (t)*a2))
#define FORCE(x, y, t) (UTRUET(x, y, t) - kappa * (UTRUEXX(x, y, t) + UTRUEYY(x, y, t)))
    string optionName = option == scalarIndexing ? "scalarIndexing" : option == arrayIndexing ? "arrayIndexing "
                                                                  : option == cIndexing       ? "cIndexing     "
                                                                  : option == fortranRoutine  ? "fortranRoutine"
                                                                                              : "uknown";

    // store two time levels
    RealArray ua[2];
    ua[0].redim(Rx, Ry);
    ua[0] = 0.;
    ua[1].redim(Rx, Ry);
    ua[1] = 0.;

    // intial conditions
    RealArray &u0 = ua[0];
    Real t = 0.;

    for (i2 = nd2a; i2 <= nd2b; i2++)
    {
        for (i1 = nd1a; i1 <= nd1b; i1++)
        {
            u0(i1, i2) = UTRUE(x(i1, i2, 0), x(i1, i2, 1), t);
        }
    }

    // time-step restriction:
    // forward euler: kappa*dt*(1/dx^2 + 1/dy^2) < cfl*.5
    Real cfl = .9;
    Real dt = cfl * (.5 / kappa) / (1. / (dx[0] * dx[0]) + 1. / (dx[1] * dx[1]));

    int numSteps = ceil(tFinal / dt);
    dt = tFinal / numSteps; // adjust dt to reach final time

    Real rx = kappa * (dt / (dx[0] * dx[0]));
    Real ry = kappa * (dt / (dx[1] * dx[1]));
    int cur = 0;
    int i, n;
    Index I1 = Range(n1a, n1b);
    Index I2 = Range(n2a, n2b); // is this correct ?

    Real cpu1 = getCPU();
    for (n = 0; n < numSteps; n++)
    {
        t = n * dt; // cur time

        int next = (cur + 1) % 2;
        RealArray &u = ua[cur];
        RealArray &un = ua[next];

        if (option == scalarIndexing)
        {
            for (i2 = n2a; i2 <= n2b; i2++)
            {
                for (i1 = n1a_l; i1 <= n1b_l; i1++)
                {
                    un(i1, i2) = dt * FORCE(x(i1, i2, 0), x(i1, i2, 1), t) + u(i1, i2) + rx * (u(i1 + 1, i2) - 2. * u(i1, i2) + u(i1 - 1, i2)) + ry * (u(i1, i2 + 1) - 2. * u(i1, i2) + u(i1, i2 - 1));
                }
            }
        }

        //============= boundary condition =============
        for (int axis = 0; axis <= 1; axis++)
        {
            for (int side = 0; side <= 1; side++)
            {
                int is = 1 - 2 * side; // is = 1 on left, is=-1 on right
                if (boundaryCondition(side, axis) == dirichlet)
                {
                    if (axis == 0)
                    {
                        if ((rank == 0 && side == 0) || (rank == num_processors - 1 && side == 1))
                        {
                            // left or right side
                            i1 = gridIndexRange(side, axis);

                            int i1g = i1 - is; // index of ghost point
                            for (i2 = nd2a; i2 <= nd2b; i2++)
                            {
                                un(i1, i2) = UTRUE(x(i1, i2, 0), x(i1, i2, 1), t + dt);
                                un(i1g, i2) = 3. * un(i1, i2) - 3. * u(i1 + is, i2) + un(i1 + 2 * is, i2); // extrap ghost
                            }
                        }
                        else if (side == 0)
                        {
                            i1 = n1a_l;

                            int i1g = i1 - 1; // index of ghost point

                            Real *buff = new Real[nd2b - nd2a + 1];
                            Real *buff2 = new Real[nd2b - nd2a + 1];
                            for (i2 = nd2a; i2 <= nd2b; i2++)
                            {
                                buff[i2 - nd2a] = un(i1, i2);
                            }
                            MPI_Send(buff, nd2b - nd2a + 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
                            MPI_Recv(buff2, nd2b - nd2a + 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            for (i2 = nd2a + 1; i2 <= nd2b - 1; i2++)
                            {

                                un(i1g, i2) = buff2[i2 - nd2a];
                            }
                        }
                        else if (side == 1)
                        {
                            i1 = n1b_l;

                            int i1g = i1 + 1; // index of ghost point

                            Real *buff = new Real[nd2b - nd2a + 1];
                            Real *buff2 = new Real[nd2b - nd2a + 1];
                            for (i2 = nd2a; i2 <= nd2b; i2++)
                            {
                                buff[i2 - nd2a] = un(i1, i2);
                            }
                            MPI_Send(buff, nd2b - nd2a + 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                            MPI_Recv(buff2, nd2b - nd2a + 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            for (i2 = nd2a + 1; i2 <= nd2b - 1; i2++)
                            {

                                un(i1g, i2) = buff2[i2 - nd2a];
                            }
                        }
                    }
                    else
                    {
                        // bottom or top
                        i2 = gridIndexRange(side, axis);
                        int i2g = i2 - is; // index of ghost point
                        for (i1 = nd1a; i1 <= nd1b; i1++)
                        {
                            un(i1, i2) = UTRUE(x(i1, i2, 0), x(i1, i2, 1), t + dt);
                            un(i1, i2g) = 3. * un(i1, i2) - 3. * u(i1, i2 + is) + un(i1, i2 + 2 * is); // extrap ghost
                        }
                    }
                }
                else
                {
                    abort();
                }
            }
        } // end for axis
        cur = next;
    }
    Real cpuTimeStep = getCPU() - cpu1;

    //============= compute errors ==============
    RealArray &uc = ua[cur];
    RealArray err(Rx, Ry);

    Real maxErr = 0., maxNorm = 0.;
    for (i2 = n2a; i2 <= n2b; i2++)
    {
        for (i1 = n1a_l; i1 <= n1b_l; i1++)
        {
            err(i1, i2) = fabs(uc(i1, i2) - UTRUE(x(i1, i2, 0), x(i1, i2, 1), tFinal));
            maxErr = max(err(i1, i2), maxErr);
            maxNorm = max(uc(i1, i2), maxNorm);
        }
    }
    maxErr /= max(maxNorm, REAL_MIN); // relative error

    if (nx <= 10)
    {
        uc.display("ua[cur]");
        err.display("err");
    }

    MPI_Finalize();
    return 0;
}
