//
// Solvetheheatequationinone-dimension
//

// RUN USING:

/// home/henshw/software/mpich-3.3.1-install/bin/mpiexec -n 2 /home/mucelj/Math6890/HW8/heat1d 200 output/temp

#include "getLocalIndexBounds.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <mpi.h>
#include <iostream>

// define anewtype"Real"whichisequivalenttoa"double "
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
    MPI_Init(&argc, &argv);

    int num_processors;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

    string debug_file_name;
    debug_file_name = "debug/" + to_string(rank) + ".txt";
    FILE *debug_file = fopen(debug_file_name.c_str(), "w");
    fprintf(debug_file, "Usage:heat1d[Nx][matlabFileName.m]\n"
                        "Nx=numberofgridcells.\n"
                        "matlabFileName.m:saveresultstothisfile.\n");
#define TRIG_DD1
#define TRIG_NN2
#define POLY_DD3
#define POLY_NN4
//=====Choosethesolutionhereorcompilewith-DSOLUTION=[1|2|3|4]=====
#ifndef SOLUTION
// #define SOLUTION TRIG_DD
// #define SOLUTION TRIG_NN
#define SOLUTION POLY_DD
// #define SOLUTION POLY_NN
#endif

    const Real pi = M_PI;

    int debug = 0; // setto1fordebuginfo
    Real xa = 0., xb = 1.;
    Real kappa = .1;
    Real tFinal = .0001;
    Real cfl = .9; // time-stepsafetyfactor

    int Nx = 10; // default
    string matlabFileName = "heat1d.m";

    if (argc >= 2) // readanycommandlinearguments
    {
        Nx = atoi(argv[1]);
        fprintf(debug_file, "SettingNx=%d\n", Nx);
        if (argc >= 3)
        {
            matlabFileName = argv[2];
            fprintf(debug_file, "SettingmatlabFileName=[%s]\n", matlabFileName.c_str());
        }
    }

    bool is_blocked = atoi(argv[3]);

    int nx_l;
    int n1a_l;
    int n1b_l;

    getLocalIndexBounds(rank, num_processors, Nx, nx_l, n1a_l, n1b_l);

    //=============Gridandindexing==============
    // xaxb
    // G---X---+---+---+---+--...---+---X---G
    // Nx
    // n1an1b
    // nd1and1b
    // Cindex:3...

    Real dx = (xb - xa) / Nx;
    const int numGhost = 1;
    const int n1a = n1a_l;
    const int n1b = n1b_l;
    const int nd1a = n1a - numGhost;
    const int nd1b = n1b + numGhost;
    const int nd1 = nd1b - nd1a + 1; // totalnumberofgridpoints;

    // Createanarrayofgridpoints:
    Real *x_p = new Real[nd1];
#define x(i) x_p[i - nd1a]

    for (int i = nd1a; i <= nd1b; i++)
        x(i) = xa + (i)*dx;

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

    // #if SOLUTION == TRIG_DD
    //     // TruesolutionfordirichletBC’s
    //     boundaryCondition(0, 0) = dirichlet;
    //     boundaryCondition(1, 0) = dirichlet;

    //     const char solutionName[] = "trueDD";

    // #define UTRUE(x, t) sin(kxPi *(x)) * exp(-kappaPiSq *(t))
    // #define UTRUEX(x, t) kxPi *cos(kxPi *(x)) * exp(-kappaPiSq *(t))
    // #define FORCE(x, t) (0.)

    // #elif SOLUTION == TRIG_NN

    //     // TruesolutionforNeumannBC’s
    //     boundaryCondition(0, 0) = neumann;
    //     boundaryCondition(1, 0) = neumann;
    //     const char solutionName[] = "trueNN";

    // #define UTRUE(x, t) cos(kxPi *(x)) * exp(-kappaPiSq *(t))
    // #define UTRUEX(x, t) -kxPi *sin(kxPi *(x)) * exp(-kappaPiSq *(t))
    // #define FORCE(x, t) (0.)

#if (SOLUTION == POLY_DD) || (SOLUTION == POLY_NN)

// polynomialmanufacturedsolution
#if SOLUTION == POLY_DD
    const char solutionName[] = "polyDD";
    boundaryCondition(0, 0) = dirichlet;
    boundaryCondition(1, 0) = dirichlet;
#else
    const charsolutionName[] = "polyNN";
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
    printf("ERROR:unknownsolution");
    abort();
#endif

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
    int numSteps = ceil(tFinal / dt);
    dt = tFinal / numSteps; // adjustdttoreachthefinaltime
    const Real rx = kappa * dt / dx2;

    fprintf(debug_file, "-------------------Solvetheheatequationin1Dsolution=%s---------------------\n",
            solutionName);
    fprintf(debug_file, "numGhost=%d,n1a=%d,n1b=%d,nd1a=%d,nd1b=%d\n", numGhost, n1a, n1b, nd1a, nd1b);
    fprintf(debug_file, "numSteps=%d,Nx=%d,kappa=%g,tFinal=%g,boundaryCondition(0,0)=%d,boundaryCondition(1,0)=%d\n",
            numSteps, Nx, kappa, tFinal, boundaryCondition(0, 0), boundaryCondition(1, 0));

    //----------TIME-STEPPINGLOOP--------
    Real cpu0 = getCPU();

    MPI_Request send_left;
    MPI_Request receive_left;

    MPI_Request send_right;
    MPI_Request receive_right;

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
        // for (int side = 0; side <= 1; side++)
        // {
        //     const int i = side == 0 ? n1a : n1b; // boundaryindex
        //     const int is = 1 - 2 * side;         // is=1onleft,-1onright
        //     if (boundaryCondition(side, 0) == dirichlet)
        //     {
        //         un(i) = UTRUE(x(i), t + dt);
        //         un(i - is) = 3. * un(i) - 3. * un(i + is) + un(i + 2 * is); // extrapolateghost
        //     }
        //     else
        //     {
        //         // NeumannBC
        //         un(i - is) = un(i + is) - 2. * is * dx * UTRUEX(x(i), t + dt);
        //     }
        // }

        // ---- new boundary contitions ----
        if (rank == 0)
        {
            int side = 0;
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
            if (num_processors > 1)
                if (is_blocked)
                {
                    // Send to rank 1
                    MPI_Send(&un(n1b), 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
                    if (debug > 0)
                    {
                        fprintf(debug_file, "Sending un(%d) = %f to rank %d\n", n1b, un(n1b), 1);
                        fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1b), t + dt));
                    }

                    // Receive from rank 1
                    MPI_Recv(&un(n1b + 1), 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if (debug > 0)
                    {
                        fprintf(debug_file, "Receiving un(%d) = %f from rank %d\n", n1b + 1, un(n1b + 1), 1);
                        fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1b + 1), t + dt));
                    }
                }
                else
                {
                    // Send to rank 1
                    MPI_Isend(&un(n1b), 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &send_right);
                    if (debug > 0)
                    {
                        fprintf(debug_file, "Sending un(%d) = %f to rank %d\n", n1b, un(n1b), 1);
                        fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1b), t + dt));
                    }

                    // Receive from rank 1
                    MPI_Irecv(&un(n1b + 1), 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &receive_right);
                    if (debug > 0)
                    {
                        fprintf(debug_file, "Receiving un(%d) = %f from rank %d\n", n1b + 1, un(n1b + 1), 1);
                        fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1b + 1), t + dt));
                    }
                }
        }
        if (rank == num_processors - 1)
        {
            int side = 1;
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
            if (num_processors > 1)
                if (is_blocked)
                {
                    // Receive from rank num_p-2
                    MPI_Recv(&un(n1a - 1), 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if (debug > 0)
                    {
                        fprintf(debug_file, "Receiving un(%d) = %f from rank %d\n", n1a - 1, un(n1a - 1), rank - 1);
                        fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1a - 1), t + dt));
                    }

                    // Send to rank num_p-2
                    MPI_Send(&un(n1a), 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
                    if (debug > 0)
                    {
                        fprintf(debug_file, "Sending un(%d) = %f to rank %d\n", n1a, un(n1a), rank - 1);
                        fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1a), t + dt));
                    }
                }
                else
                {
                    // Receive from rank num_p-2
                    MPI_Irecv(&un(n1a - 1), 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &receive_left);
                    if (debug > 0)
                    {
                        fprintf(debug_file, "Receiving un(%d) = %f from rank %d\n", n1a - 1, un(n1a - 1), rank - 1);
                        fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1a - 1), t + dt));
                    }

                    // Send to rank num_p-2
                    MPI_Isend(&un(n1a), 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &send_left);
                    if (debug > 0)
                    {
                        fprintf(debug_file, "Sending un(%d) = %f to rank %d\n", n1a, un(n1a), rank - 1);
                        fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1a), t + dt));
                    }
                }
        }
        if (rank != 0 && rank != num_processors - 1)
        {
            if (is_blocked)
            {
                // Send to rank "right"
                MPI_Send(&un(n1b), 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                if (debug > 0)
                {
                    fprintf(debug_file, "Sending un(%d) = %f to rank %d\n", n1b, un(n1b), rank + 1);
                    fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1b), t + dt));
                }

                // Receive from rank "left"
                MPI_Recv(&un(n1a - 1), 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (debug > 0)
                {
                    fprintf(debug_file, "Receiving un(%d) = %f from rank %d\n", n1a - 1, un(n1a - 1), rank - 1);
                    fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1a - 1), t + dt));
                }

                // Send to rank "left"

                MPI_Send(&un(n1a), 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
                if (debug > 0)
                {
                    fprintf(debug_file, "Sending un(%d) = %f to rank %d\n", n1a, un(n1a), rank - 1);
                    fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1a), t + dt));
                }

                // Receive from rank "right"
                MPI_Recv(&un(n1b + 1), 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (debug > 0)
                {
                    fprintf(debug_file, "Receiving un(%d) = %f from rank %d\n", n1b + 1, un(n1b + 1), rank + 1);
                    fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1b + 1), t + dt));
                }
            }
            else
            {
                // Send to rank "right"
                MPI_Isend(&un(n1b), 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &send_right);
                if (debug > 0)
                {
                    fprintf(debug_file, "Sending un(%d) = %f to rank %d\n", n1b, un(n1b), rank + 1);
                    fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1b), t + dt));
                }

                // Receive from rank "left"
                MPI_Irecv(&un(n1a - 1), 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &receive_left);
                if (debug > 0)
                {
                    fprintf(debug_file, "Receiving un(%d) = %f from rank %d\n", n1a - 1, un(n1a - 1), rank - 1);
                    fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1a - 1), t + dt));
                }

                // Send to rank "left"
                MPI_Isend(&un(n1a), 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &send_left);
                if (debug > 0)
                {
                    fprintf(debug_file, "Sending un(%d) = %f to rank %d\n", n1a, un(n1a), rank - 1);
                    fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1a), t + dt));
                }

                // Receive from rank "right"
                MPI_Irecv(&un(n1b + 1), 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &receive_right);
                if (debug > 0)
                {
                    fprintf(debug_file, "Receiving un(%d) = %f from rank %d\n", n1b + 1, un(n1b + 1), rank + 1);
                    fprintf(debug_file, "Actuall value   = %f \n", UTRUE(x(n1b + 1), t + dt));
                }
            }
        }

        if (!is_blocked)
        {
            if (rank != num_processors - 1)
            {
                MPI_Wait(&send_right, MPI_STATUS_IGNORE);
                MPI_Wait(&receive_right, MPI_STATUS_IGNORE);
            }
            if (rank != 0)
            {
                MPI_Wait(&send_left, MPI_STATUS_IGNORE);
                MPI_Wait(&receive_left, MPI_STATUS_IGNORE);
            }
        }

        if (debug > 1)
        {
            fprintf(debug_file, "step%d:AfterupdateinteriorandrealBCs\nu=[", n + 1);
            for (int i = nd1a; i <= nd1b; i++)
                printf("%12.4e,", un(i));
            fprintf(debug_file, "]\n");
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
            fprintf(debug_file, "step=%d,t=%9.3e,maxErr=%9.2e\n", n + 1, t + dt, maxErr);
        }

    } // endtime-steppingloop
    // if (!is_blocked)
    // {
    //     if (rank != num_processors - 1)
    //     {
    //         MPI_Wait(&send_right, MPI_STATUS_IGNORE);
    //         MPI_Wait(&receive_right, MPI_STATUS_IGNORE);
    //     }
    //     if (rank != 0)
    //     {
    //         MPI_Wait(&send_left, MPI_STATUS_IGNORE);
    //         MPI_Wait(&receive_left, MPI_STATUS_IGNORE);
    //     }
    // }

    Real cpuTimeStep = getCPU() - cpu0;
    Real *cpuTimes;
    Real max_cpu = cpuTimeStep;
    int tag = 20;
    if (rank == 0)
    {
        cpuTimes = new Real[num_processors];
        cpuTimes[0] = cpuTimeStep;
        for (int i = 0; i < num_processors - 1; i++)
        {

            MPI_Recv(&max_cpu, 1, MPI_DOUBLE, i + 1, tag + i + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            fprintf(debug_file, "RECEIVED MESSAGE FROM RANK(%d) with tag(%d): %f \n", i + 1, tag + i + 1, max_cpu);
            cpuTimes[i + 1] = max_cpu;
        }
    }
    else
    {
        MPI_Send(&cpuTimeStep, 1, MPI_DOUBLE, 0, tag + rank, MPI_COMM_WORLD);
        fprintf(debug_file, "SENT MESSAGE FROM RANK(%d) with tag(%d): %f\n\ns", rank, tag + rank, cpuTimeStep);
    }
    if (rank == 0)
    {
        for (int i = 0; i < num_processors; i++)
        {
            if (max_cpu < cpuTimes[i])
            {
                max_cpu = cpuTimes[i];
            }
        }
    }
    Real maxErr;
    maxErr = 0.;

    //----checktheerror----
    t += dt; // tFinal;
    if (fabs(t - tFinal) > 1e-3 * dt / tFinal)
    {
        printf("ERROR:AFTERTIME_STEPPING:t=%16.8eISNOTEQUALtotFinal=%16.8e\n", t, tFinal);
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
    fprintf(debug_file, "numSteps=%4d,Nx=%3d,maxErr=%9.2e,cpu=%9.2e(s)\n", numSteps, Nx, maxErr, cpuTimeStep);

    if (rank == 0)
    {
        //-- -Writeafileforplottinginmatlab--
        FILE *matlabFile = fopen(matlabFileName.c_str(), "w");
        fprintf(matlabFile, "%%Filewrittenbyheat1d.C\n");
        fprintf(matlabFile, "xa=%g;xb=%g;kappa=%g;t=%g;maxErr=%10.3e;cpuTimeStep=%10.3e;\n", xa, xb, kappa, tFinal, maxErr, cpuTimeStep);

        fprintf(debug_file, "max_cpu = %f\n", max_cpu);

        fprintf(matlabFile, "Nx=%d;dx=%14.6e;numGhost=%d;n1a=%d;n1b=%d;nd1a=%d;nd1b=%d;\n", Nx, dx, numGhost, n1a, n1b, nd1a, nd1b);
        fprintf(matlabFile, "solutionName=\'%s\';\n", solutionName);

        writeMatlabVector(matlabFile, x_p, "x", nd1a, nd1b);
        writeMatlabVector(matlabFile, u_p[cur], "u", nd1a, nd1b);
        writeMatlabVector(matlabFile, error_p, "err", nd1a, nd1b);

        // fclose(matlabFile);
        fprintf(debug_file, "Wrotefile%s\n\n", matlabFileName.c_str());
        delete[] error_p;
    }
    Real *maxErrors = new Real[num_processors];
    if (rank == 0)
    {

        maxErrors[0] = maxErr;
        for (int i = 0; i < num_processors - 1; i++)
        {

            MPI_Recv(&maxErr, 1, MPI_DOUBLE, i + 1, tag + i + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            fprintf(debug_file, "RECEIVED MESSAGE FROM RANK(%d) with tag(%d): %f \n", i + 1, tag + i + 1, maxErr);
            maxErrors[i + 1] = maxErr;
        }
    }
    else
    {
        MPI_Send(&maxErr, 1, MPI_DOUBLE, 0, tag + rank, MPI_COMM_WORLD);
        fprintf(debug_file, "SENT MESSAGE FROM RANK(%d) with tag(%d): %f\n\ns", rank, tag + rank, maxErr);
    }
    if (rank == 0)
    {
        for (int i = 0; i < num_processors; i++)
        {
            if (maxErr < maxErrors[i])
            {
                maxErr = maxErrors[i];
            }
        }
        fprintf(debug_file, "Maximum Error: %f\n", maxErr);
    }
    delete[] u_p[0];
    delete[] u_p[1];
    delete[] x_p;

    delete[] boundaryCondition_p;
    MPI_Finalize();
    return 0;
}