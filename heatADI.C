//==============================================================================
// Solvethe2Dheatequationusing anADIscheme
//
//==============================================================================

#include "A++.h"

typedef double Real;
typedef doubleSerialArray RealArray;
typedef intSerialArray IntegerArray;

#include <float.h>
#include <limits.h>
#define REAL_EPSILONDBL_EPSILON
#define REAL_MINDBL_MIN

// getCPU():Returnthecurrentwall-clocktimeinseconds
#include "getCPU.h"

// include commandstoparsecommandlinearguments
#include "parseCommand.h"

// Tridiagonalfactorandsolve:
#include "tridiagonal.h"

// functiontowriteanarraytoamatlabreabablefile:
#include "writeMatlabArray.h"

// define macrosforarrayreferences
#define UN(i1, i2) un_p[(i1) + Ngx * (i2)]
#define U(i1, i2) u_p[(i1) + Ngx * (i2)]

#define AX(i, j) Ax_p[(i) + 3 * (j)]
#define AY(i, j) Ay_p[(i) + 3 * (j)]
#define RHSX(i1) rhsx_p[(i1)]
#define RHSY(i1) rhsy_p[(i1)]

#define RHS(i) rhs_p[(i)]

//=====================================================================================
// FormthetridiagonalmatrixfortheAFSscheme
//=====================================================================================
int formTridiagonalMatrix(RealArray &Ax, const Real kappa, const Real dt, const Real dx)
{
    const int n1a = Ax.getBase(1), n1b = Ax.getBound(1);

    for (int i1 = n1a + 1; i1 <= n1b - 1; i1++)
    {
        Ax(0, i1) = -(.5 * kappa * dt) * (1. / (dx * dx));             // lowerdiagonal
        Ax(1, i1) = 1. - (.5 * kappa * dt) * (-2. * (1. / (dx * dx))); // diagonal
        Ax(2, i1) = -(.5 * kappa * dt) * (1. / (dx * dx));             // upperdiagonal
    }
    Ax(0, n1a) = 0.;
    Ax(1, n1a) = 1.;
    Ax(2, n1a) = 0.; // DirichletBC
    Ax(0, n1b) = 0.;
    Ax(1, n1b) = 1.;
    Ax(2, n1b) = 0.; // DirichletBC

    return 0;
}

int main(int argc, char *argv[])
{
    printf("Usage:heatADI-nx=<i>-tFinal=<f>-saveMatlab=[0|1|2]-matlabFileName=<s>\n");

    const Real pi = M_PI; // 4.*atan2(1.,1.);

    ios::sync_with_stdio();    // SynchronizeC++andCI/Osubsystems
    Index::setBoundsCheck(on); // TurnonA++arrayboundschecking

    int debug = 0;                           // setto1fordebugging
    int plot = 0;                            // set11forplotting
    Real kappa = .05;                        // coefficientofdiffusion
    Real xa = 0., xb = 1., ya = 0., yb = 1.; // spaceintervalinterval
    Real tFinal = .5;                        // finaltime

    int nx = 100, ny = nx;

    int saveMatlab = 0; // 1=saveamatlabfile,2=savesolutiontoo
    string matlabFileName = "heatADI.m";

    const int dirichlet = 1, neumann = 2;
    const int numberOfDimensions = 2;
    IntegerArray boundaryCondition(2, numberOfDimensions);
    boundaryCondition(0, 0) = dirichlet; // left
    boundaryCondition(1, 0) = dirichlet; // right
    boundaryCondition(0, 1) = dirichlet; // bottom
    boundaryCondition(1, 1) = dirichlet; // top

    string line;
    for (int i = 1; i < argc; i++)
    {
        line = argv[i];
        printf("Input:argv[%d]=[%s]\n", i, line.c_str());

        if (parseCommand(line, "-nx=", nx))
        {
            ny = nx;
        }
        else if (parseCommand(line, "-debug=", debug))
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

    // exactsolutionfunction:
    const Real kx = 2., ky = 3;
    const Real kxp = kx * pi;
    const Real kyp = ky * pi;
#define UTRUE(x, y, t) sin(kxp *(x)) * sin(kyp *(y)) * exp(-kappa *(kxp * kxp + kyp * kyp) * (t))

    //---Setupthegrid--
    // xaxb
    //+--+--+--+--....--+--+
    // nx

    // NoghostneededforDirichletBCs
    const int numGhost = 0;

    const int n1a = numGhost;
    const int n2a = numGhost;
    const int n1b = n1a + nx;
    const int n2b = n2a + ny;
    const int nd1a = n1a - numGhost;
    const int nd2a = n2a - numGhost;
    const int nd1b = n1b + numGhost;
    const int nd2b = n2b + numGhost;
    const int nd1 = nd1b - nd1a + 1;
    const int nd2 = nd2b - nd2a + 1;

    IntegerArray gridIndexRange(2, numberOfDimensions);
    IntegerArray dimension(2, numberOfDimensions);
    gridIndexRange(0, 0) = n1a;
    gridIndexRange(1, 0) = n1b;
    gridIndexRange(0, 1) = n2a;
    gridIndexRange(1, 1) = n2b;
    dimension(0, 0) = nd1a;
    dimension(1, 0) = nd1b;
    dimension(0, 1) = nd2a;
    dimension(1, 1) = nd2b;

    //--formthe2Dgridpoints--
    Realdx[2];
    dx[0] = (xb - xa) / nx;
    dx[1] = (yb - ya) / ny;

    Range Rx(nd1a, nd1b), Ry(nd2a, nd2b);
    RealArray x(Rx, Ry, 2);
    int i1, i2;
    for (i2 = nd2a; i2 <= nd2b; i2++)
        for (i1 = nd1a; i1 <= nd1b; i1++)
        {
            x(i1, i2, 0) = xa + (i1 - n1a) * dx[0];
            x(i1, i2, 1) = ya + (i2 - n2a) * dx[1];
        }

    // allocatespaceforthesolution
    RealArray un(Rx, Ry); // holdsU_i^n
    RealArray u(Rx, Ry);  // tempspace

    Real dt = min(dx[0], dx[1]); // timestep(adjustedbelow)
    int Nt = ceil(tFinal / dt);  // numberoftime-steps
    dt = tFinal / Nt;            // adjustdttoreachtFinalexactly

    printf("-----2DHeatEquation:ADIscheme------\n");
    printf("saveMatlab=%d,matlabFileName=%s\n", saveMatlab, matlabFileName.c_str());
    printf("kappa=%.3g,nx=%d,ny=%d,tFinal=%6.2f,kx=%g,ky=%g\n", kappa, nx, ny, tFinal, kx, ky);

    Real t = 0., th, tn;
    Index I1 = Rx, I2 = Ry; // allpoints

    un(I1, I2) = UTRUE(x(I1, I2, 0), x(I1, I2, 1), t); // initialconditions

    //-----Formtheimplicitmatrices-----
    // Note:thematrixissymmetricsowecouldjuststore2entries(*BUT*factorizationneedsextrastorageanyway)
    // Keepallthreeincaseweaddadvection

    const intNgx = nd1, Ngy = nd2; // sizeoftridiagonalsystems
    RealArray Ax(3, Ngx), Ay(3, Ngy);
    RealArray rhsx(Ngx), rhsy(Ngy);

    Real *rhsx_p = rhsx.getDataPointer();
    Real *rhsy_p = rhsy.getDataPointer();

    formTridiagonalMatrix(Ax, kappa, dt, dx[0]);
    factorTridiagonalMatrix(Ax);

    formTridiagonalMatrix(Ay, kappa, dt, dx[1]);
    factorTridiagonalMatrix(Ay);

    const Real rxBy2 = .5 * kappa * dt / (dx[0] * dx[0]);
    const Real ryBy2 = .5 * kappa * dt / (dx[1] * dx[1]);

    // getpointersforarrayreferencemacros:
    Real *un_p = un.getDataPointer();
    Real *u_p = u.getDataPointer();

    //============Starttime-steppingloop==========
    Real cpuTriSolves = 0; // cpufortridiagonalsolves
    Real cpu0 = getCPU();
    for (int n = 0; n < Nt; n++) // time-steploop
    {
        t = n * dt;       // currenttime
        th = t + .5 * dt; // halftime
        tn = t + dt;      // newtime

        // StageI:
        //[I+(.5*kappa*dt)*(D+xD-x)]U^*_ij=[I+(.5*kappa*dt)*(D+yD-y)]U^n_ij
        // Note:saveU^*inU
        for (i2 = n2a + 1; i2 <= n2b - 1; i2++)
            for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
            {
                U(i1, i2) = UN(i1, i2) + ryBy2 * (UN(i1, i2 + 1) - 2. * UN(i1, i2) + UN(i1, i2 - 1)); // rhs
            }

        // BoundaryconditionsW^*=g(x,t+dt/2)
        for (int side = 0; side <= 1; side++)
        {
            i1 = side == 0 ? n1a : n1b;
            for (i2 = n2a; i2 <= n2b; i2++)
                U(i1, i2) = UTRUE(x(i1, i2, 0), x(i1, i2, 1), th); // left/right
            i2 = side == 0 ? n2a : n2b;
            for (i1 = n2a; i1 <= n1b; i1++)
                U(i1, i2) = UTRUE(x(i1, i2, 0), x(i1, i2, 1), th); // bottom/top
        }

        //---tridiagonalsolvesinx-direction--
        Real cpu1 = getCPU();
        for (i2 = n2a + 1; i2 <= n2b - 1; i2++) // excludetopanbottomboundaries
        {
            for (i1 = n1a; i1 <= n1b; i1++)
                RHSX(i1) = U(i1, i2);

            solveTridiagonal(Ax, rhsx);

            for (i1 = n1a; i1 <= n1b; i1++)
                U(i1, i2) = RHSX(i1); // saveU^*inU
        }
        cpuTriSolves += getCPU() - cpu1;

        //-----StageII
        //[I+(.5*kappa*dt)*(D+yD-y)]U^{n+1}_ij=[I+(.5*kappa*dt)*(D+xD-x)]U^*_ij
        for (i2 = n2a + 1; i2 <= n2b - 1; i2++)
            for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
            {
                UN(i1, i2) = U(i1, i2) + rxBy2 * (U(i1 + 1, i2) - 2. * U(i1, i2) + U(i1 - 1, i2)); // rhs
            }

        // BoundaryconditionsW^*=g(x,t+dt/2)
        for (int side = 0; side <= 1; side++)
        {
            i1 = side == 0 ? n1a : n1b;
            for (i2 = n2a; i2 <= n2b; i2++)
                UN(i1, i2) = UTRUE(x(i1, i2, 0), x(i1, i2, 1), tn); // left/right
            i2 = side == 0 ? n2a : n2b;
            for (i1 = n2a; i1 <= n1b; i1++)
                UN(i1, i2) = UTRUE(x(i1, i2, 0), x(i1, i2, 1), tn); // bottom/top
        }

        //---tridiagonalsolvesiny-direction--
        cpu1 = getCPU();
        for (i1 = n1a + 1; i1 <= n1b - 1; i1++) // excludeleftandrightboundaries
        {
            for (i2 = n2a; i2 <= n2b; i2++)
                RHSY(i2) = UN(i1, i2);

            solveTridiagonal(Ay, rhsy);

            for (i2 = n2a; i2 <= n2b; i2++)
                UN(i1, i2) = RHSY(i2);
        }
        cpuTriSolves += getCPU() - cpu1;
    }
    //============ENDtime-steppingloop==========
    Real cpuTimeStep = getCPU() - cpu0;

    t = tn; // lasttime
    if (fabs(t - tFinal) > 1.e-12 * tFinal)
    {
        printf("...doneERROR:t=%12.4e,tFinal=%12.4e,t-tFinal=%9.2e\n", t, tFinal, t - tFinal);
    }

    //---computeerrors--
    RealArray err(Ngx, Ngy);

    Real maxErr = 0., maxNorm = 0.;
    for (i2 = n2a; i2 <= n2b; i2++)
        for (i1 = n1a; i1 <= n1b; i1++)
        {
            err(i1, i2) = fabs(un(i1, i2) - UTRUE(x(i1, i2, 0), x(i1, i2, 1), tFinal));
            maxErr = max(err(i1, i2), maxErr);
            maxNorm = max(un(i1, i2), maxNorm);
        }
    maxErr /= max(maxNorm, REAL_MIN); // relativeerror

    printf("ADI:nx=%3dny=%3dNt=%3d,maxNorm=%8.2emaxRelErr=%8.2ecpu(s):total=%9.2e,triSolves=%9.2e\n",
           nx, ny, Nt, maxNorm, maxErr, cpuTimeStep, cpuTriSolves);

    //---OPTIONALLYwriteamatlabfileforplottinginmatlab--
    if (saveMatlab)
    {
        FILE *matlabFile = fopen(matlabFileName.c_str(), "w");
        fprintf(matlabFile, "%%FilewrittenbyheatADI.C\n");
        fprintf(matlabFile, "xa=%g;xb=%g;ya=%g;yb=%g;kappa=%g;t=%g;maxErr=%10.3e;cpuTimeStep=%10.3e;cpuTriSolves=%10.3e;\n",
                xa, xb, ya, yb, kappa, tFinal, maxErr, cpuTimeStep, cpuTriSolves);

        fprintf(matlabFile, "n1a=%d;n1b=%d;nd1a=%d;nd1b=%d;\n", n1a, n1b, nd1a, nd1b);
        fprintf(matlabFile, "n2a=%d;n2b=%d;nd2a=%d;nd2b=%d;\n", n2a, n2b, nd2a, nd2b);
        fprintf(matlabFile, "dx(1)=%14.6e;dx(2)=%14.6e;numGhost=%d;\n", dx[0], dx[1], numGhost);

        if (saveMatlab > 1)
        {
            writeMatlabArray(matlabFile, x, "x", 2, dimension);
            writeMatlabArray(matlabFile, un, "u", 1, dimension);
            writeMatlabArray(matlabFile, err, "err", 1, dimension);
        }
        fclose(matlabFile);
        printf("Wrotefile[%s]\n", matlabFileName.c_str());
    }

    return 0;
}