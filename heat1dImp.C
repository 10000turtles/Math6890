//===================================================================================
//
// Solvetheheatequationinone-dimensionwithanIMPLICTMETHOD
// TRAPEZOIDALRULEINTIME
//
//===================================================================================

#include "A++.h"

// Tridiagonalfactorandsolve:
#include "tridiagonal.h"

// define sometypes
typedef double Real;
typedef double SerialArrayRealArray;
typedef int SerialArrayIntegerArray;

#include <string>
using std::max;
using std::string;

// getCPU():Returnthecurrentwall-clocktimeinseconds
#include "getCPU.h"

// include commandstoparsecommandlinearguments
#include "parseCommand.h"

// Functiontosaveavectortoamatlabfile.
#include "writeMatlabVector.h"

// Truesolutionoptions:
static const int trueDD = 1;
static const int trueNN = 2;
static const int poly = 3;

static Real kappa = .1;
static Real kx = 3.;
static Real kxPi = kx * M_PI;
static Real kappaPiSq = kappa * kxPi * kxPi;

// TruesolutionfordirichletBC’s
#define UTRUEDD(x, t) sin(kxPi *(x)) * exp(-kappaPiSq *(t))
#define UTRUEDDx(x, t) kxPi *cos(kxPi *(x)) * exp(-kappaPiSq *(t))
#define FORCEDD(x, t) (0.)

// TruesolutionforNeumannBC’s
#define UTRUENN(x, t) cos(kxPi *(x)) * exp(-kappaPiSq *(t))
#define UTRUENNx(x, t) -kxPi *sin(kxPi *(x)) * exp(-kappaPiSq *(t))
#define FORCENN(x, t) (0.)

// polynomialmanufacturedsolution
static const Real b0 = 1., b1 = .5, b2 = .25;
static const Real a0 = 1., a1 = .3, a2 = -.1;
#define POLY(x, t) (b0 + (x) * (b1 + (x)*b2)) * (a0 + (t) * (a1 + (t)*a2))
#define POLYx(x, t) (b1 + 2. * (x)*b2) * (a0 + (t) * (a1 + (t)*a2))
#define POLYT(x, t) (b0 + (x) * (b1 + (x)*b2)) * (a1 + 2. * (t)*a2)
#define POLYXx(x, t) (2. * b2) * (a0 + (t) * (a1 + (t)*a2))
#define POLYFORCE(x, t) (POLYT(x, t) - kappa * POLYXx(x, t))

//--------------------------------------------------------------------------------------
// Functiontoevaluatethetruesolution
// solutionOption(input):truesolutionoption
// t(input):time
// x(input):gridpoints
// I1(input):evaluateattheseindexvalues
// uTrue(output):uTrue(I1)=truesolution
//--------------------------------------------------------------------------------------
int getTrue(int solutionOption, Real t, RealArray &x, Index &I1, RealArray &uTrue)
{
    if (solutionOption == trueDD)
        uTrue(I1) = UTRUEDD(x(I1), t);
    else if (solutionOption == trueNN)
        uTrue(I1) = UTRUENN(x(I1), t);
    else if (solutionOption == poly)
        uTrue(I1) = POLY(x(I1), t);
    else
    {
        printf("getTrue:unknownsolutonOption=%d\n", solutionOption);
        abort();
    }

    return 0;
}

//--------------------------------------------------------------------------------------
// Functiontoevaluatethex-derivativeofthetruesolution
//--------------------------------------------------------------------------------------
int getTruex(int solutionOption, Real t, RealArray &x, Index &I1, RealArray &uTruex)
{
    if (solutionOption == trueDD)
        uTruex(I1) = UTRUEDDx(x(I1), t);
    else if (solutionOption == trueNN)
        uTruex(I1) = UTRUENNx(x(I1), t);
    else if (solutionOption == poly)
        uTruex(I1) = POLYx(x(I1), t);
    else
    {
        printf("getTrue:unknownsolutonOption=%d\n", solutionOption);
        abort();
    }
    return 0;
}

//--------------------------------------------------------------------------------------
// FunctiontoevaluatethePDEforcing
//--------------------------------------------------------------------------------------
int getForce(int solutionOption, Real t, RealArray &x, Index &I1, RealArray &force)
{
    if (solutionOption == trueDD)
        force(I1) = FORCEDD(x(I1), t);
    elseif(solutionOption == trueNN)
        force(I1) = FORCENN(x(I1), t);
    elseif(solutionOption == poly)
        force(I1) = POLYFORCE(x(I1), t);
    else
    {
        printf("getTrue:unknownsolutonOption=%d\n", solutionOption);
        abort();
    }
    return 0;
}

int main(int argc, char *argv[])
{

    printf("Usage:heat1dImp-Nx=<i>-tFinal=<f>-sol=[true|poly]-bc1=[d|n]-bc2=[d|n]-debug=<i>matlabFileName=<s>\n");

    const Realpi = M_PI;

    const intnumberOfDimensions = 1;
    int Nx = 10;

    Real xa = 0., xb = 1.;
    Real tFinal = 1.;

    // setupboundaryconditionarray
    const intdirichlet = 1, neumann = 2;
    IntegerArrayboundaryCondition(2, numberOfDimensions);
    boundaryCondition(0, 0) = dirichlet; // left
    boundaryCondition(1, 0) = dirichlet; // right

    int debug = 0;
    string matlabFileName = "heat1d.m";

    const int trueSolution = 0, polynomialSolution = 1;
    int sol = trueSolution;
    string solName = "true";

    string line;
    for (int i = 1; i < argc; i++)
    {
        line = argv[i];
        // printf("Input:argv[%d]=[%s]\n",i,line.c_str());

        if (parseCommand(line, "-Nx=", Nx))
        {
        }
        else if (parseCommand(line, "-debug=", debug))
        {
        }
        else if (parseCommand(line, "-tFinal=", tFinal))
        {
        }
        else if (parseCommand(line, "-matlabFileName=", matlabFileName))
        {
        }
        else if (line.substr(0, 4) == "-sol")
        {
            solName = line.substr(5); // solName=character5toend
            if (solName == "true")
                sol = trueSolution;
            else if (solName == "poly")
                sol = polynomialSolution;
            else
            {
                printf("ERROR:Unknown-sol=[%s]\n", solName.c_str());
                abort();
            }
            printf("settingsolName=[%s]\n", solName.c_str());
        }
        else if (line.substr(0, 5) == "-bc1=" ||
                 line.substr(0, 5) == "-bc2=")
        {
            int side = line.substr(0, 5) == "-bc1=" ? 0 : 1;
            string bcName = line.substr(5, 1);
            if (bcName == "d")
            {
                boundaryCondition(side, 0) = dirichlet;
                printf("SETTINGboundaryCondition(%d,0)=dirichlet\n", side);
            }
            elseif(bcName == "n")
            {
                boundaryCondition(side, 0) = neumann;
                printf("SETTINGboundaryCondition(%d,0)=neumann\n", side);
            }
            else
            {
                printf("Uknownbc:line=[%s]\n", line.c_str());
            }
        }
    }
    string bcName[2];
    for (int side = 0; side <= 1; side++)
    {
        if (boundaryCondition(side, 0) == dirichlet)
            bcName[side] = "D";
        else
            bcName[side] = "N";
    }

    // setthesolutionOption:
    intsolutionOption = trueDD;
    if (boundaryCondition(0, 0) == dirichlet && boundaryCondition(1, 0) == dirichlet)
    {
        solutionOption = sol == trueSolution ? trueDD : poly;
    }
    else if (boundaryCondition(0, 0) == neumann && boundaryCondition(1, 0) == neumann)
    {
        solutionOption = sol == trueSolution ? trueNN : poly;
    }
    else
    {
        printf("Unexpectedboundaryconditionswiththetruesolution--FINISHME\n");
        abort();
    }
    const string solutionName = solName + bcName[0] + bcName[1];

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

    int nd1 = nd1b - nd1a + 1; // totalnumberofgridpoints;

    // Gridpoints:
    Range Rx(nd1a, nd1b);
    RealArray x(Rx);

    for (int i = nd1a; i <= nd1b; i++)
        x(i) = xa + (i - n1a) * dx;

    if (debug > 1)
    {
        for (inti = nd1a; i <= nd1b; i++)
            printf("x(%2d)=%12.4e\n", i, x(i));
    }

    RealArrayu[2]; // twoarrayswillbeusedforcurrentandnewtimes
    u[0].redim(Rx);
    u[1].redim(Rx);

// Macrostodefine fortranlikearrays
#define uc(i) u[cur](i)
#define un(i) u[next](i)

    // initialconditions
    Real t = 0.;
    int cur = 0; //"current"solution,indexintou_p[]
    Index I1 = Range(nd1a, nd1b);
    getTrue(solutionOption, t, x, I1, u[cur]);

    if (debug > 0)
    {
        printf("Afterinitialconditions\nu=[");
        for (inti = nd1a; i <= nd1b; i++)
            printf("%10.4e,", uc(i));
        printf("]\n");
    }

    // Choosetime-step
    const Real dx2 = dx * dx;
    Real dt = dx; // dt,adjustedbelow
    const int numSteps = ceil(tFinal / dt);
    dt = tFinal / numSteps; // adjustdttoreachthefinaltime

    //---Tridiagonalmatrix---
    // Ax(0:2,i1):holdsthe3diagonals
    //[Ax(1,0)Ax(2,0)]
    //[Ax(0,1)Ax(1,1)Ax(2,1)]
    //[Ax(0,2)Ax(1,2)Ax(2,2)]
    //
    Range Ix(n1a, n1b); // interiorandboundarypoints
    RealArray Ax(Range(3), Ix);
    RealArray rhsx(Ix);
    Real *rhsx_p = rhsx.getDataPointer();
#define RHS(i) rhsx_p[i - n1a]

    //----Fillthetridiagonalmatrix---
    const Real rx = kappa * dt / dx2;
    for (int i1 = n1a + 1; i1 <= n1b - 1; i1++)
    {
        Ax(0, i1) = -.5 * rx; // lowerdiagonal
        Ax(1, i1) = 1. + rx;  // diagonal
        Ax(2, i1) = -.5 * rx; // upperdiagonal
    }
    for (int side = 0; side <= 1; side++)
    {
        int i1 = side == 0 ? n1a : n1b;
        if (boundaryCondition(side, 0) == dirichlet)
        { // DirichletBC
            Ax(0, i1) = 0.;
            Ax(1, i1) = 1.;
            Ax(2, i1) = 0.;
        }
        else
        { // NeumannBC
            // CombineNeumannBCwithinteriorequationontheboundary
            // toeliminatetheghostpoint
            int is = 1 - 2 * side;
            Ax(1 - is, i1) = 0.;
            Ax(1, i1) = 1. + rx;
            Ax(1 + is, i1) = -rx; // NeumannBC
        }
    }

    // FactormatrixONCE
    factorTridiagonalMatrix(Ax);

    printf("---------ImplicitSolveoftheHeatEquationin1D,solutionOption=%d,solutionName=%s-----------\n",
           solutionOption, solutionName.c_str());
    printf("numSteps=%d,Nx=%d,debug=%d,kappa=%g,tFinal=%g,\n"
           "boundaryCondition(0,0)=%s,boundaryCondition(1,0)=%s\n",
           numSteps, Nx, debug, kappa, tFinal, bcName[0].c_str(), bcName[1].c_str());

    RealArray uTrue(Rx); // storeuTruehere
    RealArray fn[2];     // saveforcingattwotimelevels
    fn[0].redim(Rx);
    fn[1].redim(Rx);

    getForce(solutionOption, t, x, I1, fn[0]);

    //----------TIME-STEPPINGLOOP--------
    Real cpu0 = getCPU();
    for (int n = 0; n < numSteps; n++)
    {
        t = n * dt; // currenttime

        int cur = n % 2;        // currenttimelevel
        int next = (n + 1) % 2; // nexttimelevel

        //---AssigntheRHS---

        // getforcingatt+dt
        getForce(solutionOption, t + dt, x, I1, fn[next]);

        Real *uc_p = u[cur].getDataPointer();
#define UC(i) uc_p[i - nd1a]
        Real *un_p = u[next].getDataPointer();
#define UN(i) un_p[i - nd1a]

        Real *fc_p = fn[cur].getDataPointer();
#define FC(i) fc_p[i - nd1a]
        Real *fn_p = fn[next].getDataPointer();
#define FN(i) fn_p[i - nd1a]

        for (int i = n1a; i <= n1b; i++)
        {
            RHS(i) = UC(i) + (.5 * rx) * (UC(i + 1) - 2. * UC(i) + UC(i - 1)) + (.5 * dt) * (FC(i) + FN(i));
        }

        //----RHSforboundaryconditions---
        for (int side = 0; side <= 1; side++)
        {
            const inti1 = side == 0 ? n1a : n1b;
            const intis = 1 - 2 * side;
            IndexIb = Range(i1, i1);

            if (boundaryCondition(side, 0) == dirichlet)
            {
                getTrue(solutionOption, t + dt, x, Ib, uTrue);
                RHS(i1) = uTrue(i1);
            }
            else
            {
                // TheNeumannBC,(U(i1+1)-U(i1-1))/(2*dx)=g,
                // iscombinewiththeequationonthboundarytoeliminatethe
                // ghostpointvalue.Thisresultsinanadjustmenttotherhs.(seeclassnotes)
                getTruex(solutionOption, t + dt, x, Ib, uTrue); // storeuxinuTrue(i1)
                RHS(i1) = RHS(i1) - (is * dx * rx) * uTrue(i1);
            }
        }

        solveTridiagonal(Ax, rhsx);

        // fillinthesolution
        for (int i = n1a; i <= n1b; i++)
            UN(i) = RHS(i);

        // Fillghostpointvalues
        for (int side = 0; side <= 1; side++)
        {
            const int i1 = side == 0 ? n1a : n1b;
            const int is = 1 - 2 * side;
            IndexIb = Range(i1, i1);
            if (boundaryCondition(side, 0) == dirichlet)
            {
                // Useextrapolation:
                UN(i1 - is) = 3. * UN(i1) - 3. * UN(i1 + is) + UN(i1 + 2 * is);
            }
            else
            {
                // Neumann:(U(i1+1)-U(i1-1))/(2*dx)=g
                // getTruex(solutionOption,t+dt,x,Ib,uTrue);//noneedtore-evaluate
                UN(i1 - is) = UN(i1 + is) - 2. * is * dx * uTrue(i1);
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
            Real maxNorm = 0.;

            getTrue(solutionOption, t + dt, x, I1, uTrue);

            for (int i = nd1a; i <= nd1b; i++)
            {
                Real err = fabs(un(i) - uTrue(i));
                // Realerr=fabs(un(i)-uTrue(x(i),t+dt));
                maxErr = max(maxErr, err);
                maxNorm = max(maxNorm, abs(un(i)));
            }
            maxErr /= maxNorm;
            printf("step=%d,t=%9.3e,maxNorm=%9.2e,maxRelErr=%9.2e\n", n + 1, t + dt, maxNorm, maxErr);
        }
    }

    Real cpuTimeStep = getCPU() - cpu0;

    // checktheerror:
    t += dt; // tFinal;
    if (fabs(t - tFinal) > 1e-3 * dt / tFinal)
    {
        printf("ERROR:AFTERTIME_STEPPING:t=%16.8eISNOTEQUALtotFinal=%16.8e\n", t, tFinal);
    }

    RealArray err(Rx);

    cur = numSteps % 2;
    Real maxErr = 0.;
    Real maxNorm = 0.;

    getTrue(solutionOption, t, x, I1, uTrue);

    // Computefinalerror:donotinclude unusedghostpointonDirichletboundaries
    const int n1ae = boundaryCondition(0, 0) == dirichlet ? n1a : n1a - 1;
    const int n1be = boundaryCondition(1, 0) == dirichlet ? n1b : n1b + 1;
    for (int i = n1ae; i <= n1be; i++)
    {
        err(i) = uc(i) - uTrue(i);
        maxErr = max(maxErr, abs(err(i)));
        maxNorm = max(maxNorm, abs(uc(i)));
    }
    maxErr /= maxNorm;

    printf("numSteps=%4d,Nx=%3d,maxNorm=%9.2e,maxRelErr=%9.2e,cpu=%9.2e(s)\n",
           numSteps, Nx, maxNorm, maxErr, cpuTimeStep);

    //---Writeamatlabfileforplottinginmatlab--
    FILE *matlabFile = fopen(matlabFileName.c_str(), "w");
    fprintf(matlabFile, "%%Filewrittenbyheat1dImp.C\n");
    fprintf(matlabFile, "xa=%g;xb=%g;kappa=%g;t=%g;maxErr=%10.3e;cpuTimeStep=%10.3e;\n", xa, xb, kappa, tFinal, maxErr, cpuTimeStep);
    fprintf(matlabFile, "Nx=%d;dx=%14.6e;numGhost=%d;n1a=%d;n1b=%d;nd1a=%d;nd1b=%d;\n", Nx, dx, numGhost, n1a, n1b, nd1a, nd1b);
    fprintf(matlabFile, "solutionName=\’%s\’;\n", solutionName.c_str());

    writeMatlabVector(matlabFile, x, "x", nd1a, nd1b);
    writeMatlabVector(matlabFile, u[cur], "u", nd1a, nd1b);
    writeMatlabVector(matlabFile, err, "err", nd1a, nd1b);

    fclose(matlabFile);
    printf("Wrotefile[%s]\n", matlabFileName.c_str());

    return 0;
}