//=====================================================================
//
// POISSONEQUATIONINTWODIMENSIONS
//-Delta(u)=f
//
// SolvewithA++arrays
//
//=====================================================================
#include "A++.h"

typedef double Real;
typedef doubleSerialArray RealArray;
typedef intSerialArray IntegerArray;

// include commandstpparsecommandlinearguments
#include "parseCommand.h"

// getCPU():Returnthecurrentwall-clocktimeinseconds
#include "getCPU.h"
#include <omp.h>
#include "math.h"

// enumforBC’s
enum BoundaryConditionsEnum
{
    periodic = -1,
    dirichlet = 1,
    neumann = 2
};

// storeparametersinthisclass
class PoissonParameters
{
public:
    Real dx[2];
    IntegerArray gridIndexRange;
    IntegerArray dimension;
    IntegerArray boundaryCondition;
    RealArray x;
    Real tol;
    RealArray uTrue;
    int maxIterations;
    int debug;
    int intervalToCheckResidual;
};

//-----------------------------------------------------------------------
// Returnthemax-normresidual
//-----------------------------------------------------------------------
Real getMaxResidual(RealArray &u, RealArray &f, PoissonParameters &par, int thread)
{
    const IntegerArray &gridIndexRange = par.gridIndexRange;
    const Real(&dx)[2] = par.dx; // Note:referencetoanarray

    const Real dxSq = dx[0] * dx[0];
    const Real dySq = dx[1] * dx[1];

    const int n1a = gridIndexRange(0, 0), n1b = gridIndexRange(1, 0);
    const int n2a = gridIndexRange(0, 1), n2b = gridIndexRange(1, 1);

    Real maxRes = 0.;
    int i1, i2;

#pragma omp parallel default(shared) num_threads(thread)
    {
#pragma omp for private(i2, i1) reduction(max : maxRes)
        for (i2 = n2a + 1; i2 <= n2b - 1; i2++)
            for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
            {
                Real res = f(i1, i2) + ((u(i1 + 1, i2) - 2. * u(i1, i2) + u(i1 - 1, i2)) / dxSq + (u(i1, i2 + 1) - 2. * u(i1, i2) + u(i1, i2 - 1)) / dySq);
                maxRes = max(maxRes, res);
            }
    }
    return maxRes;
}

//----------------------------------------------------------------------
// Returnthemax-normerror
//----------------------------------------------------------------------
Real getMaxError(RealArray &u, RealArray &err, PoissonParameters &par, int thread)
{
    const IntegerArray &gridIndexRange = par.gridIndexRange;
    const IntegerArray &dimension = par.dimension;
    const IntegerArray &boundaryCondition = par.boundaryCondition;
    const RealArray &x = par.x;
    const Real &tol = par.tol;
    const RealArray &uTrue = par.uTrue;
    const Real(&dx)[2] = par.dx; // Note:referencetoanarray
    const int &maxIterations = par.maxIterations;
    const int &debug = par.debug;

    const int n1a = gridIndexRange(0, 0), n1b = gridIndexRange(1, 0);
    const int n2a = gridIndexRange(0, 1), n2b = gridIndexRange(1, 1);

    Real maxErr = 0.;
    int i1, i2;
#pragma omp parallel default(shared) num_threads(thread)
    {
#pragma omp for private(i2, i1) reduction(max : maxErr)
        for (i2 = n2a; i2 <= n2b; i2++)
            for (i1 = n1a; i1 <= n1b; i1++)
            {
                Real xi = x(i1, i2, 0);
                Real yi = x(i1, i2, 1);
                err(i1, i2) = fabs(u(i1, i2) - uTrue(i1, i2));
                maxErr = max(err(i1, i2), maxErr);
            }
    }
    if (n1b - n1a + 1 <= 10)
    {
        u.display("u");
        err.display("err");
    }

    return maxErr;
}

#define F(i1, i2) f_p[(i1) + nd1 * (i2)]
#define U(i1, i2) u_p[(i1) + nd1 * (i2)]
#define UN(i1, i2) un_p[(i1) + nd1 * (i2)]

//--------------------------------------------
// Jacobiiteration
//--------------------------------------------
Real jacobiIteration(RealArray &u, RealArray &f, PoissonParameters &par, int thread)
{
    const IntegerArray &gridIndexRange = par.gridIndexRange;
    const IntegerArray &dimension = par.dimension;
    const IntegerArray &boundaryCondition = par.boundaryCondition;
    const RealArray &x = par.x;
    const Real &tol = par.tol;
    const RealArray &uTrue = par.uTrue;
    const Real(&dx)[2] = par.dx; // Note:referencetoanarray
    const int &maxIterations = par.maxIterations;
    const int &debug = par.debug;
    const int &intervalToCheckResidual = par.intervalToCheckResidual;

    const int n1a = gridIndexRange(0, 0), n1b = gridIndexRange(1, 0);
    const int n2a = gridIndexRange(0, 1), n2b = gridIndexRange(1, 1);
    const int nd1a = dimension(0, 0), nd1b = dimension(1, 0);
    const int nd2a = dimension(0, 1), nd2b = dimension(1, 1);
    const int nd1 = nd1b - nd1a + 1;

    const Real pi = M_PI;

    Range Rx(nd1a, nd1b), Ry(nd2a, nd2b);

    doubleArray ua[2];
    ua[0].redim(Rx, Ry);
    ua[0] = 0.;
    ua[1].redim(Rx, Ry);
    ua[1] = 0.;

    int current = 0;
    ua[current] = u; // initialguess

    int n;
    Real h = dx[0]; // assumedx=dy
    assert(dx[0] == dx[1]);
    Real omega = 1.;

    // Realres0=getMaxResidual(ua[current],f,h);
    Real res0 = getMaxResidual(ua[current], f, par, thread);
    Real maxRes = res0, maxResOld = res0;
    Real CR = 1.;

    const Real *f_p = f.getDataPointer();

    Real cpu1 = getCPU();
    for (n = 0; n < maxIterations; n++)
    {
        const int next = (current + 1) % 2;
        doubleArray &u = ua[current];
        doubleArray &un = ua[next];

        const Real *u_p = u.getDataPointer();
        Real *un_p = un.getDataPointer();

        // omega-JacobiIteration
        int i1, i2;
#pragma omp parallel default(shared) num_threads(thread)
        {
#pragma omp for private(i2, i1)
            for (i2 = n2a + 1; i2 <= n2b - 1; i2++)
                for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
                {
                    Real z = .25 * (h * h * F(i1, i2) + U(i1 + 1, i2) + U(i1 - 1, i2) + U(i1, i2 + 1) + U(i1, i2 - 1));
                    UN(i1, i2) = U(i1, i2) + omega * (z - U(i1, i2));
                }
        }
        current = next;

        if ((n % intervalToCheckResidual) == 0 || n == (maxIterations - 1))
        {
            // checkforconvergence
            maxRes = getMaxResidual(u, f, par, thread);
            CR = pow((maxRes / maxResOld), 1. / intervalToCheckResidual);
            maxResOld = maxRes;

            if (false)
                printf("n=%6d,maxRes=%9.3e\n", n, maxRes);

            if (maxRes < tol)
                break;
        }
    }
    Real cpu = getCPU() - cpu1;

    const int numIterations = n;
    // printf("numIterations=%d,nx=%d,ny=%d,cputime=%9.2e(s)\n",numIterations,nx,ny,cpu);
    if (maxRes < tol)
    {
        // printf("CONVERGENCE:maxRes=%8.2e<tol=%9.2e.\n",maxRes,tol);
    }
    else
        printf("Jacobi:ERROR:maxRes=%8.2e>tol=%9.2e:NOCONVERGENCE\n", maxRes, tol);

    //---computeerrors--
    u = ua[current];
    RealArray err(Rx, Ry);

    Real maxErr = getMaxError(u, err, par, thread);
    maxRes = getMaxResidual(u, f, par, thread);

    // Averageconvergencerate:
    Real aveCR = pow((maxRes / res0), 1. / max(1, numIterations));
    // Asymptoticconvergencerate:
    Real ACR = fabs(1. - omega * (1. - cos(pi * h)));
    printf("Jac:omega=%4.2fIts=%6dres=%8.2eerr=%8.2eCR=%7.5faveCR=%7.5fACR=%7.5fcpu=%7.1escpu/it=%7.1es\n",
           omega, numIterations, maxRes, maxErr, CR, aveCR, ACR, cpu, cpu / numIterations);

    return cpu;
}

//--------------------------------------------
// Gauss-SeidelorSOR
//--------------------------------------------
Real gaussSeidelIteration(RealArray &u, RealArray &f, PoissonParameters &par, int thread)
{
    const IntegerArray &gridIndexRange = par.gridIndexRange;
    const IntegerArray &dimension = par.dimension;
    const IntegerArray &boundaryCondition = par.boundaryCondition;
    const RealArray &x = par.x;
    const Real &tol = par.tol;
    const RealArray &uTrue = par.uTrue;
    const Real(&dx)[2] = par.dx; // Note:referencetoanarray
    const int &maxIterations = par.maxIterations;
    const int &debug = par.debug;
    const int &intervalToCheckResidual = par.intervalToCheckResidual;

    const int n1a = gridIndexRange(0, 0), n1b = gridIndexRange(1, 0);
    const int n2a = gridIndexRange(0, 1), n2b = gridIndexRange(1, 1);
    const int nd1a = dimension(0, 0), nd1b = dimension(1, 0);
    const int nd2a = dimension(0, 1), nd2b = dimension(1, 1);
    const int nd1 = nd1b - nd1a + 1;
    Range Rx(nd1a, nd1b), Ry(nd2a, nd2b);

    const Real pi = M_PI;

    int n;
    Real h = dx[0]; // assumedx=dy
    assert(dx[0] == dx[1]);

    Real omega = 2. / (1. + sin(pi * h)); // optimalomega

    Real res0 = getMaxResidual(u, f, par, thread);
    Real maxRes = res0, maxResOld = maxRes, CR = 1.;

    const Real *f_p = f.getDataPointer();
    Real *u_p = u.getDataPointer();

    Real cpu1 = getCPU();
    for (n = 0; n < maxIterations; n++)
    {
        // omega-GSIteration
        int i1, i2;
#pragma omp parallel default(shared) num_threads(thread)
        {
#pragma omp for private(i2, i1)
            for (i2 = n2a + 1; i2 <= n2b - 1; i2++)
                for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
                {
                    Real z = .25 * (h * h * F(i1, i2) + U(i1 + 1, i2) + U(i1 - 1, i2) + U(i1, i2 + 1) + U(i1, i2 - 1));
                    U(i1, i2) = U(i1, i2) + omega * (z - U(i1, i2));
                }
        }
        if ((n % intervalToCheckResidual) == 0 || n == (maxIterations - 1))
        {
            // checkforconvergence
            maxRes = getMaxResidual(u, f, par, thread);

            CR = pow((maxRes / maxResOld), 1. / intervalToCheckResidual);
            maxResOld = maxRes;
            if (false)
            {
                printf("GS:n=%6d,maxRes=%9.3e,CR=%7.4f\n", n, maxRes, CR);
            }

            if (maxRes < tol)
                break;
        }
    }
    Real cpu = getCPU() - cpu1;

    const int numIterations = n;
    // printf("GS:numIterations=%d,nx=%d,ny=%d,cputime=%9.2e(s)\n",numIterations,nx,ny,cpu);
    if (maxRes < tol)
    {
        // printf("GS:CONVERGENCE:maxRes=%8.2e<tol=%9.2e.\n",maxRes,tol);
    }
    else
        printf("GS:ERROR:maxRes=%8.2e>tol=%9.2e:NOCONVERGENCE\n", maxRes, tol);

    //---computeerrors--
    RealArray err(Rx, Ry);
    Real maxErr = getMaxError(u, err, par, thread);
    maxRes = getMaxResidual(u, f, par, thread);

    // Averageconvergencerate:
    Real aveCR = pow((maxRes / res0), 1. / max(1, numIterations));
    // Asymptoticconvergencerate:
    Real ACR = omega - 1.;
    printf("GS: omega=%4.2f Its=%6d res=%8.2e err=%8.2e CR=%7.5f aveCR=%7.5f ACR=%7.5f cpu=%7.1e scpu/it=%7.1e s\n",
           omega, numIterations, maxRes, maxErr, CR, aveCR, ACR, cpu, cpu / numIterations);

    return cpu;
}

//--------------------------------------------
// Red-BlackGauss-Seideliteration
//--------------------------------------------
// intredBlackIteration(RealArray&u,RealArray&f,Realxa,Realya,Realdx[],intmaxIterations,Realtol,RealArray&uTrue)
Real redBlackIteration(RealArray &u, RealArray &f, PoissonParameters &par, int thread)
{

    const IntegerArray &gridIndexRange = par.gridIndexRange;
    const IntegerArray &dimension = par.dimension;
    const IntegerArray &boundaryCondition = par.boundaryCondition;
    const RealArray &x = par.x;
    const Real &tol = par.tol;
    const RealArray &uTrue = par.uTrue;
    const Real(&dx)[2] = par.dx; // Note:referencetoanarray
    const int &maxIterations = par.maxIterations;
    const int &debug = par.debug;
    const int &intervalToCheckResidual = par.intervalToCheckResidual;

    const int n1a = gridIndexRange(0, 0), n1b = gridIndexRange(1, 0);
    const int n2a = gridIndexRange(0, 1), n2b = gridIndexRange(1, 1);
    const int nd1a = dimension(0, 0), nd1b = dimension(1, 0);
    const int nd2a = dimension(0, 1), nd2b = dimension(1, 1);
    const int nd1 = nd1b - nd1a + 1;
    Range Rx(nd1a, nd1b), Ry(nd2a, nd2b);

    const Real pi = M_PI;

    doubleArray ua[2];
    ua[0].redim(Rx, Ry);
    ua[0] = 0.;
    ua[1].redim(Rx, Ry);
    ua[1] = 0.;

    int current = 0;
    ua[0] = u; // initialguess
    ua[1] = u; // initialguess

    int n;
    Real h = dx[0]; // assumedx=dy
    assert(dx[0] == dx[1]);

    // thisistheoptimalomegaforANYordering(Owlbookp32)
    const Real omega = 2. / (1. + sin(pi * h));

    Real res0 = getMaxResidual(ua[current], f, par, thread);
    Real maxRes = res0, maxResOld = res0;
    Real CR = 1.;

    const Real *f_p = f.getDataPointer();

    Real cpu1 = getCPU();
    for (n = 0; n < maxIterations; n++)
    {
        const int next = (current + 1) % 2;
        doubleArray &u = ua[current];
        doubleArray &un = ua[next];

        Real *u_p = u.getDataPointer();
        Real *un_p = un.getDataPointer();

        //---omegaRed-Black--
        // Question:whatisthefastestwaytoloopoverredandblackpoints?
        //(1)Useifstatement
        //(2)Useloopswithstrideof2
        // Note:loopsarefasterifLHSarrayisdifferentfromRHSarrays
        int i1, i2;
#pragma omp parallel default(shared) num_threads(thread)
        {
#pragma omp for private(i2, i1)
            for (i2 = n2a + 1; i2 <= n2b - 1; i2++)
                for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
                {
                    if ((i1 + i2) % 2 == 0) // redpoints
                    {
                        Real z = .25 * (h * h * F(i1, i2) + U(i1 + 1, i2) + U(i1 - 1, i2) + U(i1, i2 + 1) + U(i1, i2 - 1));
                        UN(i1, i2) = U(i1, i2) + omega * (z - U(i1, i2));
                    }
                    else
                    {
                        UN(i1, i2) = U(i1, i2);
                    }
                }
        }
#pragma omp parallel default(shared) num_threads(thread)
        {
#pragma omp for private(i2, i1)
            for (int i2 = n2a + 1; i2 <= n2b - 1; i2++)
                for (int i1 = n1a + 1; i1 <= n1b - 1; i1++)
                {
                    if ((i1 + i2) % 2 == 1) // blackpoints
                    {
                        Real z = .25 * (h * h * F(i1, i2) + UN(i1 + 1, i2) + UN(i1 - 1, i2) + UN(i1, i2 + 1) + UN(i1, i2 - 1));
                        U(i1, i2) = UN(i1, i2) + omega * (z - UN(i1, i2));
                    }
                    else
                    {
                        U(i1, i2) = UN(i1, i2);
                    }
                }
        }
        if ((n % intervalToCheckResidual) == 0 || n == (maxIterations - 1))
        {
            // checkforconvergence
            maxRes = getMaxResidual(u, f, par, thread);
            CR = pow((maxRes / maxResOld), 1. / intervalToCheckResidual);
            maxResOld = maxRes;

            if (false)
                printf("RB:n=%6d,maxRes=%9.3e,CR=%7.5f\n", n, maxRes, CR);

            if (maxRes < tol)
                break;
        }
    }
    Real cpu = getCPU() - cpu1;

    const int numIterations = n;
    // printf("numIterations=%d,nx=%d,ny=%d,cputime=%9.2e(s)\n",numIterations,nx,ny,cpu);
    if (maxRes < tol)
    {
        // printf("CONVERGENCE:maxRes=%8.2e<tol=%9.2e.\n",maxRes,tol);
    }
    else
        printf("Red-Black:ERROR:maxRes=%8.2e>tol=%9.2e:NOCONVERGENCE\n", maxRes, tol);

    //---computeerrors--
    u = ua[current];
    RealArray err(Rx, Ry);
    Real maxErr = getMaxError(u, err, par, thread);
    maxRes = getMaxResidual(u, f, par, thread);

    // Averageconvergencerate:
    Real aveCR = pow((maxRes / res0), 1. / max(1, numIterations));
    // Asymptoticconvergencerate:
    Real ACR = omega - 1.;
    printf("RB: omega=%4.2f Its=%6d res=%8.2e err=%8.2e CR=%7.5faveCR=%7.5fACR=%7.5fcpu=%7.1escpu/it=%7.1es\n",
           omega, numIterations, maxRes, maxErr, CR, aveCR, ACR, cpu, cpu / numIterations);

    return cpu;
}

//--------------------------------------------
// ConjugateGradientiteration
//--------------------------------------------
Real conjugateGradientIteration(RealArray &u, RealArray &f, PoissonParameters &par, int thread)
{
    const IntegerArray &gridIndexRange = par.gridIndexRange;
    const IntegerArray &dimension = par.dimension;
    const IntegerArray &boundaryCondition = par.boundaryCondition;
    const RealArray &x = par.x;
    const Real &tol = par.tol;
    const RealArray &uTrue = par.uTrue;
    const Real(&dx)[2] = par.dx; // Note:referencetoanarray
    const int &maxIterations = par.maxIterations;
    const int &debug = par.debug;
    const int &intervalToCheckResidual = par.intervalToCheckResidual;

    const int n1a = gridIndexRange(0, 0), n1b = gridIndexRange(1, 0);
    const int n2a = gridIndexRange(0, 1), n2b = gridIndexRange(1, 1);
    const int nd1a = dimension(0, 0), nd1b = dimension(1, 0);
    const int nd2a = dimension(0, 1), nd2b = dimension(1, 1);
    const int nd1 = nd1b - nd1a + 1;

    const Real dxSq = dx[0] * dx[0];
    const Real dySq = dx[1] * dx[1];

    Range Rx(nd1a, nd1b), Ry(nd2a, nd2b);

    // Temporaryarrays
    RealArray z(Rx, Ry), r(Rx, Ry), p(Rx, Ry);

    const Real *f_p = f.getDataPointer();

    Real *u_p = u.getDataPointer();
    Real *z_p = z.getDataPointer();
    Real *r_p = r.getDataPointer();
    Real *p_p = p.getDataPointer();
#define Z(i1, i2) z_p[(i1) + nd1 * (i2)]
#define R(i1, i2) r_p[(i1) + nd1 * (i2)]
#define P(i1, i2) p_p[(i1) + nd1 * (i2)]

    Real cpu1 = getCPU();

    u = 0.; // initialguess--***FIXMEfornon-zeroinitialguess***
    r = 0;  // initialresidualgoeshere

    int i1, i2;
#pragma omp parallel default(shared) num_threads(thread)
    {
#pragma omp for private(i2, i1)

        for (i2 = n2a + 1; i2 <= n2b - 1; i2++)
            for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
            {
                R(i1, i2) = F(i1, i2); // initialresidual
                P(i1, i2) = R(i1, i2); // initialsearchdirection
            }
    }
    Real alpha, beta, rNormSquared, rNormSquaredNew, pz;

    // rNormSquared=r^Tr
    rNormSquared = 0.;
#pragma omp parallel default(shared) num_threads(thread)
    {
#pragma omp for private(i2, i1) reduction(+ : rNormSquared)
        for (i2 = n2a + 1; i2 <= n2b - 1; i2++)
            for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
            {
                rNormSquared += R(i1, i2) * R(i1, i2);
            }
    }
    Real res0 = sqrt(rNormSquared);
    Real maxRes = res0, maxResOld = res0;
    Real CR = 1.;

    // printf("CG:res0=%9.2e,rNormSquared=%9.2e\n",res0,rNormSquared);

    int n, nOld = -1;
    for (n = 0; n < maxIterations; n++)
    {

        // Note:Someloopshavebeencombinedforspeed

        // z=Ap
        // pz=p^Tz
        // Note:A=-Delta
        pz = 0.;
#pragma omp parallel default(shared) num_threads(thread)
        {
#pragma omp for private(i2, i1) reduction(+ : pz)
            for (i2 = n2a + 1; i2 <= n2b - 1; i2++)
                for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
                {
                    Z(i1, i2) = -((P(i1 + 1, i2) - 2. * P(i1, i2) + P(i1 - 1, i2)) / dxSq + (P(i1, i2 + 1) - 2. * P(i1, i2) + P(i1, i2 - 1)) / dySq);
                    pz += P(i1, i2) * Z(i1, i2);
                }
        }
        ////pz=p^Tz
        // pz=0.;
        // for(inti2=n2a+1;i2<=n2b-1;i2++)
        // for(inti1=n1a+1;i1<=n1b-1;i1++)
        //{
        // pz+=P(i1,i2)*Z(i1,i2);
        // }

        alpha = rNormSquared / pz; // steplength

        // Mergingtheseintooneloopisfasterinserial
        // u+=alpha*p
        // r-=alpha*z
        // rNormSquaredNew=r^Tr
        rNormSquaredNew = 0.;
#pragma omp parallel default(shared) num_threads(thread)
        {
#pragma omp for private(i2, i1) reduction(+ : rNormSquaredNew)
            for (int i2 = n2a + 1; i2 <= n2b - 1; i2++)
                for (int i1 = n1a + 1; i1 <= n1b - 1; i1++)
                {
                    U(i1, i2) += alpha * P(i1, i2); // newsolution
                    R(i1, i2) -= alpha * Z(i1, i2); // residual
                    rNormSquaredNew += R(i1, i2) * R(i1, i2);
                }
        }
        ////r-=alpha*z
        // rNormSquaredNew=0.;
        // for(inti2=n2a+1;i2<=n2b-1;i2++)
        // for(inti1=n1a+1;i1<=n1b-1;i1++)
        //{
        // R(i1,i2)-=alpha*Z(i1,i2);//residual
        // rNormSquaredNew+=R(i1,i2)*R(i1,i2);
        // }

        ////rNormSquaredNew=r^Tr*MERGEWITHPREVIOUSLOOP??
        // rNormSquaredNew=0.;
        // for(inti2=n2a+1;i2<=n2b-1;i2++)
        // for(inti1=n1a+1;i1<=n1b-1;i1++)
        //{
        // rNormSquaredNew+=R(i1,i2)*R(i1,i2);
        // }

        beta = rNormSquaredNew / rNormSquared; // improvementthisstep

        rNormSquared = rNormSquaredNew;

        // p=r+beta*p
#pragma omp parallel default(shared) num_threads(thread)
        {
#pragma omp for private(i2, i1)
            for (i2 = n2a + 1; i2 <= n2b - 1; i2++)
                for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
                {
                    P(i1, i2) = R(i1, i2) + beta * P(i1, i2); // newsearchdirection
                }
        }
        if ((n % intervalToCheckResidual) == 0 || n == (maxIterations - 1))
        {
            // checkforconvergence****WESHOULDREALLYUSEEXISTING2-normRESIDUAL****
            // Dothisfornowtobeconsistentwithotherschemes.
            maxRes = getMaxResidual(u, f, par, thread);
            CR = pow((maxRes / maxResOld), 1. / max(1, n - nOld));
            maxResOld = maxRes;
            nOld = n;

            if (debug > 0)
                printf("CG:n=%3d:pz=%9.2e,alpha=%9.2e,beta=%9.2e,maxRes=%8.2e,CR=%7.5f\n", n, pz, alpha, beta, maxRes, CR);

            if (maxRes < tol)
                break;
        }
    }
    Real cpu = getCPU() - cpu1;

    const int numIterations = n;
    // printf("numIterations=%d,nx=%d,ny=%d,cputime=%9.2e(s)\n",numIterations,nx,ny,cpu);
    if (maxRes < tol)
    {
        // printf("CONVERGENCE:maxRes=%8.2e<tol=%9.2e.\n",maxRes,tol);
    }
    else
        printf("CG:ERROR:maxRes=%8.2e>tol=%9.2e:NOCONVERGENCE\n", maxRes, tol);

    //---computeerrors--
    RealArray err(Rx, Ry);
    Real maxErr = getMaxError(u, err, par, thread);
    maxRes = getMaxResidual(u, f, par, thread);

    // Averageconvergencerate:
    Real aveCR = pow((maxRes / res0), 1. / max(1, numIterations));
    // Asymptoticconvergencerate:
    const Real pi = M_PI;
    Real ACR = 1. - pi * dx[0]; //**CHECKME**
    printf("CG:Its=%6dres=%8.2eerr=%8.2eCR=%7.5faveCR=%7.5fACR=%7.5fcpu=%7.1escpu/it=%7.1es\n",
           numIterations, maxRes, maxErr, CR, aveCR, ACR, cpu, cpu / numIterations);

    return cpu;
}

int main(int argc, char *argv[])
{

    printf("Usage:poisson-nx=<i>-tol=<f>-maxIterations=<i>-debug=<i>\n"
           "nx=numberofgridpointsinxandydirections\n"
           "tol=convergencetolerance\n"
           "maxIterations=maxnumberofiterations\n");

    const Real pi = M_PI; // 4.*atan2(1.,1.);

    // Parametersarestoredhere:
    PoissonParameters par;

    // Makereferencestoparametersforclarity
    IntegerArray &gridIndexRange = par.gridIndexRange;
    IntegerArray &dimension = par.dimension;
    IntegerArray &boundaryCondition = par.boundaryCondition;
    RealArray &x = par.x;
    Real &tol = par.tol;
    RealArray &uTrue = par.uTrue;
    Real(&dx)[2] = par.dx; // Note:referencetoanarray
    int &maxIterations = par.maxIterations;
    int &debug = par.debug;
    int &intervalToCheckResidual = par.intervalToCheckResidual;

    intervalToCheckResidual = 50; // checkresidualeverythismanyiterations

    const int numberOfDimensions = 2;
    Real xa = 0., xb = 1.; // domainis[xa,xb]X[ya,yb]
    Real ya = 0., yb = 1.;

    debug = 0;
    maxIterations = 1000;
    tol = 1.e-3;
    int nx = 100, ny = nx;
    int thread = 1;
    int scallingTest = 0;

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
        else if (parseCommand(line, "-maxIterations=", maxIterations))
        {
        }
        else if (parseCommand(line, "-tol=", tol))
        {
        }
        else if (parseCommand(line, "-threads=", thread))
        {
        }
        else if (parseCommand(line, "-scallingTest=", scallingTest))
        {
        }
    }

    printf("-----------------------------------------------------------------\n");
    printf("---------SolvethePoissonEquationintwodimensions----------\n");
    printf("-Delta(u)=f\n");
    printf("DirichletBCs\n");
    printf("nx=%d, ny=%d, maxIterations=%d, tol=%9.2e\n", nx, ny, maxIterations, tol);
    printf("-----------------------------------------------------------------\n");

    const Real kx = 1., ky = 1.;
    const Real kxp = kx * pi;
    const Real kyp = ky * pi;

    // ThistruesolutionistooeasyforCGsinceuisamultipleoff:
    // #define UTRUE(x,y)sin(kxp*(x))*sin(kyp*(y))
    ////f=-Delta(u)
    // #define FORCE(x,y)(kxp*kxp+kyp*kyp)*UTRUE(x,y)

#define UTRUE(x, y) (x) * (1. - (x)) * (y) * (1. - (y))
// f=-Delta(u)
#define FORCE(x, y) 2. * ((x) * (1. - (x)) + (y) * (1. - (y)))

    const int numGhost = 0;
    const int n1a = numGhost;
    const int n1b = n1a + nx;
    const int nd1a = n1a - numGhost;
    const int nd1b = n1b + numGhost;
    const int nd1 = nd1b - nd1a + 1;

    const int n2a = numGhost;
    const int n2b = n2a + ny;
    const int nd2a = n2a - numGhost;
    const int nd2b = n2b + numGhost;
    const int nd2 = nd2b - nd2a + 1;

    gridIndexRange.redim(2, numberOfDimensions);
    dimension.redim(2, numberOfDimensions);
    boundaryCondition.redim(2, numberOfDimensions);

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
    x.redim(Rx, Ry, 2);

    // Realdx[2];
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
    RealArray f(Rx, Ry);
    uTrue.redim(Rx, Ry);

    RealArray u(Rx, Ry);

    Real xi, yi;
#pragma omp parallel default(shared) num_threads(thread)
    {
#pragma omp for private(i2, i1, xi, yi)
        for (i2 = nd2a; i2 <= nd2b; i2++)
            for (i1 = nd1a; i1 <= nd1b; i1++)
            {
                xi = x(i1, i2, 0);
                yi = x(i1, i2, 1);
                f(i1, i2) = FORCE(xi, yi);
                uTrue(i1, i2) = UTRUE(xi, yi);
            }
    }
    if (scallingTest == 0)
    {
        //=====JACOBI=====
        u = 0.; // initialguess
        jacobiIteration(u, f, par, thread);

        //=====Gauss-Seidel=====
        u = 0.; // initialguess
        gaussSeidelIteration(u, f, par, thread);

        //=====Red-BlackGauss-Seidel=====
        u = 0.; // initialguess
        redBlackIteration(u, f, par, thread);

        //====ConjugateGradient=========
        u = 0.; // initialguess
        conjugateGradientIteration(u, f, par, thread);
    }
    if (scallingTest > 0)
    {
        Real firstTime;
        char s[1000];
        for (int i = 0; i < 7; i++)
        {
            Real time = 0;

            u = 0.;
            if (scallingTest == 1)
            {
                time = jacobiIteration(u, f, par, int(pow(2, i)));
            }
            if (scallingTest == 2)
            {
                time = gaussSeidelIteration(u, f, par, int(pow(2, i)));
            }
            if (scallingTest == 3)
            {
                time = redBlackIteration(u, f, par, int(pow(2, i)));
            }
            if (scallingTest == 4)
            {
                time = conjugateGradientIteration(u, f, par, int(pow(2, i)));
            }
            if (i == 0)
                firstTime = time;

            sprintf(s, "%s\n%d  %f  %f  %f", s, int(pow(2, i)), time, firstTime / time, firstTime / time / int(pow(2, i)));
        }
        printf("%s\n", s);
    }

    return 0;
}