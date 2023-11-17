#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <ctime>
#include <string>
#include <iostream>

using std::max;
using std::string;
using namespace std;

typedef double Real;

#define x(i) x_p[i - nd1a]
#define boundaryCondition(side, axis) boundaryCondition_p[(side) + 2 * (axis)]
#define UTRUE(x, t) (b0 + (x) * (b1 + (x) * b2)) * (a0 + (t) * (a1))
#define UTRUEX(x, t) (b1 + 2. * (x) * b2) * (a0 + (t) * (a1))
#define UTRUET(x, t) (b0 + (x) * (b1 + (x) * b2)) * (a1)
#define UTRUEXX(x, t) (2. * b2) * (a0 + (t) * (a1))
#define FORCE(x, t) (UTRUET(x, t) - kappa * UTRUEXX(x, t))

static void HandleError(cudaError_t err, const char *file, int line)
{
    if (err != cudaSuccess)
    {
        printf("%s in %s atline %d\n", cudaGetErrorString(err), file, line);
        exit(EXIT_FAILURE);
    }
}

#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))

__global__ void time_step(Real *u_c, Real *u_n, Real *x_p, Real *params)
{
#define uc(i) u_c[i - nd1a]
#define un(i) u_n[i - nd1a]

    Real b0 = 1., b1 = .5, b2 = .25;
    Real a0 = 1., a1 = .3, kappa = .1;

    Real rx = params[0];
    Real t = params[1];
    Real dt = params[2];
    int n1a = params[3];
    int n1b = params[4];
    int nd1a = params[5];
    int nd1b = params[6];

    int i = n1a + threadIdx.x + blockIdx.x * blockDim.x;
    int n = n1b - n1a + 1;

    if (i < n)
    {
        un(i) = uc(i) + rx * (uc(i + 1) - 2. * uc(i) + uc(i - 1)) + dt * FORCE(x(i), t);
    }
}

__global__ void boundary_conditions(Real *u_c, Real *u_n, Real *x_p, Real *params)
{
    Real b0 = 1., b1 = .5, b2 = .25;
    Real a0 = 1., a1 = .3, kappa = .1;

    Real rx = params[0];
    Real t = params[1];
    Real dt = params[2];
    int n1a = params[3];
    int n1b = params[4];
    int nd1a = params[5];
    int nd1b = params[6];

    for (int side = 0; side <= 1; side++)
    {
        const int i = side == 0 ? n1a : n1b; // boundaryindex
        const int is = 1 - 2 * side;         // is=1onleft,-1onright

        un(i) = UTRUE(x(i), t + dt);
        un(i - is) = UTRUE(x(i - is), t + dt); // extrapolateghost
    }
    params[1] = t + dt;

    return;
}

__global__ void update_current(Real *u_c, Real *u_n, Real *params)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int nd1 = params[7];
    if (i < nd1)
    {
        u_c[i] = u_n[i];
    }
}

#define uc(i) u_p[cur][i - nd1a]
#define un(i) u_p[next][i - nd1a]

int main(int argc, char *argv[])
{
    // Setup

    const Real pi = M_PI;

    int debug = 0; // setto1fordebuginfo
    Real xa = 0., xb = 1.;
    Real kappa = .1;
    Real tFinal = atof(argv[2]);
    Real cfl = .9; // time-stepsafetyfactor

    int Nx = 10; // default
    string matlabFileName = "heat1d.m";

    if (argc >= 2) // readanycommandlinearguments
    {
        Nx = atoi(argv[1]);
        printf("SettingNx=%d\n", Nx);
    }

    Real dx = (xb - xa) / Nx;
    const int numGhost = 1;
    const int n1a = 0;
    const int n1b = Nx;
    const int nd1a = n1a - numGhost;
    const int nd1b = n1b + numGhost;
    const int nd1 = nd1b - nd1a + 1; // totalnumberofgridpoints;

    Real *x_p = new Real[nd1];

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

    const Real kx = 3.;
    const Real kxPi = kx * pi;
    const Real kappaPiSq = kappa * kxPi * kxPi;

    const char solutionName[] = "polyDD";
    boundaryCondition(0, 0) = dirichlet;
    boundaryCondition(1, 0) = dirichlet;

    Real b0 = 1., b1 = .5, b2 = .25;
    Real a0 = 1., a1 = .3;

    Real *u_p[2];
    u_p[0] = new Real[nd1];
    u_p[1] = new Real[nd1];

    // initialconditions
    Real t = 0.;
    int cur = 0; //"current"solution,indexintou_p[]
    for (int i = nd1a; i <= nd1b; i++)
        uc(i) = UTRUE(x(i), t);

    const Real dx2 = dx * dx;
    Real dt = cfl * .5 * dx2 / kappa; // dt,adjustedbelow
    const int numSteps = ceil(tFinal / dt);
    dt = tFinal / numSteps; // adjustdttoreachthefinaltime
    const Real rx = kappa * dt / dx2;

    // Declare GPU vars
    Real *gpu_u_p_c;
    Real *gpu_u_p_n;
    Real *gpu_x_p;

    Real params[8];
    params[0] = rx;
    params[1] = t;
    params[2] = dt;
    params[3] = n1a;
    params[4] = n1b;
    params[5] = nd1a;
    params[6] = nd1b;
    params[7] = nd1;

    Real *gpu_params;

    HANDLE_ERROR(cudaMalloc((void **)&gpu_u_p_c, nd1 * sizeof(Real)));
    HANDLE_ERROR(cudaMalloc((void **)&gpu_u_p_n, nd1 * sizeof(Real)));
    HANDLE_ERROR(cudaMalloc((void **)&gpu_x_p, nd1 * sizeof(Real)));
    HANDLE_ERROR(cudaMalloc((void **)&gpu_params, 8 * sizeof(Real)));

    HANDLE_ERROR(cudaMemcpy(gpu_u_p_c, u_p[0], nd1 * sizeof(Real), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(gpu_u_p_n, u_p[1], nd1 * sizeof(Real), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(gpu_x_p, x_p, nd1 * sizeof(Real), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(gpu_params, params, 8 * sizeof(Real), cudaMemcpyHostToDevice));

    cudaEvent_t start;
    cudaEvent_t end;
    float time = 0;

    cudaEventCreate(&start);
    cudaEventCreate(&end);

    cudaEventRecord(start, 0);

    // Run Time Stepping Loop
    for (int n = 0; n < numSteps; n++)
    {
        t = n * dt;

        const int cur = n % 2;
        const int next = (n + 1) % 2;

        int Nt = 128;
        long int Nb = ceil((float)Nx / Nt);

        time_step<<<Nb, Nt>>>(gpu_u_p_c, gpu_u_p_n, gpu_x_p, gpu_params);
        cudaDeviceSynchronize();
        boundary_conditions<<<1, 1>>>(gpu_u_p_c, gpu_u_p_n, gpu_x_p, gpu_params);
        cudaDeviceSynchronize();
        Nt = 128;
        Nb = ceil((float)nd1 / Nt);
        update_current<<<Nb, Nt>>>(gpu_u_p_c, gpu_u_p_n, gpu_params);
        cudaDeviceSynchronize();
    }
    t += dt;

    cudaEventRecord(end, 0);

    cudaDeviceSynchronize();
    cudaEventElapsedTime(&time, start, end);

    // Bring back data

    Real *u_p_new = new Real[nd1];

    HANDLE_ERROR(cudaMemcpy(u_p_new, gpu_u_p_n, nd1 * sizeof(Real), cudaMemcpyDeviceToHost));

    // Ending stuff

    Real cpuTimeStep = time / 1000;
    Real *error_p = new Real[nd1];
#define error(i) error_p[i - nd1a]

    cur = 0;
    Real maxErr = 0.;
    for (int i = nd1a; i <= nd1b; i++)
    {
        error(i) = u_p_new[i + 1] - UTRUE(x(i), t);
        maxErr = max(maxErr, abs(error(i)));
    }

    printf("numSteps=%4d,Nx=%3d,maxErr=%9.2e,gpu=%9.2e(s)\n", numSteps, Nx, maxErr, cpuTimeStep);
}