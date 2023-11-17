
// export PATH=/usr/local/cuda/bin:$PATH

#include <stdio.h>
#include <iostream>
#include <ctime>

using namespace std;

static void HandleError(cudaError_t err, const char *file, int line)
{
    if (err != cudaSuccess)
    {
        printf("%s in %s atline %d\n", cudaGetErrorString(err), file, line);
        exit(EXIT_FAILURE);
    }
}
#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))

__global__ void add_vector(long int *a, long int *b, long int *c, long int *n)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    i = i % *n;
    if (i < *n)
        c[i] = a[i] + b[i];
}

int main(int args, char **argv)
{
    long int n = atoi(argv[1]);
    long int Nt = 1;
    float first_time;

    for (int count = 0; count <= 10; count++)
    {
        long int Nb = ceil((float)n / Nt);

        long int *a = new long int[n];
        long int *b = new long int[n];
        long int *c = new long int[n];

        for (long int i = 0; i < n; i++)
        {
            a[i] = -i;
            b[i] = i * i;
            c[i] = 0;
        }

        // Cuda code
        long int *gpu_a;
        long int *gpu_b;
        long int *gpu_c;
        long int *gpu_n;

        HANDLE_ERROR(cudaMalloc((void **)&gpu_a, n * sizeof(long int)));

        HANDLE_ERROR(cudaMalloc((void **)&gpu_b, n * sizeof(long int)));
        HANDLE_ERROR(cudaMalloc((void **)&gpu_c, n * sizeof(long int)));
        HANDLE_ERROR(cudaMalloc((void **)&gpu_n, sizeof(long int)));

        HANDLE_ERROR(cudaMemcpy(gpu_a, a, n * sizeof(long int), cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(gpu_b, b, n * sizeof(long int), cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(gpu_c, c, n * sizeof(long int), cudaMemcpyHostToDevice));
        HANDLE_ERROR(cudaMemcpy(gpu_n, &n, sizeof(long int), cudaMemcpyHostToDevice));

        cudaEvent_t start;
        cudaEvent_t end;
        float time = 0;

        cudaEventCreate(&start);
        cudaEventCreate(&end);

        cudaEventRecord(start, 0);
        add_vector<<<Nb, Nt>>>(gpu_a, gpu_b, gpu_c, gpu_n);
        cudaEventRecord(end, 0);

        cudaDeviceSynchronize();
        cudaEventElapsedTime(&time, start, end);

        HANDLE_ERROR(cudaMemcpy(c, gpu_c, n * sizeof(long int), cudaMemcpyDeviceToHost));

        long int error = 0;

        if (count == 0)
            first_time = (1.0 * std::clock()) / CLOCKS_PER_SEC;

        for (long int i = 0; i < n; i++)
        {
            c[i] = c[i] - a[i] - b[i];
            error = error + c[i];
        }

        if (count == 0)
            first_time = (1.0 * std::clock()) / CLOCKS_PER_SEC - first_time;

        cout << Nb << " " << Nt << " " << Nb * Nt << " " << n << " " << error << " " << time << " " << first_time / time * 1000 << endl;

        free(a);
        free(b);
        free(c);

        cudaFree(gpu_a);
        cudaFree(gpu_b);
        cudaFree(gpu_c);

        Nt = Nt * 2;
    }
    return 0;
}
