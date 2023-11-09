
#include <stdio.h>

int main(int args, char **)
{
    int n = 4;

    int *a = new int[n];
    int *b = new int[n];
    int *c = new int[n];

    for (int i = 0; i < n; i++)
    {
        a[i] = -i;
        b[i] = i * i;
    }

    for (int i = 0; i < n; i++)
        c[i] = a[i] + b[i];

    return 0;
}