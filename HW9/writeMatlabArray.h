#ifndef WRITE_MATLAB_ARRAY_H
#define WRITE_MATLAB_ARRAY_H

//-----------------------------------------------------------------------------------
// Saveanarraytoamatlabfile.
//-----------------------------------------------------------------------------------
int writeMatlabArray(FILE *matlabFile, RealArray &u, const char *name, int numberOfComponents, IntegerArray &dimension)
{
    const int nd1a = dimension(0, 0);
    const int nd1b = dimension(1, 0);
    const int nd2a = dimension(0, 1);
    const int nd2b = dimension(1, 1);
    const int nd1 = nd1b - nd1a + 1;
    const int nd2 = nd2b - nd2a + 1;

    const int numPerLine = 8; // numberofentriesperline
    if (numberOfComponents == 1)
    {
        fprintf(matlabFile, "%s=zeros(%d,%d);\n", name, nd1, nd2);

        int count = 0;
        for (int i2 = nd1a; i2 <= nd1b; i2++)
            for (int i1 = nd1a; i1 <= nd1b; i1++)
            {
                fprintf(matlabFile, "%s(%3d,%3d)=%12.5e;", name, i1 - nd1a + 1, i2 - nd2a + 1, u(i1, i2));
                if (count % numPerLine == numPerLine - 1)
                    fprintf(matlabFile, "\n"); // newline
                count++;
            }
        fprintf(matlabFile, "\n");
    }
    else
    {
        fprintf(matlabFile, "%s=zeros(%d,%d,%d);\n", name, nd1, nd2, numberOfComponents);

        int count = 0;
        for (int m = 0; m < numberOfComponents; m++)
            for (int i2 = nd1a; i2 <= nd1b; i2++)
                for (int i1 = nd1a; i1 <= nd1b; i1++)
                {
                    fprintf(matlabFile, "%s(%3d,%3d,%d)=%12.5e;", name, i1 - nd1a + 1, i2 - nd2a + 1, m + 1, u(i1, i2, m));
                    if (count % numPerLine == numPerLine - 1)
                        fprintf(matlabFile, "\n"); // newline
                    count++;
                }
        fprintf(matlabFile, "\n");
    }
    return 0;
}

#endif