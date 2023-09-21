#ifndef WRITE_MATLAB_VECTOR_H
#define WRITE_MATLAB_VECTOR_H
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
#undef u
}

//-----------------------------------------------------------------------------------
// Saveavectortoamatlabfile.
//-----------------------------------------------------------------------------------
int writeMatlabVector(FILE *matlabFile, RealArray &u, cons tchar *name, intnd1a, intnd1b)
{
    const int numPerLine = 8; // numberofentriesperline
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
#endif