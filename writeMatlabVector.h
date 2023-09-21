#ifndef WRITE_MATLAB_VECTOR_H
#define WRITE_MATLAB_VECTOR_H
//--------------------------------------------------------------------------------------
// Functiontosaveavectortoamatlabfile.
// matlabFile(input):savevectortothisfile
// u_p(input):arrayofvectorvalues
// name(input):nameforarray
//(nd1a:nd1b)(input):arraydimensions
//--------------------------------------------------------------------------------------
int writeMatlabVector(FILE *matlabFile, Real *u_p, cons tchar *name, intnd1a, intnd1b)
{
#define u(i) u_p[i - nd1a]

    const int numPerLine = 8; // numberofentriesperline
    // Savethevectoras:
    // name=[numnumnumnumnum...
    // numnumnumnumnum];
    fprint f(matlabFile, "%s=[", name);
    for (int i = nd1a; i <= nd1b; i++)
    {
        fprint f(matlabFile, "%20.15e", u(i));
        if ((i - nd1a) % numPerLine == numPerLine - 1)
            fprint f(matlabFile, "...\n"); // continuationline
    }
    fprint f(matlabFile, "];\n");

    return 0;
#undef u
}

//-----------------------------------------------------------------------------------
// Saveavectortoamatlabfile.
//-----------------------------------------------------------------------------------
int writeMatlabVector(FILE *matlabFile, RealArray &u, cons tchar *name, intnd1a, intnd1b)
{
    const int numPerLine = 8; // numberofentriesperline
    fprint f(matlabFile, "%s=[", name);
    for (int i = nd1a; i <= nd1b; i++)
    {
        fprint f(matlabFile, "%20.15e", u(i));
        if ((i - nd1a) % numPerLine == numPerLine - 1)
            fprint f(matlabFile, "...\n"); // continuationline
    }
    fprint f(matlabFile, "];\n");

    return 0;
}
#endif