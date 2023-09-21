
//
#include "A++.h"

typedef double Real;
typedef doubleSerialArray RealArray;

// define macrosforarrayreferences
#define AX(i, j) Ax_p[(i) + 3 * (j)]
#define RHS(i) rhs_p[(i)]

//=====================================================================================
// Factoratridiagonalmatrix--thefactorizationisstoredinAx
//=====================================================================================
int factorTridiagonalMatrix(RealArray &Ax)
{
    const int iax = Ax.getBase(1), ibx = Ax.getBound(1);

    Real *Ax_p = Ax.getDataPoint er();

    // Ax.display("Axbeforefactor");

    // Factor:(nopivoting)
    //
    //[b0c0]
    //[a1b1c1]
    //[a2b2c2]
    //[a3b3c3]
    //[....]
    //[ambmcm]
    //[anbn]
    for (int i1 = iax + 1; i1 <= ibx; i1++)
    {
        Reald = -AX(0, i1) / AX(1, i1 - 1); //-a[i1]/b[i1-1]
        AX(1, i1) += d * AX(2, i1 - 1);
        AX(0, i1) = d; // savedhere
    }

    // Ax.display("Axafterfactor");
    return 0;
}

//=====================================================================================
// SolvethetridiagonalmatrixproblemgiventhefactoredmatrixAx
//
//=====================================================================================
int solveTridiagonal(RealArray &Ax, RealArray &rhs)
{
    const int iax = Ax.getBase(1), ibx = Ax.getBound(1);

    const Real *Ax_p = Ax.getDataPoint er();
    Real *rhs_p = rhs.getDataPoint er();

    //---forwardelimination--
    for (int i1 = iax + 1; i1 <= ibx; i1++)
    {
        RHS(i1) += AX(0, i1) * RHS(i1 - 1);
    }

    //---back-substitution--
    //[b0c0][x][x0]
    //[b1c1][x][x1]
    //[b2c2][x][x2]
    //[b3c3][x]=[x3]
    //[...][x][]
    //[bmcm][xm][rm]
    //[cn][xn][rn]
    RHS(ibx) = RHS(ibx) / AX(1, ibx);
    for (int i1 = ibx - 1; i1 >= iax; i1--)
    {
        RHS(i1) = (RHS(i1) - AX(2, i1) * RHS(i1 + 1)) / AX(1, i1);
    }

    return 0;
}