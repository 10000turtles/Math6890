#ifndef TRIDIAGONAL_H
#define TRIDIAGONAL_H "tridiagonal.h"

//
// Routinesforsolvingtridiagonallinearsystems
//
#include "A++.h"

typedef double Real;
typedef doubleSerialArray RealArray;

int factorTridiagonalMatrix(RealArray &Ax);

int solveTridiagonal(RealArray &Ax, RealArray &rhs);

#endif