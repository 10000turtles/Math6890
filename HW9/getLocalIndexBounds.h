#ifndef GET_LOCAL_INDEX_BOUNDS_H
#define GET_LOCAL_INDEX_BOUNDS_H

//=======================================================================
// Distribute an index dimension across npx processors:
// rank(input):rank of process,0<=rank<npx
// npx(input):distribute across this many processors
// nx(input):total number of grid cells on global grid
// nx_l,n1a_l,n1b_l(output):localvalues
//
// Note:nx+1=totalgridpointsacrossallprocessors
//========================================================================
int getLocalIndexBounds(const int rank, const int npx, const int nx,
                        int &nx_l, int &n1a_l, int &n1b_l)
{
    //----Example:distributegridpoints-----
    // Global:
    // X--X--X--X--X--X--X--X--X--X--X
    // nx=10
    // n1an1b
    // Local:(rank=0,np=2)
    // X--X--X--X--X--X
    // nx_l=5
    // n1a_l n1b_l
    // Local:(rank=1,np=2)
    // X--X--X--X--X
    // nx_l=4
    // n1a_ln1b_l
    //
    nx_l = (nx + 1) / npx - 1; // nx_l+1=(nx_1+1)/np
    n1a_l = (nx_l + 1) * rank; // localstartingindex

    // Theremaybeextrapointsif(nx+1)isnotamultipleofnp
    int extra = nx + 1 - npx * (nx_l + 1);
    if (rank < extra)
    { // addoneextrapointtoproc’sontheleftside
        nx_l += 1;
        n1a_l += rank;
    }
    else
    {
        n1a_l += extra;
    }
    n1b_l = n1a_l + nx_l; // localendindex

    return 0;
}

#endif