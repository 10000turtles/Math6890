//
// Debugfiledemo:
// Openseparatedebugfilesoneachrank
// Examples:
// mpiexec-n2debugFileDemo-debug=1
// mpiexec-n4debugFileDemo-debug=1

// define thistoindicatewereareusing MPI
#define USE_PPP

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <float.h>
#include <assert.h>
// sleepfunctionisherewithLinux
#include <unistd.h>

// define anewtype"Real"whichisequivalenttoa"double "
typedef double Real;

#include <string>
using std::max;
using std::string;

// include commandstpparsecommandlinearguments
#include "parseCommand.h"

#include <ctime>
//--------------------------------------------
// Returnthecurrentwall-clocktimeinseconds
//--------------------------------------------
inline double getCPU()
{
#ifdef USE_PPP
    return MPI_Wtime(); // useMPItimer
#else
    return (1.0 * std::clock()) / CLOCKS_PER_SEC;
#endif
}

//======================================================================================
// Returnthemaxvalueofascalaroverallprocessorsinacommunicator
/// processor:return theresulttothisprocessor(-1equalsallprocessors)
//======================================================================================
Real getMaxValue(Real value, int processor = -1, MPI_Comm comm = MPI_COMM_WORLD)
{
    Real maxValue = value;
#ifdef USE_PPP
    if (processor == -1)
        MPI_Allreduce(&value, &maxValue, 1, MPI_DOUBLE, MPI_MAX, comm);
    else
        MPI_Reduce(&value, &maxValue, 1, MPI_DOUBLE, MPI_MAX, processor, comm);
#endif
    return maxValue;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv); // initializeMPI
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // myprocessnumber
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np); // totalnumberofproc’s

    if (myRank == 0)
        printf("Usage:debugFileDemo-debug=i\n");

    int debug = 1; // setto1fordebuginfo

    string line;
    bool echo = false; // donotechoinparseCommand
    for (int i = 1; i < argc; i++)
    {
        line = argv[i];
        if (parseCommand(line, "-debug=", debug, echo))
        {
        }
    }

    FILE *debugFile = NULL;
    if (debug > 0)
    {
        // openadebugfileoneachprocessor(inthedebugfolder)
        char debugFileName[80];
        sprintf(debugFileName, "debug/debugFileDemoNp%dProc%d.debug", np, myRank);
        debugFile = fopen(debugFileName, "w");
    }

    // WriteheaderinfotobothstdoutanddebugFile
    for (int ifile = 0; ifile <= 1; ifile++)
    {
        FILE *file = ifile == 0 ? stdout : debugFile;

        if ((ifile == 0 && myRank == 0) || // writetoterminalifmyRank==0
            (ifile == 1 && debug > 0))     // writetodebugFileifdebug>0
        {
            fprintf(file, "-------------------DebugFileDemo---------------------\n");
            fprintf(file, "np=%d,myRank=%d\n", np, myRank);
        }
    }

    Real cpu0 = getCPU();

    // Pretendtodosomecomputationsbysleepingforsometime

    // sleep(myRank*100000);//sleepforsometimeinseconds
    usleep(myRank * 100000); // sleepforsometimeinmicrosecs,1e-6(s)

    Real cpu = getCPU() - cpu0;

    if (debug > 0)
    {
        fprintf(debugFile, "Time to sleep=%9.2e(s)\n", cpu);
        fflush(debugFile); // flushoutputtofile
    }

    cpu = getMaxValue(cpu);
    if (myRank == 0)
        printf("Max time to sleep=%9.2e(s)\n", cpu);

    if (debugFile)
        fclose(debugFile);

    // closedownMPI
    MPI_Finalize();

    return 0;
}