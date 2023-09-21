#ifndef GET_CPU_H
#define GET_CPU_H

#include <ctime>
//--------------------------------------------
// Returnthecurrentwall-clocktimeinseconds
//--------------------------------------------
inline double getCPU()
{
    return (1.0 * std::clock()) / CLOCKS_PER_SEC;
}

#endif