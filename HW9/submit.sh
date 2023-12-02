#!/bin/bash -xe 
# Usage: 
# sbatch submit.sh 
# DRPcluster:2x8core 
# 
#SBATCH --job-name=heat2d #Jobname 
#SBATCH --ntasks=32 #TotalnumberofMPIranks 1
#SBATCH --nodes=4 #Numberofnodestobeallocated 
#SBATCH --time=00:05:00 #Walltimelimit(days-hrs:min:sec) 
#SBATCH --output=heat2dDemo%j.log #Paerelaorkingdirectory 
#SBATCH --error=heat2dDemo%j.err #Pathtrelativetotheworking directory 
echo "Date =$(date)" 
echo "Hostname =$(hostname-s)" 
echo "WorkingDirectory=$(pwd)"
echo ""
echo "Number ofNodesAllocated = $SLURM_JOB_NUM_NODES" 
echo "Number ofTasksAllocated = $SLURM_NTASKS" 

mpirun ./heat2d -nx=256 -tFinal=.01 -debug=1 -commOption=0