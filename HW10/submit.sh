#!/bin/bash -xe 
# Usage: 
# sbatch submit.sh 
# 
#SBATCH --job-name=addVector #Jobname 
#SBATCH --gres=gpu:1 #numberofgpuspernode 
#SBATCH --ntasks=1 #totalnumberofMPIranks 
#SBATCH --nodes=1 #Maximumnumberofnodestobeallocated 
#SBATCH --time=00:01:00 #Walltimelimit(days-hrs:min:sec) 
#SBATCH --output=addVector%j.log 
#SBATCH --error=addVector%j.err 

echo "Date =$(date)" 
echo "Hostname =$(hostname-s)" 
echo "WorkingDirectory=$(pwd)" 
echo "" 
echo "Number ofNodesAllocated = $SLURM_JOB_NUM_NODES" 
echo "Number ofTasksAllocated = $SLURM_NTASKS" 
./addVector -n=10000 
./addVector -n=1000000 
./addVector -n=100000000