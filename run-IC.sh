#!/bin/bash --login
#SBATCH -n 1                  #Number of processors in our pool
#SBATCH -o /scratch/c.c1990490/output.%J #Job output
#SBATCH -t 8:00:00               #Max wall time for entire job
#SBATCH -d singleton

#change the partition to compute if running in Swansea
#SBATCH -p htc                 #Use the High Throughput partition which is intended for serial jobs
echo "$1 $2"
make initialise
./initialise $1 $2
