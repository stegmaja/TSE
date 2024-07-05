#!/bin/bash
 
#SBATCH --mail-user=stegmaja@mpa-garching.mpg.de --mail-type=all --mem=0 --array=1-299 --nodes=1 --ntasks-per-node=40 --dependency=singleton -J TSE -t 24:00:00

module purge

module load anaconda/3/2023.03

srun="srun -N1 -n1 --exclusive"

parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog parallel_joblog_$TASK_ID"

$parallel "$srun ./runtask {1} {2}" ::: {1..10000} ::: 0.0002 0.002 0.02

sbatch $0