#!/bin/bash
 
#SBATCH --mail-user=stegmaja@mpa-garching.mpg.de --mail-type=all --mem=20G --array=1-299 --nodes=1 --ntasks-per-node=40 --dependency=singleton -J TSE -t 24:00:00

module purge

module load anaconda/3/2023.03 parallel/201807

TASK_ID=${SLURM_ARRAY_TASK_ID}
TASK_MAX=${SLURM_ARRAY_TASK_MAX}

srun="srun -N1 -n1 --exclusive"

parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog parallel_joblog_$TASK_ID"

$parallel "$srun ./runtask {1} {2}" ::: $(seq $TASK_ID $TASK_MAX 10000) ::: 0.0002 0.002 0.02

sbatch $0