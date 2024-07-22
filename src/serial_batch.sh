#!/bin/bash
#SBATCH --mail-user=stegmaja@mpa-garching.mpg.de --mail-type=all --mem=20G --array=1-100 --nodes=1 --ntasks-per-node=40 --dependency=singleton -J TSE -t 24:00:00

module purge
module load anaconda/3/2023.03 parallel/201807

TASK_ID=${SLURM_ARRAY_TASK_ID}
TASK_MAX=${SLURM_ARRAY_TASK_MAX}

# Set number of OMP threads to fit the number of available cpus, if applicable.
export OMP_NUM_THREADS=1
parallel --delay .1 -j 40 --resume-failed --joblog parallel_joblog_$TASK_ID ./runtask {1} {2} ::: $(seq $TASK_ID $TASK_MAX 50000) ::: 0.03 0.02 0.01 0.0075 0.005 0.0025 0.001 0.00075 0.0005 0.00025 0.0001

sbatch $0