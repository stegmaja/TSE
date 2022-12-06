#!/bin/bash --login
#SBATCH -n 32  
#SBATCH -o /scratch/c.c1990490/output.TSE2.0.%J #Job output
#SBATCH -t 12:00:00   
#SBATCH -p htc             
module purge
module load hpcw
module load parallel
module load anaconda/2020.02
parallel -N 1 -j $SLURM_NTASKS "python Stellartriples.py" ::: $(seq 0 50000)