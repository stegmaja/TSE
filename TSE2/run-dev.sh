#!/bin/bash --login
#SBATCH -n 64  
#SBATCH -o /scratch/c.c1718684/output.%J #Job output
#SBATCH -t 1:00:00   
#SBATCH -p dev             
module purge
module load hpcw
module load parallel
module load anaconda/2020.02
parallel -N 1 -j $SLURM_NTASKS "python Stellartriples.py" ::: $(seq 0 50000)