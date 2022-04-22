#!/bin/bash --login
#SBATCH -n 8                  #Number of processors in our pool
#SBATCH -o /scratch/c.c1990490/output.%J #Job output
#SBATCH -t 72:00:00               #Max wall time for entire job
#SBATCH -d singleton

#change the partition to compute if running in Swansea
#SBATCH -p htc                 #Use the High Throughput partition which is intended for serial jobs

# Parallelization for local machine
if [ -z "$SLURM_NTASKS" ]
then
      echo "Local Machine"
      SLURM_NTASKS=15
fi

echo "$(pwd)"

module load hpcw
module load parallel

initialise () {
  mkdir "$1-$2-$3-$4"
  touch "$1-$2-$3-$4/Massless-remnant.dat"
  touch "$1-$2-$3-$4/broken.dat"
  touch "$1-$2-$3-$4/collision.dat"
  touch "$1-$2-$3-$4/breaks.dat"
  touch "$1-$2-$3-$4/unstable.dat"
  touch "$1-$2-$3-$4/disrupt.dat"
  touch "$1-$2-$3-$4/3-RLO.dat"
  touch "$1-$2-$3-$4/CO-triple.dat"
  touch "$1-$2-$3-$4/CO-binary.dat"
  touch "$1-$2-$3-$4/no-CO-binary.dat"
  touch "$1-$2-$3-$4/CO-binary-alone.dat"
  touch "$1-$2-$3-$4/no-CO-binary-alone.dat"
  touch "$1-$2-$3-$4/merger.dat"
  touch "$1-$2-$3-$4/parallel_joblog"

  i=$(wc -l < "stable-IC.dat")

  parallel -N 1 -j $SLURM_NTASKS --resume --joblog "$1-$2-$3-$4/parallel_joblog" "./LKspin-s $1 $2 $3 $4 {1}" ::: $(seq 1 $i)
}

initialise $1 $2 $3 $4
