#!/bin/bash --login

for i in {1..10}
do
   if [ ! -d "./$i" ] 
   then
      cp -rp TSE "$i"
      cd "$i"
      sbatch --job-name="$i-IC" run-IC.sh 0.0002 1000000
      cd ..
   else
      cd "$i"
      sbatch --job-name="$i-1000000" run.sh 0.020 2 3 1
      sbatch --job-name="$i-1" run.sh 0.0002 2 3 1
      cd ..   
   fi
done
