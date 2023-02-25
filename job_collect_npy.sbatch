#!/bin/bash

#SBATCH -J   job           # Job name
#SBATCH -o   out.o%j       # Name of stdout output file
#SBATCH -e   out.e%j       # Name of stderr error file
#SBATCH -p development     # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 56              # Total # of mpi tasks
#SBATCH -t 00:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --mail-user=username@tacc.utexas.edu
#SBATCH -A DMR21001        # Project/Allocation name (req'd if you have more than 1)

pol=x

rm -r npy
mkdir npy
for z in $(seq 0.00 0.25 1.00)
do
  echo A0_$z
  cd npy
  cp ../A0_$z/*.npy .
  cd ..
done

wait

# ---------------------------------------------------
