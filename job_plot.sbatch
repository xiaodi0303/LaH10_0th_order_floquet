#!/bin/bash

#SBATCH -J   job           # Job name
#SBATCH -o   out.o%j       # Name of stdout output file
#SBATCH -e   out.e%j       # Name of stderr error file
#SBATCH -p development     # Queue (partition) name
#SBATCH -N 1               # Total # of nodes 
#SBATCH -n 56              # Total # of mpi tasks
#SBATCH -t 00:10:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --mail-user=username@tacc.utexas.edu
#SBATCH -A DMR21001        # Project/Allocation name (req'd if you have more than 1)

pol=x
m_order=2

mkdir figures
rm data.dat
rm data_max.dat

for z in $(seq 0.00 0.25 1.00)
do
  echo A0_$z
  cd A0_$z
  python plot_bands_dos.py $pol $z $m_order > output_plot &
  #python plot_pdos.py $pol $z &
  cd ..
done

wait

for z in $(seq 0.00 0.25 1.00)
do
  cd A0_$z
  echo $z  `grep rho_fermi  output_plot  | awk '{print$3}'` >> ../figures/data.txt
  echo $z  `grep rho_max  output_plot  | awk '{print$3}'` >> ../figures/data_max.txt
  cd ..
done

cd figures
mv ../A0_*/*png .
convert -resize 80% -delay 20 -loop 0 bands*.png bands_floquet.gif &
#convert -resize 80% -delay 20 -loop 0 pdos*.png pdos_floquet.gif &
#python  ../python_codes/plot_E_rho.py $pol




wait

# ---------------------------------------------------
