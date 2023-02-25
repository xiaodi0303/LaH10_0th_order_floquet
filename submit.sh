
pol=x
m_order=2
omega=3 # (eV)

sed -i "s/CC/pol=$pol/" job_plot.sbatch

cd wannier_files
grep -A200 "Final State" wannier.wout | grep "WF centre and spread" | awk '{print$7"  "$8"  "$9}' |  sed 's/,//g' > orbital_coordinates.in 
cd ..

for z in $(seq 0.00 0.25 1.00)
do
  echo "A0 = $z"
  mkdir A0_$z ; cd A0_$z
  ln -s ../python_codes/LaH10_bands_floquet.py .
  ln -s ../python_codes/LaH10_pdos.py .
  ln -s ../python_codes/LaH10_write_hrJ.py .
  ln -s ../python_codes/plot* .
  ln -s ../k-grid_files/kpoints* .
  ln -s ../wannier_files/* .
  python  LaH10_write_hrJ.py $pol $z $m_order
  
  #Below must have folder hrJ in main(hrJ_$m_order copied as "hrJ" from job_collect_hrJ.sbatch executed and copied from a lower order Floquet folder)
  cp ../hrJ/hrJ_pol_${pol}_A_${z}* .
  #Below must have folder npy in main(npy_$m_order copied as "npy" from job_collect_npy.sbatch executed and copied from a lower order Floquet folder)
  cp ../npy/hf_block_A_${z}* .

  for i in {0..7}
  do
    sed "s/AA/A=$z/" ../python_codes/job_template1  | sed "s/BB/$i/" | sed "s/CC/pol=$pol/" > job_${i}.sbatch
    #sbatch  job_${i}.sbatch
  done
  sed "s/AA/A=$z/" ../python_codes/job_template2  | sed "s/CC/pol=$pol/" | sed "s/DD/m_order=$m_order/" | sed "s/EE/omega=$omega/" > job_bands.sbatch
  sbatch  job_bands.sbatch

  cd ..
done

