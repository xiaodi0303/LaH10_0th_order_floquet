pol=x
m_order=1
omega=3

for z in $(seq 0.00 0.01 1.00)
do
    echo "A0 = $z"
    cd A0_$z
    #python  LaH10_write_hrJ.py $pol $z
    sed "s/AA/A=$z/" ../python_codes/job_template2  | sed "s/CC/pol=$pol/" | sed "s/DD/m_order=$m_order/" | sed "s/EE/omega=$omega/" > job_bands.sbatch
    sed -i "s/^#SBATCH -t.*/#SBATCH -t 01:00:00        # Run time (hh:mm:ss)/" job_bands.sbatch
    sed -i "/python LaH10_bands_floquet.py/i\python LaH10_write_hrJ.py \$pol \$A \$m_order" job_bands.sbatch
    sbatch  job_bands.sbatch
    cd ..
done

