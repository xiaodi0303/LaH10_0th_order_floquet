
pol=x
m_order=12
omega=1.5 # (eV)

for z in $(seq 0.00 0.01 1.00)
do
    echo "A0 = $z"
    cd A0_$z
    #python  LaH10_write_hrJ.py $pol $z $m_order
    #sed "s/AA/A=$z/" ../python_codes/job_template2  | sed "s/CC/pol=$pol/" | sed "s/DD/m_order=$m_order/" | sed "s/EE/omega=$omega/" > job_bands.sbatch
    sed -i "s/^pol=.*/pol=$pol/" job_bands.sbatch
    sed -i "s/^m_order=.*/m_order=$m_order/" job_bands.sbatch
    sed -i "s/^omega=.*/omega=$omega/" job_bands.sbatch
    sbatch  job_bands.sbatch
    cd ..
done

