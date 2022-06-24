set style data dots
set nokey
set xrange [0: 4.34612]
set yrange [ -3.91302 : 47.72245]
set arrow from  0.89683,  -3.91302 to  0.89683,  47.72245 nohead
set arrow from  1.99523,  -3.91302 to  1.99523,  47.72245 nohead
set arrow from  3.26354,  -3.91302 to  3.26354,  47.72245 nohead
set arrow from  3.89770,  -3.91302 to  3.89770,  47.72245 nohead
set xtics ("W"  0.00000,"L"  0.89683,"GAMMA"  1.99523,"X"  3.26354,"W"  3.89770,"K"  4.34612)
 plot "wannier_band.dat"
