
import os
import numpy as np
import matplotlib.pyplot as plt
import sys

pol = sys.argv[1]
A0  = float(sys.argv[2])

num_e = 19
totalsteps = 35001

Emin = -2
Emax = 2

omega = 500
label_list=[r'1×$La:p$',r'1×$La:d$',r'1×$La:f$',r'8×$H^1:s$',r'2×$H^2:s$']

E_axis = np.loadtxt('dos_E_axis.dat')
dim=E_axis.shape[0]

spectra=np.zeros(dim)
sum_spectra=np.zeros(dim)
pspectra=np.zeros((5,dim))

dir_list=os.listdir()
count = 0
for i in dir_list:
  if "dos_spectra.dat" in i:
    tmp = np.loadtxt(i)
    spectra += tmp
    count += 1
print("num_spectra = ", count)

for k in range(5):
  count = 0
  for i in dir_list:
    if "dos_spectra_"+str(k)+".dat" in i:
      tmp = np.loadtxt(i)
      pspectra[k] += tmp
      count += 1
  print("num_pspectra_" + str(k) + " = ", count)


dE = E_axis[1]-E_axis[0]

temp=0.0
mark = -1
for ii in range(0,totalsteps):
    temp += spectra[ii]*dE
    sum_spectra[ii] = temp
    if (mark < 0) and (sum_spectra[ii] > num_e):
        mark *= -1
        EF_index = ii-1

print("EF_index = ", EF_index)

fermi = E_axis[EF_index] + (num_e-sum_spectra[EF_index])/(sum_spectra[EF_index+1]-sum_spectra[EF_index])* dE 
rho_fermi = spectra[EF_index] +  (spectra[EF_index+1]-spectra[EF_index]) / (sum_spectra[EF_index+1]-sum_spectra[EF_index]) * (num_e-sum_spectra[EF_index])
print('fermi = ', fermi)
print('rho_fermi = ', rho_fermi)
rho_max = np.max(spectra[EF_index-50:EF_index+50])
print('rho_max = ', rho_max)
rho_min = np.min(spectra[EF_index-50:EF_index+50])
print('rho_min = ', rho_min)

win_file= 'wannier.win'
searchfile = open(win_file, "r")
for i, line in enumerate(searchfile):
    if "num_wann" in line: num_wann = int(line.split()[-1])
searchfile.close()

f = open('wannier_band.gnu','r')                                        
for l in f:
  if 'set xtics' in l:
    line=l[12:-2].split(',"')
    #print(line)


fig = plt.figure(figsize=(7.5,4.5), constrained_layout=True)

ax = fig.add_subplot(111)

ax.plot(E_axis-fermi, spectra, c= 'k',linewidth=2, label='total')
for a in range(5):
  ax.plot(E_axis-fermi, pspectra[a], linewidth=1.5, label=label_list[a])
ax.legend(prop=dict(size=12))

ax.axvline(x=0.0, color='gray', linestyle='--', linewidth=1)
ax.set_xlabel('energy (eV)', fontsize=16)
ax.set_ylabel('DOS (states/eV per unitcell)',fontsize=16)

ax.text(-0.75,1.5, r'$\rho(E_F)$ = '+'{:5.3f}'.format(rho_fermi), fontsize=12)
ax.text(-1.8,1.8, "polarization: "+pol, fontsize=12)

E0=int(float(A0)*omega)
ax.text(-1.8,1.65, "$E_0$ = " + str(E0) +" mV/Å", fontsize=12)

ax.set_xlim(Emin,Emax)
ax.set_ylim(0,2)

#ax.set_yticks([])

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)

ax.tick_params(labelsize=14)


#plt.suptitle("polarization: x\n $A_0$ = " + A0)
plt.savefig('pdos_'+"{:.2f}".format(A0)+'.png',dpi=300,bbox_inches='tight')
