
import os
import numpy as np
import matplotlib.pyplot as plt
import sys

pol = sys.argv[1]
A0  = float(sys.argv[2])
m_order = int(sys.argv[3])
m_max = m_order - 1

num_e = 19
totalsteps = 35001

fermi = 20.4331 # Manually input
Emin = -10
Emax = 6

'''
omega = 500

E_axis = np.loadtxt('dos_E_axis.dat')
dim=E_axis.shape

spectra=np.zeros(dim)
sum_spectra=np.zeros(dim)

dir_list=os.listdir()
count = 0
for i in dir_list:
  if "dos_spectra.dat" in i:
    tmp = np.loadtxt(i)
    spectra += tmp
    count += 1

print("num_spectra = ", count)

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

np.savetxt('Dos_E_axis.dat', E_axis)
np.savetxt('Dos_spectra.dat', spectra)
np.savetxt('Dos_sum_spectra.dat', sum_spectra)
'''
k_axis = np.loadtxt('bands_k_axis.dat')
eigenvals_data = np.loadtxt('bands_eigenvals.dat')

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

high_points=[]
hp_labels=[]

for i in range(len(line)):
  if line[i][0] == "G":
    hp_labels.append("Γ")
  else:
    hp_labels.append(line[i][0])
  high_points.append(float(line[i][-7:]))

f.close()

fig = plt.figure(figsize=(7.5,4.5), constrained_layout=True)
#gs = fig.add_gridspec(1, 4)
ax1 = fig.add_subplot()
#ax1 = fig.add_subplot(gs[:,0:3])
#ax2 = fig.add_subplot(gs[:,3])

for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(2)
for i in range(num_wann * (2 * m_max + 1)):
    y = eigenvals_data[:,i] - fermi
    ax1.plot(k_axis,y,c='b',linewidth=2)

ax1.axhline(y=0.0, color='gray', linestyle='--', linewidth=1)

if len(high_points)>=3:
    for i in range(len(high_points)):
        ax1.axvline(high_points[i],c='k',linewidth=1.5)

ax1.set_xticks(high_points)#,fontsize=fontsize)
ax1.set_xticklabels(hp_labels)#,fontsize=fontsize)
ax1.tick_params(labelsize=14)

ax1.set_xlabel(r'$k$ vector',fontsize=16)
ax1.set_ylabel('energy (eV)', fontsize=16)
ax1.set_xlim(high_points[0],high_points[-1])

ax1.set_ylim(Emin,Emax)
#################################################################

'''
ax2.plot(spectra, E_axis-fermi, c= 'b',linewidth=2)
ax2.axhline(y=0.0, color='gray', linestyle='--', linewidth=1)
#ax2.set_ylabel('energy (eV)', fontsize=16)
ax2.set_xlabel('DOS',fontsize=16)

ax2.text(2.0, -1, r'$\rho(E_F)$ = '+'{:5.3f}'.format(rho_fermi), fontsize=12)
ax2.text(2.0, -8.5, "polarization: "+pol, fontsize=10)

E0=int(A0*omega)
ax2.text(2.0, -9.5, "$E_0$ = " + str(E0) +" mV/Å", fontsize=10)

ax2.set_ylim(Emin,Emax)
ax2.set_xlim(0,10)

ax2.set_yticks([])

for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(2)

ax2.tick_params(labelsize=14)
'''

#plt.suptitle("polarization: x\n $A_0$ = " + A0)
plt.savefig('bands_dos_'+"{:d}".format(m_order)+'_order_'+"{:.2f}".format(A0)+'.png',dpi=300,bbox_inches='tight')
