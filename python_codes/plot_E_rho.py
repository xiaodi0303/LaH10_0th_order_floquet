
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
#import numpy.polyfit

pol = sys.argv[1]
omega = 500

data=np.loadtxt('data.txt')
data_max=np.loadtxt('data_max.txt')

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

ax.plot(data[:,0]*omega, data[:,1],'-ob', alpha=0.6,markersize=10 , label=r'$\rho_F$')
ax.plot(data_max[:,0]*omega, data_max[:,1],'-sr', alpha=0.6,markersize=10 , label=r'$\rho_{max}$')


plt.xlabel('electric field (mV/Å)',fontsize=14)
plt.ylabel('density of states (states/eV per unitcell)',fontsize=14)
plt.xticks(fontsize=14)#, rotation=90)
plt.yticks(fontsize=14)#, rotation=90)
plt.legend(fontsize=14)

#plt.axes().set_aspect(0.4/3.5)

ax = plt.gca()
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
ax.spines['right'].set_linewidth(1.5)
ax.spines['top'].set_linewidth(1.5)
plt.tick_params(axis='both', labelsize=12)
#plt.tick_params(axis='both', width=2, length=5, labelsize=12)

#Locator = MultipleLocator(2)
#ax.xaxis.set_major_locator(Locator)
ax.set_xlim(0,omega)
ax.set_ylim(0.8,1.7)

plt.title("polarization: "+pol)
plt.savefig('rho_E.png',dpi=300, bbox_inches='tight')
#plt.show()
plt.clf()



