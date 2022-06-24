
import numpy as np
import math
from math import e
from math import pi
from numpy import linalg as la
from scipy import linalg as LA
import matplotlib.pyplot as plt
import sys
from scipy import special

pol = sys.argv[1]
A0  = float(sys.argv[2])
A0x = A0
A0y = 0.0
A0z = 0.0


win_file= 'wannier.win'
searchfile = open(win_file, "r")
for i, line in enumerate(searchfile):
    if "num_wann" in line: num_wann = int(line.split()[-1])
    if "begin unit_cell_cart" in line:  line_begin_cell = i
    if "begin kpoint_path" in line:  line_begin_k = i
    if "end kpoint_path" in line:  line_end_k = i
    if "bands_num_points" in line: k_points_num_1st_path = float(line.split()[-1])
searchfile.close()

high_symm_k=[]
with open(win_file, "r") as f:
    for i, line in enumerate(f):
        if i == line_begin_cell+1: a1 = np.array(line.split()).astype(np.float)
        if i == line_begin_cell+2: a2 = np.array(line.split()).astype(np.float)
        if i == line_begin_cell+3: a3 = np.array(line.split()).astype(np.float)
        if i  >  line_begin_k and i < line_end_k: high_symm_k.append(line.split())
        if i == line_end_k: break

A_trans=np.array([a1,a2,a3])
A_inv=np.linalg.inv(A_trans.transpose())

r = np.zeros((num_wann,3))
X = np.loadtxt('orbital_coordinates.in')
for i in range(num_wann):
  r[i]=A_inv.dot(X[i])


volume = a1[0]*(a2[1]*a3[2]-a2[2]*a3[1]) + a1[1]*(a2[2]*a3[0]-a2[0]*a3[2]) + a1[2]*(a2[0]*a3[1]-a2[1]*a3[0])
b1= np.array([2*pi*(a2[1]*a3[2]-a2[2]*a3[1])/volume , 2*pi*(a2[2]*a3[0]-a2[0]*a3[2])/volume , 2*pi*(a2[0]*a3[1]-a2[1]*a3[0])/volume])
b2= np.array([2*pi*(a3[1]*a1[2]-a3[2]*a1[1])/volume , 2*pi*(a3[2]*a1[0]-a3[0]*a1[2])/volume , 2*pi*(a3[0]*a1[1]-a3[1]*a1[0])/volume])
b3= np.array([2*pi*(a1[1]*a2[2]-a1[2]*a2[1])/volume , 2*pi*(a1[2]*a2[0]-a1[0]*a2[2])/volume , 2*pi*(a1[0]*a2[1]-a1[1]*a2[0])/volume])

high_symm_k = np.array(high_symm_k)
num_of_paths = high_symm_k.shape[0]

path= np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
for i in range(num_of_paths):
    temp = np.array([])
    temp = np.append(temp,[high_symm_k[i][1:4].astype(np.float)])
    temp = np.append(temp,[high_symm_k[i][5:8].astype(np.float)])
    path = np.vstack( [path , temp] )
path = np.delete(path, (0), axis=0)

num_of_paths = high_symm_k.shape[0]
temp1 = b1[0]*(path[0][3]-path[0][0])+b2[0]*(path[0][4]-path[0][1])+b3[0]*(path[0][5]-path[0][2])
temp2 = b1[1]*(path[0][3]-path[0][0])+b2[1]*(path[0][4]-path[0][1])+b3[1]*(path[0][5]-path[0][2])
temp3 = b1[2]*(path[0][3]-path[0][0])+b2[2]*(path[0][4]-path[0][1])+b3[2]*(path[0][5]-path[0][2])
delta_k = math.sqrt(temp1 ** 2 + temp2 ** 2 + temp3 ** 2 ) / k_points_num_1st_path

#len_b1= math.sqrt(b1[0] ** 2 + b1[1] ** 2 + b1[2] ** 2)
#len_b2= math.sqrt(b2[0] ** 2 + b2[1] ** 2 + b2[2] ** 2)
#len_b3= math.sqrt(b3[0] ** 2 + b3[1] ** 2 + b3[2] ** 2)

print(path)
print(b1)
print(delta_k)

test = np.linspace(0.1,0.5,5)


k1= np.array([])
k2= np.array([])
k3= np.array([])

tot_num_of_k = 0
k_axis = np.array([])

kpt_file= 'wannier_band.kpt'
searchfile = open(kpt_file, "r")
tot_num_of_k = int(searchfile.readline())

tempx_old = 0.0
tempy_old = 0.0
tempz_old = 0.0
k_accu = 0.0

for i, line in enumerate(searchfile):
    n1 = float(line.split()[0])
    n2 = float(line.split()[1])
    n3 = float(line.split()[2])
    k1 = np.append(k1,n1)
    k2 = np.append(k2,n2)
    k3 = np.append(k3,n3)
    tempx = n1*b1[0]+ n2*b2[0] + n3*b3[0]
    tempy = n1*b1[1]+ n2*b2[1] + n3*b3[1]
    tempz = n1*b1[2]+ n2*b2[2] + n3*b3[2]
    k_accu = k_accu + math.sqrt((tempx-tempx_old) ** 2 + (tempy-tempy_old) ** 2 + (tempz-tempz_old) ** 2 )
    k_axis = np.append(k_axis,k_accu)
    tempx_old = tempx
    tempy_old = tempy
    tempz_old = tempz

searchfile.close()
k_axis= k_axis - k_axis[0]
print(k_axis)

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

hr_file = 'hrJ_pol_'+pol+'_A_'+"{:.2f}".format(A0)+'.dat'
data = np.loadtxt(hr_file)
lines = data.shape[0]
R = data[:,0:3]
m = data[:,3].astype(int)
n = data[:,4].astype(int)
H_R = data[:,5] + 1j * data[:,6]

H_k=np.zeros((num_wann , num_wann), dtype=np.complex)

m_order=0
phi=0
eigenvals_data= np.zeros((tot_num_of_k,num_wann))
for i in range(tot_num_of_k):
    print(i)
    H_k=np.zeros((num_wann , num_wann), dtype=np.complex)
    for j in range( 0 , lines ):
        d1= R[j][0]+r[n[j]-1][0]-r[m[j]-1][0]
        d2= R[j][1]+r[n[j]-1][1]-r[m[j]-1][1]
        d3= R[j][2]+r[n[j]-1][2]-r[m[j]-1][2]
        H_k[m[j]-1][n[j]-1] += H_R[j]* e ** (1j*2*pi* (k1[i]*d1 +k2[i]*d2 +k3[i]*d3) )
    H_k = H_k * e **(-1j*m_order* (phi+pi/2) )
    #eigenvals = la.eigvals(H_k).real
    #eigenvals = la.eigvalsh(H_k).real
    eigenvals = LA.eigvalsh(H_k).real
    #eigenvals = np.sort(eigenvals)
    eigenvals_data[i] = eigenvals

np.savetxt('bands_eigenvals.dat', eigenvals_data)
np.savetxt('bands_k_axis.dat', k_axis)

'''
fig = plt.figure()
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
for i in range(num_wann):
    y = eigenvals_data[:,i-1] - fermi
    ax.plot(k_axis,y,c='b',linewidth=2)




if len(high_points)>=3:
    for i in range(len(high_points)):
        ax.axvline(high_points[i],c='k',linewidth=1.5)


ax.set_xticks(high_points)#,fontsize=fontsize)
ax.set_xticklabels(hp_labels)#,fontsize=fontsize)
ax.tick_params(labelsize=14)

ax.set_ylabel('Energy (eV)', fontsize=18)
ax.set_xlim(high_points[0],high_points[-1])

ax.set_ylim(Emin,Emax)

plt.savefig("bands_0th_order.png",dpi=300,bbox_inches='tight')
'''
