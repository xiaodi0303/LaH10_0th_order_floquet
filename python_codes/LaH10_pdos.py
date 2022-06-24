
import numpy as np
import math
from math import e
from math import pi
from numpy import linalg as la
from scipy import linalg as LA
import matplotlib.pyplot as plt
import sys
from scipy import special

kk = 75
file_index1 = int(sys.argv[3])
file_index2 = int(sys.argv[4])

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
searchfile.close()

with open(win_file, "r") as f:
    for i, line in enumerate(f):
        if i == line_begin_cell+1: a1 = np.array(line.split()).astype(float)
        if i == line_begin_cell+2: a2 = np.array(line.split()).astype(float)
        if i == line_begin_cell+3: a3 = np.array(line.split()).astype(float)

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


k_grid = np.loadtxt('kpoints.dat_'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2))
num_k = k_grid.shape[0]
#print(num_k)


hr_file = 'hrJ_A_'+"{:.2f}".format(A0)+'.dat'
hr_file = 'hrJ_pol_'+pol+'_A_'+"{:.2f}".format(A0)+'.dat'
data = np.loadtxt(hr_file)
lines = data.shape[0]
R = data[:,0:3]
m = data[:,3].astype(int)
n = data[:,4].astype(int)
H_R = data[:,5] + 1j * data[:,6]

H_k=np.zeros((num_wann , num_wann), dtype=complex)

m_order=0
phi=0
eigenvals_data= np.zeros(num_k*num_wann)
vecs_weight_data= np.zeros((num_k*num_wann,5))
for i in range(num_k):
    #print(i)
    H_k=np.zeros((num_wann , num_wann), dtype=complex)
    for j in range( 0 , lines ):
        d1= R[j][0]+r[n[j]-1][0]-r[m[j]-1][0]
        d2= R[j][1]+r[n[j]-1][1]-r[m[j]-1][1]
        d3= R[j][2]+r[n[j]-1][2]-r[m[j]-1][2]
        H_k[m[j]-1][n[j]-1] += H_R[j]* e ** (1j*2*pi* (k_grid[i][0]*d1 + k_grid[i][1]*d2 + k_grid[i][2]*d3) )
    H_k = H_k * e **(-1j*m_order* (phi+pi/2) )
    #eigenvals = la.eigvals(H_k).real
    #eigenvals = la.eigvalsh(H_k).real

    #eigenvals = LA.eigvalsh(H_k).real
    #print(eigenvals)
    #eigenvals_data[i*num_wann:(i+1)*num_wann] = eigenvals

    #print(np.allclose(H_k, H_k.conj().T, rtol=1e-8, atol=1e-8))
    val, vec = LA.eigh(H_k)
    for a in range(num_wann):
        #tmp=vec[0:25,a]
        #tmp=np.vdot(tmp,tmp)
        #print('mark', a, np.absolute(tmp))
        tmp=vec[15:23,a]
        tmp=np.vdot(tmp,tmp)
        w_H1_s=np.absolute(tmp)
        tmp=vec[23:25,a]
        tmp=np.vdot(tmp,tmp)
        w_H2_s=np.absolute(tmp)
        tmp=vec[0 :3 ,a]
        tmp=np.vdot(tmp,tmp)
        w_La_p=np.absolute(tmp)
        tmp=vec[3 :8 ,a]
        tmp=np.vdot(tmp,tmp)
        w_La_d=np.absolute(tmp)
        tmp=vec[8 :15,a]
        tmp=np.vdot(tmp,tmp)
        w_La_f=np.absolute(tmp)
    vector=np.array([w_La_p, w_La_d, w_La_f, w_H1_s, w_H2_s])
    #print(w_La_p, w_La_d, w_La_f, w_H1_s, w_H2_s)
    #print(w_La_p+ w_La_d+ w_La_f+ w_H1_s+ w_H2_s)
    #print(eigen[1])
    eigenvals_data[i*num_wann:(i+1)*num_wann] = val
    vecs_weight_data[i*num_wann:(i+1)*num_wann] = vector


eigenvals_data = eigenvals_data.ravel()

np.savetxt('dos_eigenvals.dat'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2), eigenvals_data)
np.savetxt('dos_vecs_weight.dat'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2), vecs_weight_data)

'''
import sys

kk = 75
file_index1 = int(sys.argv[2])
file_index2 = int(sys.argv[3])

A0  = float(sys.argv[1])
'''

################################################################


import numpy as np
import matplotlib.pyplot as plt
import math

def lz(x,x0,a):
    rtemp1 = 1/math.pi*a/((x-x0)**2+a**2)
    return rtemp1
def gs(x,x0,a):
    rtemp1 = 1/np.sqrt(2*np.pi*a**2) * np.exp( -(x-x0)**2/(2*a**2) )
    return rtemp1

vecs_weight_data = np.loadtxt('dos_vecs_weight.dat'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2))
filename = 'dos_eigenvals.dat'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2)

data = np.loadtxt(filename)
length_data = len(data)

broadening = 0.1
L_energy = -5
R_energy = 30
totalsteps = 35001

num_e = 19

win_file= 'wannier.win'
searchfile = open(win_file, "r")
for i, line in enumerate(searchfile):
    if "num_wann" in line: num_wann = int(line.split()[-1])
searchfile.close()


spectra_all = np.zeros(totalsteps)
spectra_0 = np.zeros(totalsteps)
spectra_1 = np.zeros(totalsteps)
spectra_2 = np.zeros(totalsteps)
spectra_3 = np.zeros(totalsteps)
spectra_4 = np.zeros(totalsteps)
#spectra_5 = np.zeros(totalsteps)
E_axis = np.zeros(totalsteps)

dE = (R_energy-L_energy)/(totalsteps-1)
#wk = np.loadtxt("wk.dat"+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2))
vec_w = vecs_weight_data
wk_gamma = 2/(kk*kk*kk) 
#tmp = wk/wk_gamma
#tmp=np.rint(tmp)
#wk_new = 2*tmp/(kk*kk*kk)

for ii in range(0,totalsteps):
    energy = L_energy + dE*ii
    E_axis[ii] = energy
    for jj in range(0,length_data):
        #spectra[ii] = spectra[ii] + lz(energy,data[jj],broadening)*wk_new[int(jj/num_wann)]
        tmp = gs(energy,data[jj],broadening)*wk_gamma
        spectra_all[ii] = spectra_all[ii] + tmp 
        spectra_0[ii] = spectra_0[ii] + tmp*vec_w[int(jj/num_wann)][0]
        spectra_1[ii] = spectra_1[ii] + tmp*vec_w[int(jj/num_wann)][1]
        spectra_2[ii] = spectra_2[ii] + tmp*vec_w[int(jj/num_wann)][2]
        spectra_3[ii] = spectra_3[ii] + tmp*vec_w[int(jj/num_wann)][3]
        spectra_4[ii] = spectra_4[ii] + tmp*vec_w[int(jj/num_wann)][4]
#        spectra_5[ii] = spectra_5[ii] + tmp*vec_w[int(jj/num_wann)][5]


if file_index1 == 0:
  np.savetxt('dos_E_axis.dat', E_axis)

np.savetxt('dos_spectra.dat'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2), spectra_all)
np.savetxt('dos_spectra_0.dat'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2), spectra_0)
np.savetxt('dos_spectra_1.dat'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2), spectra_1)
np.savetxt('dos_spectra_2.dat'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2), spectra_2)
np.savetxt('dos_spectra_3.dat'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2), spectra_3)
np.savetxt('dos_spectra_4.dat'+"{:02d}".format(file_index1)+'_'+"{:02d}".format(file_index2), spectra_4)
#np.savetxt('dos_sum_spectra.dat', sum_spectra)


