
import numpy as np
import math
from math import e
from math import pi
from numpy import linalg as la
from scipy import linalg as LA
import matplotlib.pyplot as plt
import sys
from scipy import special

pol  = sys.argv[1]
A0  = float(sys.argv[2])
m_order = int(sys.argv[3])
m0_order = -1 * m_order


if pol =='x':
  A0x = A0
  A0y = 0.0
  A0z = 0.0
elif pol=='xy':
  A0x = A0/np.sqrt(2)
  A0y = A0/np.sqrt(2)
  A0z = 0.0
elif pol=='xyz':
  A0x = A0/np.sqrt(3)
  A0y = A0/np.sqrt(3)
  A0z = A0/np.sqrt(3)

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
        if i == line_begin_cell+3:
          a3 = np.array(line.split()).astype(float) 
          break

A_trans=np.array([a1,a2,a3])
A_inv=np.linalg.inv(A_trans.transpose())

r = np.zeros((num_wann,3))
X = np.loadtxt('orbital_coordinates.in')
for i in range(num_wann):
  r[i]=A_inv.dot(X[i])

hr_file = 'wannier_hr.dat'
with open(hr_file, "r") as f:
  f.readline()
  f.readline()
  num_cells=int(f.readline())
  
skip = int(num_cells/15)+4
#print(num_cells,skip)
data = np.loadtxt(hr_file,skiprows=skip)

lines = data.shape[0]
R = data[:,0:3]
m = data[:,3].astype(int)
n = data[:,4].astype(int)
H_R_real = data[:,5]
H_R_img = data[:,6]

mn = data[:,3:5].astype(int)

H_R_new = np.zeros((lines,2))
H_R_new0 = np.zeros((lines,2))
H_R_complex = np.zeros(lines, dtype = complex)
H_R_complex_new = np.zeros(lines, dtype = complex) # Create array H(m)
H_R_complex_new0 = np.zeros(lines, dtype = complex) # Create array H(-m)

for j in range( 0 , lines ):
    d1= R[j][0]+r[n[j]-1][0]-r[m[j]-1][0]
    d2= R[j][1]+r[n[j]-1][1]-r[m[j]-1][1]
    d3= R[j][2]+r[n[j]-1][2]-r[m[j]-1][2]
    dx = d1*a1[0] + d2*a2[0] + d3*a3[0]
    dy = d1*a1[1] + d2*a2[1] + d3*a3[1]
    dz = d1*a1[2] + d2*a2[2] + d3*a3[2]
    A_dot_d = A0x*dx + A0y*dy + A0z*dz
    H_R_complex[j] = H_R_real[j] + 1j * H_R_img[j] # Convert H_R from data to complex numbers
    H_R_complex_new[j] = H_R_complex_new[j] + 1j ** m_order * special.jv(m_order, A_dot_d) * H_R_complex[j] #* e ** (1j * m_order * omega * time)
    H_R_complex_new0[j] = H_R_complex_new0[j] + 1j ** m0_order * special.jv(m0_order, A_dot_d) * H_R_complex[j] #* e ** (1j * m_order * omega * time)
    #H_R_complex_new[j] = H_R_complex_new[j] + 1j ** m_order * H_R_complex[j]
    H_R_new[j][0] = H_R_complex_new[j].real
    H_R_new[j][1] = H_R_complex_new[j].imag
    H_R_new0[j][0] = H_R_complex_new0[j].real
    H_R_new0[j][1] = H_R_complex_new0[j].imag
    #H_R_new[j][0] = 1j ** m_order * special.jv(m_order, A_dot_d) * H_R_real[j]
    #H_R_new[j][1] = 1j ** m_order * special.jv(m_order, A_dot_d) * H_R_img[j]
    
'''
    if (R[j].all == 0) and (m[j] == n[j]):
        H_R_complex_new[j] = H_R_complex[j] # Only do 0th order and no more higher order Floquet calculations on E_site
    else:
        H_R_complex_new[j] = H_R_complex_new[j] + 1j ** m_order * special.jv(m_order, A_dot_d) * H_R_complex[j] #* e ** (1j * m_order * omega * time)
'''

new_data = np.hstack((R, mn, H_R_new))
new0_data = np.hstack((R, mn, H_R_new0))
np.savetxt('hrJ_pol_'+pol+'_A_'+'{:.2f}'.format(A0)+'_order_'+'{:d}'.format(m_order)+'.dat', new_data, fmt= 5*"%6d" + "%24.12e"*2)
np.savetxt('hrJ_pol_'+pol+'_A_'+'{:.2f}'.format(A0)+'_order_'+'{:d}'.format(m_order)+'_0.dat', new0_data, fmt= 5*"%6d" + "%24.12e"*2)
