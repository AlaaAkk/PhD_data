import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.transforms as tr
from matplotlib import rc
from matplotlib import rcParams
import numpy as np

from ase import Atoms
from ase.io import read, write
from ase.build import molecule
def split_line(lines):
    """Split input line"""
    line_array = np.array(lines.strip().split(" "))
    line_vals = line_array[line_array != ""]
    return line_vals
rc('text', usetex=True)
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 20
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
params = {
        'axes.labelsize': 20,
        'font.size': 20,
        'xtick.labelsize': 20,
        'ytick.labelsize': 20,
        'mathtext.fontset': 'stix',
        #'font.family': 'sans-serif',
        #'font.sans-serif': 'Arial',
        'axes.linewidth': 2.0
    }
rcParams.update(params)
fig = plt.figure(figsize=(8.1, 8.1))
gridx, gridy, gridz = np.loadtxt('C6H6.data', unpack=True)
gridx=gridx#*0.52917721
gridy=gridy#*0.52917721
N = int(len(gridz)**.5)
Z = gridz.reshape(N, N)
Z=Z.T
fg=plt.imshow(Z, origin='lower',extent=(np.amin(gridx), np.amax(gridx), np.amin(gridy),  np.amax(gridy)),
        cmap=plt.cm.jet,aspect='equal', interpolation='bicubic')
xVS=(-3.000000e-05)-1.6097968737314778
yVS=0.9295554564124517-1.696604
#plt.plot(xVS, yVS, 'wo', markersize=1, mew=1, color='white')
if os.path.exists('geometry_00.in'):
    print(" geometry file found")
    geometry=open('geometry_00.in','r')
    molc = read('geometry_00.in',format='aims')
    COM = molc.get_center_of_mass(scaled=False)
    print('COM',COM)
    V=Atoms(positions=[COM])
    tip=[0, 0, -4.614046e+00]
    molc.translate(tip-V.positions)
    molc.write('geometry_00.in',format='aims')
    print('COM',COM)
    geometry.close
n_line=0
lines=geometry.readlines()
ii=0
coord=np.zeros((155,3))
for line in lines:
    if line.rfind('atom')!=-1:
       coord[ii,:]= np.float64(split_line(line)[1:4])
       ii=ii+1
if coord is not None:
   x = coord[:, 0]
   y = coord[:, 1]
   z = coord[:, 2]
   a = np.argsort(z)
   plt.plot(y-1.696604,x-0.000030, 'wo', markersize=1, mew=1, color='k')
##################################
plt.title('TERS Image' )
plt.xlabel('Distance, \AA')
#plt.xlim(-5,3.2)
plt.ylabel('Distance, \AA')
cbar = plt.colorbar(fg, ticks=[Z.min(), Z.max()])
cbar.set_label('Intensity', rotation=270,  labelpad=20, y=0.5)
cbar.ax.set_yticklabels(['low', 'high'])  # vertically oriented colorbar
plt.tight_layout()
plt.savefig('VS.png',transparent=True, dpi=400)
plt.show()
