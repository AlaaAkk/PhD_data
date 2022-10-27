import  pandas as pd
import math
import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.io import read , write
from numpy import array, reshape, zeros, append, arange, ones
import matplotlib
import matplotlib as mpl
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import sys

### conversions and constants ###
Pi=np.pi
cmhz=29979245800.0*2*Pi #cm^-1 to Hz
amu2g=1.660539e-24
ang2bohr=1.8897261
ang2cm=1e-08
q=1.602e-19

p0=1013250  # atm to g/(cm s^2)

kk=1.380649e-16 # erg/k (cm^2.g/ks^2)
k=8.617333262145e-05 # ev/k
h=6.62607015e-27  # erg.s
hb=6.582119569e-16 # eV.s
hbar=1.054571817e-27 # erg.s


# In[3]:

#__all__ = ['main']
# Spliting function
def split_line(lines):
    """Split input line"""
    line_array = np.array(lines.strip().split(" "))
    line_vals = line_array[line_array != ""]
    return line_vals


# In[4]:
# reading total energies in aims
def read_energy(filename):
     energy=0
     with open(filename) as out:

        for line in out:
           if line.rfind('| Total energy of the DFT / Hartree-Fock s.c.f. calculation      : ')!=-1:
              energy = float(split_line(line)[-2]) # Periodic/cluster
              return energy
def Form_S(E,E0,mu_S):
   # chemical_potential(E_MoS2,mu_Mobcc,mu_S8)
    Ef=[]
    fermi=np.arange(-7,-3)
    for i in fermi:
       Ef.append(1*(E-E0-i+mu_S-7.02338194))
    return Ef
def Form_S0(E,E0,mu_S):
   # chemical_potential(E_MoS2,mu_Mobcc,mu_S8)
    Ef=[]
    fermi=np.zeros(4)
    for i in fermi:
       Ef.append(E-E0-i+mu_S)
    return Ef

EM=read_energy('Mo2.out')
E_MX2=read_energy('primitive.out')
mu_Mo=EM/2
E0=read_energy('pristine_charged.out')
E1=read_energy('VS_charged.out')
E00=read_energy('pristine.out')
E10=read_energy('VS.out')
ES8=read_energy('S8.out')
mu_Si=(E_MX2-mu_Mo)/2
mu_Sf=ES8/8
EE=Form_S(E1,E0,mu_Sf)
Ei=Form_S(E1,E0,mu_Si)
EE0=Form_S0(E10,E00,mu_Sf)
Ei0=Form_S0(E10,E00,mu_Si)
fermi=np.arange(-7,-3)
#plt.plot(fermi,EE)
fig = plt.figure()
ax1 = fig.add_subplot(111)
for axis in ["top", "bottom", "left", "right"]:
    ax1.spines[axis].set_linewidth(3.0)
ax1.plot(fermi, EE, "b", label="rich env. charged")
ax1.plot(fermi, Ei, "r", label="poor env. charged")
ax1.plot(fermi, EE0, "b", label="rich env. neutral")
ax1.plot(fermi, Ei0, "r", label="poor env. neutral")
ax1.set_xlabel('Fermi Energy eV')
ax1.set_ylabel("Formation energies eV")
ax1.legend(loc=4)
plt.tight_layout()
plt.savefig("GvsEf.pdf")
plt.show()
