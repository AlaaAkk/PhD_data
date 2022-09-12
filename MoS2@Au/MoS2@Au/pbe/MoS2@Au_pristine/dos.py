
''' Plotting Intensities vs frequencies using Vlines '''
import  pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os #like in the terminal

f = open('aims.out', 'r')
FermiValues = list()
VBMValues = list()
gapValues = list()
CBMValues = list()
lines = f.readlines()
for line in lines:
    if "(Fermi" in line.split(" "):
        try:
            FermiValues.append(float(line.split(" ")[-2]))
        except:
            continue
    if "(VBM)" in line.split(" "):
        try:
            VBMValues.append(float(line.split(" ")[-6]))
        except:
            continue
    if "HOMO-LUMO" in line.split(" "):
        try:
            gapValues.append(float(line.split(" ")[-12]))
        except:
            continue
    if "(CBM)" in line.split(" "):
        try:
            CBMValues.append(float(line.split(" ")[-6]))
        except:
            continue
VBM = VBMValues[-1]
fermi = FermiValues[-1]
print('fermi',fermi)
print('VBM',VBM)
Eg = gapValues[-1]
CBM = CBMValues[-1]
print('CBM',CBM)
print('Eg',Eg)
#print(Eg)
# Or custom limits
if (VBM !=0) and (Eg !=0):
  # Adding lines and labels for VBM and CBM
  plt.axvline(x=VBM-VBM, color='m', linestyle='--')
 # plt.axvline(x=Eg, color='m', linestyle='--')
  plt.axvline(x=fermi-VBM, color='k', linestyle='--')
  plt.axvline(x=CBM-VBM, color='c', linestyle='--')

d0=pd.read_csv( "MoS2.dat" ,sep='\s+', header=0) # total MoS2 8x8 Pristine
d2=pd.read_csv( "Mo_dn.dat" ,sep='\s+', header=0) # total MoS2 8x8 Pristine
d1=pd.read_csv( "Mo_up.dat" ,sep='\s+', header=0) # total MoS2 8x8 Pristine
d3=pd.read_csv( "S_dn.dat" ,sep='\s+', header=0) # total MoS2 8x8 Pristine
d4=pd.read_csv( "S_up.dat" ,sep='\s+', header=0) # total MoS2 8x8 Pristine
#d3['up-up1'] = d0.up - d1.up1
up=d1['up']+d4['up']
dn=d2['down']+d3['down']
plt.plot(d0['Energy']-VBM, d0['up'], 'm-', label='MoS2 spin up')
plt.plot(d0['Energy']-VBM, -d0['down'], 'c-', label='MoS2 spin down')

f = open('aims_2.out', 'r')
FermiValues = list()
lines = f.readlines()
for line in lines:
    if "(Fermi" in line.split(" "):
        try:
            FermiValues.append(float(line.split(" ")[-2]))
        except:
            continue
#VBM2 = VBMValues[-1]
fermi = FermiValues[-1]
print('fermi',fermi)
#print(Eg)
# Or custom limits
#if (VBM !=0) and (Eg !=0):
#  # Adding lines and labels for VBM and CBM
#  plt.axvline(x=VBM-VBM, color='b', linestyle='--')
# # plt.axvline(x=Eg, color='m', linestyle='--')
#  plt.axvline(x=fermi-VBM, color='k', linestyle='--')
#  plt.axvline(x=CBM-VBM, color='r', linestyle='--')

plt.plot(d1['Energy']-fermi+1.3, up, 'r-', label='MoS2@Au')
plt.plot(d1['Energy']-fermi+1.3, -dn, 'b-', label='MoS2@Au')
plt.legend(fontsize=10)
plt.xlabel('E-$\mu$', fontsize=12)
plt.ylabel('Energy [eV]', fontsize=12)
plt.tight_layout() #better
plt.xlim([-2, 2])
plt.ylim([-400, 400])
plt.show()
plt.savefig(os.path.join( 'total_dos.png'), dpi=400)
