import ase
from ase.io import read, write
import os
import pandas as pd
import matplotlib.pyplot as plt
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
#  plt.axvline(x=fermi-VBM, color='k', linestyle='--')
  plt.axvline(x=CBM-VBM, color='c', linestyle='--')
workdir = os.getcwd()
colnames=['Index Energy DOSup DOSdown']
total_dos = pd.read_csv(os.path.join(workdir, "KS_DOS_total_raw_tetrahedron.dat"),
                                                skiprows=4,
                                                skipinitialspace = True,
                                                header=None,
                                                sep="\s+",
                                                names = ['Energy', 'DOS(spin up)', 'DOS(spin down)'])
elem_df = pd.DataFrame()
elem_df["Energy"] = total_dos["Energy"]
elem_df["Dosup"] = total_dos["DOS(spin up)"]
elem_df["Dosdown"] = total_dos["DOS(spin down)"]
plt.plot(elem_df["Energy"]-VBM,elem_df["Dosup"] , label="MoS2 up", linewidth=1)
plt.plot(elem_df["Energy"]-VBM,-elem_df["Dosdown"] , label="MoS2 down", linewidth=1)
plt.legend(loc="upper left")
plt.xlim([-4, 6])
plt.ylim([-400, 400])
plt.legend(loc="upper left")
plt.xlabel('E-$\mu$', fontsize=14)
plt.ylabel('Energy [eV]', fontsize=14)
plt.tight_layout() #better
plt.show()
plt.savefig(os.path.join( 'pdos.png'), dpi=400)
