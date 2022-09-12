import ase
from ase.io import read, write
import os
import pandas as pd
import matplotlib.pyplot as plt
f = open('aims_MoS2.out', 'r')
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
  plt.axvline(x=VBM-VBM, color='b', linestyle='--')
 # plt.axvline(x=Eg, color='m', linestyle='--')
#  plt.axvline(x=fermi-VBM, color='k', linestyle='--')
  plt.axvline(x=CBM-VBM, color='orange', linestyle='--')
workdir = os.getcwd()
colnames=['Index Energy DOSup DOSdown']
total_dos = pd.read_csv(os.path.join(workdir, "KS_DOS_total_raw_tetrahedron.dat"),
                                                skiprows=4,
                                                skipinitialspace = True,
                                                header=None,
                                                sep="\s+",
                                                names = ['Energy', 'DOSup', 'DOSdown'])
elem_df = pd.DataFrame()
elem_df["Energy"] = total_dos["Energy"]
elems = ['Mo','S', 'Au']
for elem in elems:
        proj_up = pd.read_csv(os.path.join(workdir, elem+"_l_proj_dos_tetrahedron_spin_up_raw.dat"),
                                                      skiprows=4,
                                                      skipinitialspace = True,
                                                      header=None,
                                                      sep="\s+")
        proj_down = pd.read_csv(os.path.join(workdir, elem+"_l_proj_dos_tetrahedron_spin_dn_raw.dat"),
                                                        skiprows=4,
                                                        skipinitialspace = True,
                                                        header=None,
                                                        sep="\s+")
        elem_df["{}_up".format(elem)] = proj_up[1]
        elem_df["{}_dn".format(elem)] = proj_down[1]
plt.plot(elem_df["Energy"]-VBM+0.0,
              elem_df["Mo_up"]+elem_df["S_up"], label="MoS2@Au pristine up", linewidth=1)
plt.plot(elem_df["Energy"]-VBM+0.0,
              -elem_df["Mo_dn"]-elem_df["S_dn"], label="MoS2@Au pristine down", linewidth=1)
plt.xlim([-3, 4])
plt.ylim([-300, 300])
plt.legend(loc="upper left")
plt.xlabel('E-$\mu$', fontsize=14)
plt.ylabel('Energy [eV]', fontsize=14)
plt.tight_layout() #better
plt.savefig(os.path.join( 'pdos.png'), dpi=400)
plt.show()
