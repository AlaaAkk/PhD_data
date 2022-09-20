import pandas as pd
import math
import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.io import read, write
from numpy import array, reshape, zeros, append, arange, ones
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors
from read_inputs import *
from chemical_pot import *
from Ef_0K import *
from free_energy import *
from Ef_T_p import *


mpl.rcParams["font.family"] = "Sans"
mpl.rcParams["lines.linewidth"] = 3
mpl.rcParams["lines.markersize"] = 18
mpl.rcParams["font.size"] = 22  # <-- change fonsize globally
mpl.rcParams["legend.fontsize"] = 18
mpl.rcParams["axes.titlesize"] = 22
mpl.rcParams["axes.labelsize"] = 18
mpl.rcParams["xtick.major.size"] = 4
mpl.rcParams["ytick.major.size"] = 4
mpl.rcParams["xtick.major.width"] = 1.5
mpl.rcParams["ytick.major.width"] = 1.5
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "out"
mpl.rcParams["figure.titlesize"] = 14
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["axes.labelweight"] = "bold"
# mpl.rc('text', usetex=True)

mpl.rcParams["text.latex.preamble"] = [
    r"\renewcommand{\familydefault}{\sfdefault}",
    #    r'\usepackage[scaled=1]{helvet}',
    #  r'\usepackage[helvet]{sfmath}',
    #     r'\everymath={\sf}'
]

def main(filename, Ti, Tf, tstep, pi, pf, pstep):

    """main routine"""
    # Reading Energies
    E0=read_energy('pristine')
    E1=read_energy('addX')
    E2=read_energy('VX')
    E3=read_energy('VX2')
    E4=read_energy('VX22')
    E5=read_energy('VM')
    E_MX2=read_energy('primitive')
    EX8=read_energy('EX8')
    EM=read_energy('EM')
    # Chemical Potential
    mu_MX2=E0/75
    mu_X=EX8/8
    mu_Mbcc=EM/2
    # pbe:
    E0_pbe=read_energy('pristine_pbe')
    E1_pbe=read_energy('addX_pbe')
    E2_pbe=read_energy('VX_pbe')
    E3_pbe=read_energy('VX2_pbe')
    E4_pbe=read_energy('VX22_pbe')
    E5_pbe=read_energy('VM_pbe')
    E_MX2_pbe=read_energy('primitive_pbe')
    EX8_pbe=read_energy('EX8_pbe')
    EM_pbe=read_energy('EM_pbe')
    # Chemical Potential
    mu_MX2_pbe=E0_pbe/75
    mu_X_pbe=EX8_pbe/8
    mu_Mbcc_pbe=EM_pbe/2
    # Reading the frequencies
    wM=read_omega('wM')
    wX8=read_omega('wX8')
    w0=read_omega('w0')
    w1=read_omega('waddX')
    w2=read_omega('wVX')
    w3=read_omega('wVX2')
    w4=read_omega('wVX22')
    w5=read_omega('wVM')
    # Calculation Formation Energies
    mu_X, mu_M, x_X, x_M= chemical_pot_1(E_MX2, mu_Mbcc, mu_X)
    mu_X_pbe, mu_M_pbe, x_X_pbe, x_M_pbe= chemical_pot_1(E_MX2_pbe,mu_Mbcc_pbe,mu_X_pbe)


    Ef_addX = Form(mu_X,mu_M,E1, E0, -1, 1, False)
    Ef_VX = Form(mu_X,mu_M,E2, E0, +1, 1, False)
    Ef_VX2 = Form(mu_X,mu_M,E3, E0, +1, 2, False)
    Ef_VX22 = Form(mu_X,mu_M,E4, E0, +1, 2, False)
    Ef_M = Form(mu_X,mu_M,E5, E0, +1, 1, True)

    Ef_addX_pbe=Form_pbe(mu_X_pbe,mu_M_pbe,E1_pbe,E0_pbe,-1,1,False)
    Ef_VX_pbe=Form_pbe(mu_X_pbe,mu_M_pbe,E2_pbe,E0_pbe,+1,1,False)
    Ef_VX2_pbe=Form_pbe(mu_X_pbe,mu_M_pbe,E3_pbe,E0_pbe,+1,2,False)
    Ef_VX22_pbe=Form_pbe(mu_X_pbe,mu_M_pbe,E4_pbe,E0_pbe,+1,2,False)
    Ef_M_pbe=Form_pbe(mu_X_pbe,mu_M_pbe,E5_pbe,E0_pbe,+1,1,True)

    # Calculating Chemical potential as function of temperature and pressure
    T = [float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])]
    pS = pressure(float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]))
    p = [float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])]
    Ef_addX_Tp = formation_energy(E_MX2,
        EX8,p[0],p[1],p[2],T[0],T[1],T[2],8,E1,E0,-1,1,w1,w0,wX8,mat=False,
    )
    Ef_VX_Tp = formation_energy(E_MX2,
        EX8, p[0], p[1], p[2], T[0], T[1], T[2], 8, E2, E0, 1, 1, w2, w0, wX8, mat=False
    )
    Ef_VX2_Tp = formation_energy(E_MX2,
        EX8, p[0], p[1], p[2], T[0], T[1], T[2], 8, E3, E0, 1, 2, w3, w0, wX8, mat=False
    )
    Ef_VX22_Tp = formation_energy(E_MX2,
        EX8, p[0], p[1], p[2], T[0], T[1], T[2], 8, E4, E0, 1, 2, w4, w0, wX8, mat=False
    )
    Ef_VM_Tp = formation_energy(E_MX2,
        EX8, p[0], p[1], p[2], T[0], T[1], T[2], 8, E5, E0, 1, 1, w5, w0, wX8, mat=True
    )
    def plotting():
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        for axis in ["top", "bottom", "left", "right"]:
            ax1.spines[axis].set_linewidth(3.0)
        ax1.plot(x_X, Ef_addX, "b", label="addX")
        ax1.plot(x_X, Ef_VX22, "g", label='VX22')
        ax1.plot(x_X, Ef_VX, "r", label="VX")
        ax1.plot(x_X, Ef_VX2, "k", label="VX2")
        ax1.plot(x_X_pbe, Ef_addX_pbe, "b", linestyle="dashed")
        ax1.plot(x_X_pbe, Ef_VX_pbe, "r", linestyle="dashed")
        ax1.plot(x_X_pbe, Ef_VX22_pbe, "g", linestyle="dashed")
        ax1.plot(x_X_pbe, Ef_VX2_pbe, "k", linestyle="dashed")
        ax1.set_xlabel(r"$\Delta \mu_{X}$")
        ax1.set_ylabel("Formation energies eV")
        boxProps = dict(
            boxstyle="square",
            facecolor="white",
            alpha=0.0,
            linewidth=0.5,
            fill=True,
            pad=0.3,
        )
        ax2 = ax1.twiny()
        ax2.plot(x_M, Ef_M, "magenta", label="VM")
        ax2.plot(x_M_pbe, Ef_M_pbe, "magenta", linestyle="dashed")
        ax2.set_xlabel(r"$\frac{1}{2}\Delta \mu_{M}$", labelpad=20)
        ax1.legend(loc=4)
        ax2.legend(loc=1)
        plt.tight_layout()
        plt.savefig("Ef_0K.png", bbox_inches="tight", dpi=400)
        plt.savefig("Ef_0K.pdf")
        plt.show()

    def plotting_T(Ti, Tf, tstep):
        #       #pS=pressure(Ti,Tf,tstep)
        #    #===============
        #   #  First subplot
        #   #===============
        #   # set up the axes for the first plot
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        for axis in ["top", "bottom", "left", "right"]:
            ax.spines[axis].set_linewidth(3.0)
        T = arange(Ti, Tf, tstep)
        pindx = 0
        ax.plot(T, Ef_addX_Tp[pindx, :].T, "b", label="addX")
        ax.plot(T, Ef_VX_Tp[pindx, :].T, "r", label="VX")
        ax.plot(T, Ef_VX22_Tp[pindx, :].T, "g", label="VX22")
        ax.plot(T, Ef_VX2_Tp[pindx, :].T, "k", label="VX2")
        ax.plot(T, Ef_VM_Tp[pindx, :].T, "magenta", label="VM")
        ax.legend(loc=4)
        ax.set_xlabel(" Temperature [K]")
        ax.set_ylabel(" Formation Energy [eV]")
        boxProps = dict(
            boxstyle="square",
            facecolor="white",
            alpha=0.0,
            linewidth=0.5,
            fill=True,
            pad=0.3,
        )
       # ax.text(-1.2, 8.3, str(filename), bbox=boxProps)
        plt.tight_layout()
        plt.savefig("Ef_T_" + str(sys.argv[1]) + ".pdf")
        plt.show()

    def plotting_P(pi, pf, pstep):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for axis in ["top", "bottom", "left", "right"]:
            ax.spines[axis].set_linewidth(3.0)
        p = arange(pi, pf, pstep)
        Tindx = 0
        ax.plot(p, Ef_addX_Tp[:, Tindx].T, "b", label="addX")
        ax.plot(p, Ef_VX_Tp[:, Tindx].T, "r", label="VX")
        ax.plot(p, Ef_VX22_Tp[:, Tindx].T, "g", label="VX22")
        ax.plot(p, Ef_VX2_Tp[:, Tindx].T, "k", label="VX2")
        ax.plot(p, Ef_VM_Tp[:, Tindx].T, "magenta", label="VM")
        ax.legend(loc=4)
        ax.set_xlabel(" Pressure [atm]")
        ax.set_ylabel(" Formation Energy [eV]")
        plt.savefig("Temp_" + str(sys.argv[1]) + ".pdf")
        plt.show()


    plotting()
    plotting_T(float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]))
    plotting_P(float(sys.argv[5]),float(sys.argv[6]),float(sys.argv[7]))
if __name__ == "__main__":
    main(str(sys.argv[1]),
        float(sys.argv[2]),
        float(sys.argv[3]),
        int(sys.argv[4]),
        float(sys.argv[5]),
        float(sys.argv[6]),
        float(sys.argv[7]),)
