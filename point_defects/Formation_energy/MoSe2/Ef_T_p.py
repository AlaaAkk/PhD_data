import pandas as pd
import math
import numpy as np
from numpy import array, reshape, zeros, append, arange, ones
from Ef_0K import *
from read_inputs import *
from chemical_pot import *
from free_energy import *

### conversions and constants ###
Pi = np.pi
cmhz = 29979245800.0 * 2 * Pi  # cm^-1 to Hz
amu2g = 1.660539e-24
ang2bohr = 1.8897261
ang2cm = 1e-08


p0 = 1013250  # atm to g/(cm s^2)

kk = 1.380649e-16  # erg/k (cm^2.g/ks^2)
k = 8.617333262145e-05  # ev/k
h = 6.62607015e-27  # erg.s
hb = 6.582119569e-16  # eV.s
hbar = 1.054571817e-27  # erg.s


def formation_energy(E_MX2,
    ES, pi, pf, pstep, Ti, Tf, tstep, sigma, E, E0, sign, n, w, w0, wX8, mat=False
):

    pX = pressure(pi, pf, pstep)
    T = np.arange(Ti, Tf, tstep)
    Ef = np.zeros((len(pX), len(T)))
    concent = np.zeros((len(pX), len(T)))
    mu_X, mu_M = mu(Ti, Tf, tstep, ES, sigma, pi, pf, pstep,wX8,E_MX2)
    pindx = 0
    for pindx, p in enumerate(list(pX * p0)):
        # print(pindx)
        if mat == False:
            Ef[pindx, :] = np.array(
                [
                    E - E0 + (sign * n * a) + b
                    for a, b in zip(mu_X[pindx, :], DeltaF(w, w0, Ti, Tf, tstep))
                ]
            )
        # concent[pindx,:] = np.array(conc(Ef[pindx,:],1,Ti,Tf,tstep))
        else:
            Ef[pindx, :] = np.array(
                [
                    E - E0 + (sign * n * a) + b
                    for a, b in zip(mu_M[pindx, :], DeltaF(w, w0, Ti, Tf, tstep))
                ]
            )
        #  concent[pindx,:] = np.array(conc(Ef[pindx,:],1,Ti,Tf,tstep))
    #    yield concent
    return Ef

