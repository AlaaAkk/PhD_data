import numpy as np
from numpy import array, reshape, zeros, append, arange, ones
from read_inputs import *
from chemical_pot import *


def Form(mu_X,mu_M,E, E0, sign, n, mat=False):
    Ef = []
    if mat == False:
        for i in mu_X:
            Ef.append(E - E0 + (sign * n * i))
    else:
        for i in mu_M:
            Ef.append(E - E0 + (sign * n * i))
    return Ef


def Form_pbe(mu_X_pbe,mu_M_pbe,E, E0, sign, n, mat=False):
    Ef = []
    if mat == False:
        for i in mu_X_pbe:
            Ef.append(E - E0 + (sign * n * i))
    else:
        for i in mu_M_pbe:
            Ef.append(E - E0 + (sign * n * i))
    return Ef



