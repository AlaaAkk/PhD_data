import pandas as pd
import math
import numpy as np
from numpy import array, reshape, zeros, append, arange, ones

### conversions and constants ###
k = 8.617333262145e-05  # ev/k
hb = 6.582119569e-16  # eV.s


def free_energy(omega, Ti, Tf, step):
    F = []

    omega = np.array(omega)
    for T in arange(Ti, Tf, step):
        if T == 0:
            temp3 = np.array([(hb * i / 2) for i in omega])
            F.append(np.sum(temp3))
        else:
            temp3 = np.array(
                [
                    (hb * i / 2 + k * T * np.log(1 - math.exp(-(hb * i) / (k * T))))
                    for i in omega
                ]
            )
            F.append(np.sum(temp3))
    return F


def DeltaF(X, Y, Ti, Tf, step):
    deltaF = []
    F1 = free_energy(X, Ti, Tf, step)
    F2 = free_energy(Y, Ti, Tf, step)
    zip_object = zip(F1, F2)
    for i, j in zip_object:
        deltaF.append(i - j)
    return deltaF
