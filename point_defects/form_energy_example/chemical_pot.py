''' Chemical Potential of X8 '''
import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.io import read , write
from numpy import array, reshape, zeros, append, arange, ones
from read_inputs import *
import sys
import math

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


def chemical_pot_1(x, y, z):
    '''Boundaries of Chemical Potential'''
    # x is mu_MX2
    # y is \mu_M
    # z is \mu_X8

    mu_X_i = 0.5 * (x - y)
    mu_X_f = z
    mu_M_i = x - 2 * mu_X_i
    mu_M_f = x - 2 * mu_X_f

    mu_X = arange(mu_X_i, mu_X_f, 0.05)
    mu_M = arange(mu_M_f, mu_M_i, 0.05)
    var_X = mu_X - mu_X_f
    var_M = (mu_M - mu_M_f) / 2
    yield np.array(mu_X)
    yield np.array(mu_M)
    yield var_X
    yield var_M
  #  yield mu_X_i
  #  yield mu_X_f

def moment():
    ''' Moment of Inertia of X8'''
    molc = read("X8.in", format="aims")
    mass = molc.get_masses()
    mass = sum(mass) * amu2g
    Is = molc.get_moments_of_inertia()
    Is = Is * amu2g * ang2cm * ang2cm
    I = np.sqrt(Is[0]) * np.sqrt(Is[1]) * np.sqrt(Is[2])
    yield Is
    yield I
    yield mass

def pressure(pi, pf, pstep):
    pX = np.arange(pi, pf, pstep)
    return pX

def translational(Ti, Tf, tstep):
    '''Translational Part of Chemical Potential'''
    A = []
    trans = []
    Is, I, m = moment()
    Tindx = 0
    for Tindx, T in enumerate(list(arange(Ti, Tf, tstep))):
        if T == 0:
            A = 0
            trans.append(-k * T * A)
        else:
            A = np.log(
                (((2 * Pi * m) ** (3 / 2)) * ((kk * T) ** (5 / 2)))
                / (p0 * (h ** 3))
            )
            trans.append(-k * T * A)
    return np.array(trans)

def rotational(Ti, Tf, tstep, sigma):
    '''Rotational Part of Chemical Potential'''
    B = []
    rot = []
    Is, I, m = moment()
    Tindx = 0
    for Tindx, T in enumerate(list(arange(Ti, Tf, tstep))):
        if T == 0:
            B = 0
            rot.append(-k * T * B)
        else:
            B = np.log(np.sqrt(Pi) / sigma) + np.log(
                (((8 * Pi * kk * T) / (h ** 2)) ** (3 / 2)) * I
            )
            rot.append(-k * T * B)
    # print('rot', np.shape(rot))
    return np.array(rot)

def vibrational(Ti, Tf, tstep, elec,wX8):
    '''Vibrational Part of Chemical Potential'''
    C = []
    D = []
    vib = []
    Is, I, m = moment()
    Tindx = 0
    for Tindx, T in enumerate(list(arange(Ti, Tf, tstep))):
        temp2 = np.array([(hb * i) / (2) for i in wX8])
        D = np.sum(temp2)
        if T == 0:
            vib.append(D + elec)
        else:
            temp = np.array(
                [(np.log(1 - math.exp(-(hbar * i) / (kk * T)))) for i in wX8]
            )
            C = np.sum(temp)
            vib.append(k * T * C + D + elec)
    return vib

def mu_0(Ti, Tf, tstep, elec, sigma,wX8):
    '''Chemical Potential Tempreture Dependent Terms'''
    mu_0 = []
    rot = rotational(Ti, Tf, tstep, sigma)
    vib = vibrational(Ti, Tf, tstep, elec,wX8)
    trans = translational(Ti, Tf, tstep)
    mu_0 = rot + vib + trans
    return mu_0

def press(pi, pf, pstep, Ti, Tf, tstep):
    '''Pressure Term of the Chemical Potential'''
    pindx = 0
    pX = pressure(pi, pf, pstep)
    T = np.arange(Ti, Tf, tstep)
    F = np.zeros((len(pX), len(T)))
    for pindx, p in enumerate(list(pX * p0)):
        if p == 0:
            F[pindx, :] = np.zeros(len(T))
        else:
            if Ti == 0:
                F[pindx, 0] = 0
                F[pindx, 1:] = np.array(
                    [
                        k * T * np.log(p / p0)
                        for T in arange(Ti, Tf, tstep)
                        if T != 0
                    ]
                )
            else:
                F[pindx, :] = np.array(
                    [k * T * np.log(p / p0) for T in arange(Ti, Tf, tstep)]
                )
    return F

def mu(Ti, Tf, tstep, elec, sigma, pi, pf, pstep,wX8,mu_MX2):
    '''Chemical Potential as a function of Temperature and Pressure '''
    pindx = 0
    pX = pressure(pi, pf, pstep)
    T = np.arange(Ti, Tf, tstep)
    mu_10 = np.zeros(len(T))
    mu_10 = mu_0(Ti, Tf, tstep, elec, sigma,wX8)
    p_dep = np.zeros((len(pX), len(T)))
    p_dep = press(pi, pf, pstep, Ti, Tf, tstep)
    mu_X = np.zeros((len(pX), len(T)))
    mu_X8 = np.zeros((len(pX), len(T)))
    mu_M = np.zeros((len(pX), len(T)))
    for pindx, p in enumerate(list(pX * p0)):
        mu_X8[pindx, :] = np.array(mu_10 + p_dep[pindx, :])
        mu_X[pindx, :] = mu_X8[pindx, :] / 8
        mu_M[pindx, :] = mu_MX2 - 2 * mu_X[pindx, :]
    yield mu_X
    yield mu_M

