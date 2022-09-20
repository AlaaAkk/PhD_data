# Libraries
import numpy as np
from numpy import array, reshape, zeros, append, arange, ones

### conversions and constants ###
Pi = np.pi
cmhz = 29979245800.0 * 2 * Pi  # cm^-1 to Hz

# Spliting function
def split_line(lines):
    """Split input line"""
    line_array = np.array(lines.strip().split(" "))
    line_vals = line_array[line_array != ""]
    return line_vals


# reading total energies in aims
def read_energy(filename):
    energy = 0
    with open(filename) as out:

        for line in out:
            if (
                line.rfind(
                    "| Total energy of the DFT / Hartree-Fock s.c.f. calculation      : "
                )
                != -1
            ):
                energy = float(split_line(line)[-2])  # Periodic/cluster
                return energy


# reading total energies in aims
def read_omega(filename):
    """Reads frequencies in cm-1 and returns them in hz"""
    ftemp = open(filename)
    fw = np.array([float(line.split()[0]) for line in ftemp]) * cmhz
    ftemp.close()
    return fw
