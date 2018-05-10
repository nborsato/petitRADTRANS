""" nat_cst.py: provides natural constants and functions such as the black body intensity for the petitRADTRANS module. """

__author__ = "Paul Molliere"
__copyright__ = "Copyright 2016-2018, Paul Molliere"
__maintainer__ = "Paul Molliere"
__email__ = "molliere@strw.leidenunivl.nl"
__status__ = "Development"


import numpy as np

# Natural constants
# Everything is in cgs!
c = 2.99792458e10
G = 6.67384e-8
kB = 1.3806488e-16
sigma = 5.670373e-5
bar = 1e6
atm = 1.01325e6
losch = 2.68676e19
h = h=6.62606957e-27
r_sun = 6.955e10
r_jup=7.1492e9
r_jup_mean=6.9911e9
m_jup = 1.89813e30
m_sun = 1.9891e33
AU = 1.49597871e13
pc = 3.08567758e18
amu = 1.66053892e-24
nA = 6.0221413e23
R = 8.3144598
m_earth = 5.9722e27
r_earth = 637813660.

# Molecular weights in amu
molecular_weight = {}
molecular_weight['H2O'] = 18.
molecular_weight['CH4'] = 16.
molecular_weight['CO2'] = 44.
molecular_weight['CO'] = 28.
molecular_weight['H2'] = 2.
molecular_weight['He'] = 4.

def b(T,nu):
    retVal = 2.*h*nu**3./c**2.
    retVal = retVal / (np.exp(h*nu/kB/T)-1.)
    return retVal
