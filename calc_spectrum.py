#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pylab as plt
import radtrans as rt
import nat_cst as nc
import sys as sys
plt.ion()

### Create RT object
line_species = ['H2O']
rayleigh_species = ['H2','He']
ret = rt.radtrans(line_species=line_species,rayleigh_species=rayleigh_species)
press = np.logspace(-6,3.,90)
ret.setup_opa_structure(press)

temp = np.ones_like(press)*960./np.sqrt(2.)
abunds = {}
abunds['H2O'] = np.ones_like(press)*1e-4
abunds['H2'] = np.ones_like(press)*0.75
abunds['He'] = np.ones_like(press)*0.24
grav = 1e1**2.8
MMW = np.ones_like(press)*2.3
r_pl = 0.936*nc.r_jup_mean

ret.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl)
plt.plot(nc.c/ret.freq/1e-4,ret.transm_rad/nc.r_earth)
print(1)

ret.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl,haze_factor=3.)
plt.plot(nc.c/ret.freq/1e-4,ret.transm_rad/nc.r_earth)
print(2)

ret.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl,Pcloud=1e-2)
plt.plot(nc.c/ret.freq/1e-4,ret.transm_rad/nc.r_earth)
print(3)

plt.xlabel('Wavelength (micron)')
plt.ylabel(r'Radius (R$_{\rm Earth}$)')
plt.xticks([0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.],['0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.','2.','3.','4.','5.','6.','7.','8.','9.','10.'])
plt.xscale('log')
input('Press Enter!')
