#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pylab as plt
import radtrans as rt
import nat_cst as nc
import sys as sys
#plt.ion()

'''

### Create RT object
line_species = ['CO_main_iso']
rayleigh_species = ['H2','He']
ret = rt.radtrans(line_species=line_species,rayleigh_species=rayleigh_species,mode='lbl')
ret2 = rt.radtrans(line_species=['CO'],rayleigh_species=rayleigh_species,mode='c-k')
press = np.logspace(-6,3.,90)
ret.setup_opa_structure(press)
ret2.setup_opa_structure(press)

temp = np.ones_like(press)*960./np.sqrt(2.)
abunds = {}
abunds['CO_main_iso'] = np.ones_like(press)*1e-4
abunds['CO'] = np.ones_like(press)*1e-4
abunds['H2'] = np.ones_like(press)*0.75
abunds['He'] = np.ones_like(press)*0.24
grav = 1e1**2.8
MMW = np.ones_like(press)*2.3
r_pl = 0.936*nc.r_jup_mean

ret.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl)
ret2.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl)
plt.plot(nc.c/ret.freq/1e-4,ret.transm_rad/nc.r_earth,label='Clear',color='black')
plt.plot(nc.c/ret2.freq/1e-4,ret2.transm_rad/nc.r_earth,color='black',linestyle=':')
print(1)

ret.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl,haze_factor=3.)
ret2.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl,haze_factor=3.)
plt.plot(nc.c/ret.freq/1e-4,ret.transm_rad/nc.r_earth,label=r'Haze (H$_2$/He Rayleigh scattering $\times$ 3)',color='orange')
plt.plot(nc.c/ret2.freq/1e-4,ret2.transm_rad/nc.r_earth,color='orange',linestyle=':')
print(2)

ret.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl,Pcloud=1e-2)
ret2.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl,Pcloud=1e-2)
plt.plot(nc.c/ret.freq/1e-4,ret.transm_rad/nc.r_earth,label='Cloud deck at 0.01 bar',color='blue')
plt.plot(nc.c/ret2.freq/1e-4,ret2.transm_rad/nc.r_earth,color='blue',linestyle=':')
print(3)

plt.xlabel('Wavelength (micron)')
plt.ylabel(r'Radius (R$_{\rm Earth}$)')
#plt.xticks([0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.],['0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.','2.','3.','4.','5.','6.','7.','8.','9.','10.'])
plt.xlim([2.28,2.35])
#input('Press Enter!')
plt.legend(loc='best',frameon=False)
plt.title(r'CO main isotopologue, solid: $\lambda /\Delta \lambda = 10^6$, dotted: $\lambda /\Delta \lambda = 10^3$')
plt.savefig('clouds.pdf',bbox_inches=0.)
plt.show()
plt.clf()



line_species = ['CO_main_iso']
rayleigh_species = ['H2','He']
ret = rt.radtrans(line_species=line_species,rayleigh_species=rayleigh_species,mode='lbl')
ret2 = rt.radtrans(line_species=['CO_all_iso'],rayleigh_species=rayleigh_species,mode='lbl')
press = np.logspace(-6,3.,90)
ret.setup_opa_structure(press)
ret2.setup_opa_structure(press)

temp = np.ones_like(press)*960./np.sqrt(2.)
abunds = {}
abunds['CO_main_iso'] = np.ones_like(press)*1e-4
abunds['CO_all_iso'] = np.ones_like(press)*1e-4
abunds['H2'] = np.ones_like(press)*0.75
abunds['He'] = np.ones_like(press)*0.24
grav = 1e1**2.8
MMW = np.ones_like(press)*2.3
r_pl = 0.936*nc.r_jup_mean

ret.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl)
ret2.calc_transm(temp,abunds,grav,MMW,1e-2,r_pl)
plt.plot(nc.c/ret2.freq/1e-4,ret2.transm_rad/nc.r_earth,color='red',label='All CO isotopologues')
plt.plot(nc.c/ret.freq/1e-4,ret.transm_rad/nc.r_earth,color='black',label='Main CO isotopologue')

plt.xlabel('Wavelength (micron)')
plt.ylabel(r'Radius (R$_{\rm Earth}$)')
#plt.xticks([0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.],['0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.','2.','3.','4.','5.','6.','7.','8.','9.','10.'])
plt.xlim([2.28,2.35])
#input('Press Enter!')
plt.legend(loc='best',frameon=False)
plt.savefig('iso_comp.pdf',bbox_inches=0.)
plt.show()
plt.clf()

'''

rad_trans_params = {}
rad_trans_params['log_delta'] = -5.
rad_trans_params['log_gamma'] = -0.4
rad_trans_params['t_int'] = 600.
rad_trans_params['t_equ'] = 2300.
rad_trans_params['log_p_trans'] = -2.
rad_trans_params['alpha'] = 0.5

press, temp = rt.make_press_temp(rad_trans_params)
plt.plot(temp,press)
plt.yscale('log')
plt.ylim([1e3,1e-6])
plt.show()
plt.clf()

wlen_bords_micron = [2.2,3.]

ret = rt.radtrans(line_species=['CO_main_iso'],H2H2CIA=True,H2HeCIA=True,mode='lbl',wlen_bords_micron=wlen_bords_micron)
ret2 = rt.radtrans(line_species=['CO_all_iso'],H2H2CIA=True,H2HeCIA=True,mode='lbl',wlen_bords_micron=wlen_bords_micron)
ret3 = rt.radtrans(line_species=['CO'],H2H2CIA=True,H2HeCIA=True,mode='c-k',wlen_bords_micron=wlen_bords_micron)

ret.setup_opa_structure(press)
ret2.setup_opa_structure(press)
ret3.setup_opa_structure(press)

abunds = {}
abunds['CO_main_iso'] = np.ones_like(press)*1e-4
abunds['CO_all_iso'] = np.ones_like(press)*1e-4
abunds['CO'] = np.ones_like(press)*1e-4
abunds['H2'] = np.ones_like(press)*0.75
abunds['He'] = np.ones_like(press)*0.24
grav = 1e1**2.8
MMW = np.ones_like(press)*2.3
r_pl = 0.936*nc.r_jup_mean

ret.calc_flux(temp,abunds,grav,MMW)
ret2.calc_flux(temp,abunds,grav,MMW)
ret3.calc_flux(temp,abunds,grav,MMW)
plt.plot(nc.c/ret2.freq/1e-4,ret2.flux,color='red',label='All CO isotopologues')
plt.plot(nc.c/ret.freq/1e-4,ret.flux,color='black',label='Main CO isotopologue')
plt.plot(nc.c/ret3.freq/1e-4,ret3.flux,color='blue',label='Low res.',linestyle=':')

plt.xlabel('Wavelength (micron)')
plt.ylabel(r'Surface flux (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$)')
plt.xlim([2.28, 2.35])
plt.legend(loc='best',frameon=False)
plt.show()
plt.clf()

wlen = nc.c/ret2.freq/1e-4
wlen = wlen - np.min(wlen)
x_axis = 1./wlen
x_axis = x_axis[::-1]

transfm2 = np.fft.fft(ret2.flux)
plt.plot(x_axis,np.abs(transfm2),color='red',label='All CO isotopologues')
transfm = np.fft.fft(ret.flux)
plt.plot(x_axis,np.abs(transfm),color='black',label='Main CO isotopologue')

plt.yscale('log')
#plt.xscale('log')
plt.xlim([0,100])
plt.show()
plt.clf()

import copy as cp
ind_0 = x_axis < 18.5

transfm2_new = cp.copy(transfm2)
transfm2_new[ind_0] = 0.
flux_filter2 = np.fft.ifft(transfm2_new)
transfm_new = cp.copy(transfm)
transfm_new[ind_0] = 0.
flux_filter = np.fft.ifft(transfm_new)
plt.plot(nc.c/ret2.freq/1e-4,flux_filter2,color='red',label='All CO isotopologues')
plt.plot(nc.c/ret.freq/1e-4,flux_filter,color='black',label='Main CO isotopologue')

plt.xlabel('Wavelength (micron)')
plt.ylabel(r'Surface flux (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$)')
plt.xlim([2.28,2.35])
plt.legend(loc='best',frameon=False)
plt.show()
plt.clf()

ind_0 = x_axis < 60.

transfm2_new = cp.copy(transfm2)
transfm2_new[ind_0] = 0.
flux_filter2 = np.fft.ifft(transfm2_new)
transfm_new = cp.copy(transfm)
transfm_new[ind_0] = 0.
flux_filter = np.fft.ifft(transfm_new)
wlen = nc.c/ret2.freq/1e-4

plt.plot(nc.c/ret2.freq/1e-4,flux_filter2,color='red',label='All CO isotopologues')
plt.plot(nc.c/ret.freq/1e-4,flux_filter,color='black',label='Main CO isotopologue')

plt.xlabel('Wavelength (micron)')
plt.ylabel(r'Surface flux (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$)')
plt.xlim([2.28,2.35])
plt.legend(loc='best',frameon=False)
plt.show()
plt.clf()

wlen = nc.c/ret2.freq/1e-4
ind = wlen > 2.29
y_subs = 1.43e-7+(9.9e-8-1.43e-7)/(2.34-2.29)*(wlen-2.29)

y_subs2 = 1.8e-9+(5.6e-9-1.8e-9)/(2.34-2.29)*(wlen-2.29)

flux_filter2 = flux_filter2 - y_subs - y_subs2
flux_filter = flux_filter - y_subs - y_subs2

plt.plot(nc.c/ret2.freq[ind]/1e-4,flux_filter2[ind],color='red',label='All CO isotopologues')
plt.plot(nc.c/ret.freq[ind]/1e-4,flux_filter[ind],color='black',label='Main CO isotopologue')

plt.xlabel('Wavelength (micron)')
plt.ylabel(r'Surface flux (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$)')
plt.xlim([2.28,2.35])
plt.legend(loc='best',frameon=False)
plt.show()
plt.clf()

wlen = nc.c/ret2.freq/1e-4
ind = wlen > 2.29
from scipy.signal import correlate
corr = correlate(ret2.flux[ind],ret.flux[ind])
corr2 = correlate(ret.flux[ind],ret.flux[ind])
wlen = wlen[ind]

plotrange = 120.

xval = np.linspace(-(np.max(wlen)-np.min(wlen))/np.mean(wlen)*nc.c/1e5,(np.max(wlen)-np.min(wlen))/np.mean(wlen)*nc.c/1e5,int(len(corr)))
plt.plot(xval[np.abs(xval)<plotrange],corr[np.abs(xval)<plotrange]/np.max(corr[np.abs(xval)<plotrange]),label='Correlate all with main isotopolgues',color='black')
plt.plot(xval[np.abs(xval)<plotrange],corr2[np.abs(xval)<plotrange]/np.max(corr[np.abs(xval)<plotrange]),label='Correlate main with main isotopolgues',color='red')

#plt.yscale('log')
plt.xlabel('Shift (km/s)')
plt.xlim([-plotrange,plotrange])
plt.ylabel(r'Normalized cross-correlation value')
plt.legend(loc='best',frameon=False)
plt.show()
plt.clf()

'''
wlen = nc.c/ret3.freq/1e-4
wlen_smooth = rt.box_car_conv(wlen,10)
spec_smooth = rt.box_car_conv(ret3.flux,10)

from scipy.interpolate import interp1d
from scipy.signal import correlate

smooth = interp1d(wlen_smooth, spec_smooth)
ret.flux = ret.flux/smooth(nc.c/ret.freq/1e-4)
ret2.flux = ret2.flux/smooth(nc.c/ret2.freq/1e-4)

plt.plot(nc.c/ret2.freq/1e-4,ret2.flux,color='red',label='All CO isotopologues')
plt.plot(nc.c/ret.freq/1e-4,ret.flux,color='black',label='Main CO isotopologue')

plt.xlabel('Wavelength (micron)')
plt.ylabel(r'Normalized, filtered flux')
plt.xlim([2.2,2.4])
plt.legend(loc='best',frameon=False)
plt.show()
plt.clf()

corr = correlate(ret2.flux,ret.flux)

plt.plot(np.linspace(-(np.max(nc.c/ret2.freq/1e-4)-np.min(nc.c/ret2.freq/1e-4))/np.mean(nc.c/ret2.freq/1e-4)*nc.c/1e5,(np.max(nc.c/ret2.freq/1e-4)-np.min(nc.c/ret2.freq/1e-4))/np.mean(nc.c/ret2.freq/1e-4)*nc.c/1e5,int(len(corr))),correlate(ret2.flux,ret.flux))
#plt.yscale('log')
plt.xlabel('Shift (km/s)')
plt.xlim([-200,200])
plt.ylabel(r'Cross-correlation value')

plt.show()
plt.clf()
'''
