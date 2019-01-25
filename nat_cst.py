from __future__ import print_function

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
    ''' Returns the Planck function :math:`B_{\\nu}(T)` in units of
    :math:`\\rm erg/s/cm^2/Hz/steradian`.

    Args:
        T (float):
            Temperature in K.
        nu:
            numpy array containing the frequency in Hz.
    '''
    
    retVal = 2.*h*nu**3./c**2.
    retVal = retVal / (np.exp(h*nu/kB/T)-1.)
    return retVal

def guillot_global(P,kappa_IR,gamma,grav,T_int,T_equ):
    ''' Returns a temperature array, in units of K,
    of the same dimensions as the pressure P
    (in bar). For this the temperature model of Guillot (2010)
    is used (his Equation 29).

    Args:
        P:
            numpy array of floats, containing the input pressure in bars.
        kappa_IR (float):
            The infrared opacity in units of :math:`\\rm cm^2/s`.
        gamma (float):
            The ratio between the visual and infrated opacity.
        grav (float):
            The planetary surface gravity in units of :math:`\\rm cm/s^2`.
        T_int (float):
            The planetary internal temperature (in units of K).
        T_equ (float):
            The planetary equilibrium temperature (in units of K).
    '''
    tau = P*1e6*kappa_IR/grav
    T_irr = T_equ*np.sqrt(2.)
    T = (0.75 * T_int**4. * (2. / 3. + tau) + \
      0.75 * T_irr**4. / 4. * (2. / 3. + 1. / gamma / 3.**0.5 + \
      (gamma / 3.**0.5 - 1. / 3.**0.5 / gamma)* \
      np.exp(-gamma * tau *3.**0.5)))**0.25
    return T

##################################################################
### Radtrans utility for retrieval temperature model computation
##################################################################

### Box car conv. average, found on stackoverflow somewhere
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

### Global Guillot P-T formula with kappa/grav replaced by delta
def guillot_global_ret(P,delta,gamma,T_int,T_equ):

    delta = np.abs(delta)
    gamma = np.abs(gamma)
    T_int = np.abs(T_int)
    T_equ = np.abs(T_equ)
    tau = P*1e6*delta
    T_irr = T_equ*np.sqrt(2.)
    T = (0.75*T_int**4.*(2./3.+tau) + \
      0.75*T_irr**4./4.*(2./3.+1./gamma/3.**0.5+ \
                         (gamma/3.**0.5-1./3.**0.5/gamma)* \
                             np.exp(-gamma*tau*3.**0.5)))**0.25
    return T

### Modified Guillot P-T formula
def guillot_modif(P,delta,gamma,T_int,T_equ,ptrans,alpha):
    return guillot_global_ret(P,np.abs(delta),np.abs(gamma), \
                                  np.abs(T_int),np.abs(T_equ))* \
                                  (1.-alpha*(1./(1.+ \
                                                np.exp((np.log(P/ptrans))))))

### Function to make temp
def make_press_temp(rad_trans_params):

    press_many = np.logspace(-8,5,260)
    t_no_ave = guillot_modif(press_many, \
        1e1**rad_trans_params['log_delta'],1e1**rad_trans_params['log_gamma'], \
        rad_trans_params['t_int'],rad_trans_params['t_equ'], \
        1e1**rad_trans_params['log_p_trans'],rad_trans_params['alpha'])

    # new
    press_many_new = 1e1**running_mean(np.log10(press_many), 25)
    t_new          = running_mean(t_no_ave  , 25)
    index_new      = (press_many_new <= 1e3) & (press_many_new >= 1e-6)
    temp_new       = t_new[index_new][::2]
    press_new      = press_many_new[index_new][::2]
    
    return press_new, temp_new

