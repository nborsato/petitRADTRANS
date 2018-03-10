from . import fort_input as fi
from . import fort_spec as fs
from . import nat_cst as nc

import numpy as np
import copy as cp
import os
import sys

class radtrans:
    """ Class carrying out spectral calcs for a given set of opacities """
    
    def __init__(self,line_species=[],rayleigh_species=[],cloud_species=[], \
                     H2H2CIA=False,H2HeCIA=False,wlen_bords_micron=[0.05,300.], \
                     mode='c-k'):

        # Line-by-line or corr-k
        self.mode = mode

        # Line opacity species to be considered
        self.line_species = line_species

        # Rayleigh scattering species to be considered
        self.rayleigh_species = rayleigh_species

        # Cloud species to be considered
        self.cloud_species = cloud_species

        # Include CIA?
        self.H2H2CIA = H2H2CIA
        self.H2HeCIA = H2HeCIA

        # Get path to all input data (opacities, grids, etc.)
        f = open(os.path.dirname(__file__)+'/path.txt')
        lines = f.readlines()
        self.path = os.path.dirname(__file__)+'/'+lines[1].rstrip()
        f.close()

        # Read in frequency grid
        if self.mode == 'c-k':
            # For correlated-k

            # Get dimensions of molecular opacity arrays for a given P-T point,
            # they define the resolution.
            self.freq_len, self.g_len = fi.get_freq_len(self.path)
            freq_len_full = cp.copy(self.freq_len)
            # Read in the frequency range of the opcity data
            self.freq = fi.get_freq(self.path,self.freq_len)
            arr_min, arr_max = -1, -1
            
        elif self.mode == 'lbl':
            # For high-res line-by-line radiative transfer
            path_length = self.path+'/opacities/lines/line_by_line/'+ \
                                line_species[0]+'/wlen.dat'
            # Get dimensions of opacity arrays for a given P-T point
            # arr_min, arr_max denote where in the large opacity files
            # the required wavelength range sits.
            self.freq_len, arr_min, arr_max = \
              fi.get_arr_len_array_bords(wlen_bords_micron[0]*1e-4, \
                                         wlen_bords_micron[1]*1e-4, \
                                         path_length)

            self.g_len = 1
            freq_len_full = self.freq_len
            # Read in the frequency range of the opcity data
            wlen = fi.read_wlen(arr_min, arr_max, self.freq_len, path_length)
            self.freq = nc.c/wlen

        if self.mode == 'c-k':
            # Cut opacity data to required range
            index = (nc.c/self.freq > wlen_bords_micron[0]*1e-4) & \
              (nc.c/self.freq < wlen_bords_micron[1]*1e-4)

            self.freq = np.array(self.freq[index],dtype='d',order='Fortran')
            self.freq_len = len(self.freq)

        ###########################
        # Some necessary definitions, also prepare arrays for fluxes, transmission radius...
        ###########################        

        self.lambda_angstroem = np.array(nc.c/self.freq/1e-8,dtype='d',order='Fortran')
        self.flux = np.array(np.zeros(self.freq_len),dtype='d',order='Fortran')
        self.transm_rad = np.array(np.zeros(self.freq_len),dtype='d',order='Fortran')

        # Define frequency bins around grid for later interpolatioin purposes when including
        # clouds...
        self.border_freqs = np.array(nc.c/self.calc_borders(nc.c/self.freq),dtype='d',order='Fortran')

        self.Pcloud = None
        self.haze_factor = None
        self.gray_opacity = None

        ###########################
        # Read in opacities
        ###########################
        
        # Molecules:
        # First get the P-Ts where the grid is defined.
        buffer = np.genfromtxt(self.path+'/opa_input_files/opa_PT_grid.dat')
        self.line_TP_grid = np.zeros_like(buffer)
        self.line_TP_grid[:,0] = buffer[:,1]
        self.line_TP_grid[:,1] = buffer[:,0]
        # Convert from bars to cgs
        self.line_TP_grid[:,1] = 1e6*self.line_TP_grid[:,1]
        self.line_TP_grid = np.array(self.line_TP_grid.reshape(len(self.line_TP_grid[:,1]),2),dtype='d',order='Fortran')

        # Read actual opacities....
        # line_grid_kappas has the shape g_len,freq_len,len(line_species),len(line_TP_grid[:,0])
        if len(self.line_species) > 0:

            tot_str = ''
            for sstring in self.line_species:
                tot_str = tot_str + sstring + ':'

            self.line_grid_kappas = fi.read_in_molecular_opacities(self.path,tot_str,freq_len_full,self.g_len, \
                                        len(self.line_species),len(self.line_TP_grid[:,0]),self.mode, arr_min, arr_max)
            if self.mode == 'c-k':
                self.line_grid_kappas = np.array(self.line_grid_kappas[:,index,:,:],dtype='d',order='Fortran')
            else:
                self.line_grid_kappas = np.array(self.line_grid_kappas,dtype='d',order='Fortran')
            
        # Read in g grid for correlated-k
        if self.mode == 'c-k':
            buffer = np.genfromtxt(self.path+'/opa_input_files/g_comb_grid.dat')
            self.g_gauss, self.w_gauss = buffer[:,0], buffer[:,1]
            self.g_gauss,self.w_gauss = np.array(self.g_gauss,dtype='d',order='Fortran'),np.array(self.w_gauss,dtype='d',order='Fortran')
        elif self.mode == 'lbl':
            self.g_gauss, self.w_gauss = np.ones(1), np.ones(1)
            self.g_gauss,self.w_gauss = np.array(self.g_gauss,dtype='d',order='Fortran'),np.array(self.w_gauss,dtype='d',order='Fortran')
        
        # Read in the angle (mu) grid for the emission spectral calculations.
        buffer = np.genfromtxt(self.path+'/opa_input_files/mu_points.dat')
        self.mu, self.w_gauss_mu = buffer[:,0], buffer[:,1]

        ###########################
        # Read in continuum opacities
        ###########################

        # Clouds
        if len(self.cloud_species) > 0:
            self.read_cloud_opas()

        # CIA
        if H2H2CIA:
          print('  Read CIA opacities for H2-H2...')
          self.cia_h2h2_lambda, self.cia_h2h2_temp, self.cia_h2h2_alpha_grid = fi.cia_read('H2H2',self.path)
          self.cia_h2h2_alpha_grid = np.array(self.cia_h2h2_alpha_grid,dtype='d',order='Fortran')
        if H2HeCIA:
          print('  Read CIA opacities for H2-He...')
          self.cia_h2he_lambda, self.cia_h2he_temp, self.cia_h2he_alpha_grid = fi.cia_read('H2He',self.path)
          self.cia_h2he_alpha_grid = np.array(self.cia_h2he_alpha_grid,dtype='d',order='Fortran')
        if H2H2CIA or H2HeCIA:
          print(' Done.')
          print()

    def calc_borders(self,x):
        ''' Return bin borders for midpoints. '''
        xn = []
        xn.append(x[0]-(x[1]-x[0])/2.)
        for i in range(int(len(x))-1):
            xn.append(x[i]+(x[i+1]-x[i])/2.)
        xn.append(x[int(len(x))-1]+(x[int(len(x))-1]-x[int(len(x))-2])/2.)
        return np.array(xn)

    def read_cloud_opas(self):
        ''' Function to read cloud opacities '''
        self.cloud_species_mode = []
        for i in range(int(len(self.cloud_species))):
            splitstr = self.cloud_species[i].split('_')
            self.cloud_species_mode.append(splitstr[1])
            self.cloud_species[i] = splitstr[0]

        # Prepare single strings delimited by ':' which are then
        # put into Fortran routines
        tot_str_names = ''
        for sstring in self.cloud_species:
            tot_str_names = tot_str_names + sstring + ':'

        tot_str_modes = ''
        for sstring in self.cloud_species_mode:
            tot_str_modes = tot_str_modes + sstring + ':'

        self.N_cloud_lambda_bins = int(len(np.genfromtxt(self.path + \
                    '/opacities/continuum/clouds/MgSiO3_c/amorphous/mie/opa_0001.dat')[:,0]))

        # Actual reading of opacities
        rho_cloud_particles, cloud_specs_abs_opa, cloud_specs_scat_opa, cloud_aniso, \
          cloud_lambdas, cloud_rad_bins, cloud_radii \
          = fi.read_in_cloud_opacities(self.path,tot_str_names,tot_str_modes, \
                            len(self.cloud_species),self.N_cloud_lambda_bins)

        self.rho_cloud_particles = np.array(rho_cloud_particles,dtype='d',order='Fortran')
        self.cloud_specs_abs_opa = np.array(cloud_specs_abs_opa,dtype='d',order='Fortran')
        self.cloud_specs_scat_opa = np.array(cloud_specs_scat_opa,dtype='d',order='Fortran')
        self.cloud_aniso = np.array(cloud_aniso,dtype='d',order='Fortran')
        self.cloud_lambdas = np.array(cloud_lambdas,dtype='d',order='Fortran')
        self.cloud_rad_bins = np.array(cloud_rad_bins,dtype='d',order='Fortran')
        self.cloud_radii = np.array(cloud_radii,dtype='d',order='Fortran')
        
    # Preparing structures
    def setup_opa_structure(self,P):
        ''' Setup opacity array for atmospheruc structure dimensions, as well as pressure array '''
        # bar to cgs
        self.press = P*1e6
        self.continuum_opa = np.array(np.zeros(self.freq_len*len(P)).reshape(self.freq_len, \
                        len(P)),dtype='d',order='Fortran')
        self.continuum_opa_scat = np.array(np.zeros(self.freq_len*len(P)).reshape(self.freq_len, \
                        len(P)),dtype='d',order='Fortran')
        self.contr_em = np.array(np.zeros(self.freq_len*len(P)).reshape(len(P),self.freq_len), \
                                     dtype='d',order='Fortran')
        self.contr_tr = np.array(np.zeros(self.freq_len*len(P)).reshape(len(P),self.freq_len), \
                                     dtype='d',order='Fortran')
        if len(self.line_species) > 0:
            self.line_struc_kappas = np.array(np.zeros(self.g_len*self.freq_len*len(P)* \
                len(self.line_species)).reshape(self.g_len,self.freq_len, \
                len(self.line_species),len(P)),dtype='d',order='Fortran')
            self.total_tau = np.array(np.zeros_like(self.line_struc_kappas),dtype='d',order='Fortran')
            self.line_abundances = np.array(np.zeros(len(self.press)*len(self.line_species)).reshape(len(self.press), \
                                                                    len(self.line_species)),dtype='d',order='Fortran')
        else: # If there are no specified line species then we need at least an array to contain the continuum opas
              # I'll (mis)use the line_struc_kappas array for that
            self.line_struc_kappas = np.array(np.zeros(self.g_len*self.freq_len*len(P)).reshape(self.g_len,self.freq_len, \
                1,len(P)),dtype='d',order='Fortran')
            self.total_tau = np.array(np.zeros_like(self.line_struc_kappas),dtype='d',order='Fortran')
            self.line_abundances = np.array(np.zeros(len(self.press)).reshape(len(self.press),1),dtype='d',order='Fortran')
        
        self.mmw = np.zeros_like(self.press)

        if len(self.cloud_species) > 0:
            self.cloud_mass_fracs = np.array(np.zeros(int(len(self.press))* \
                        int(len(self.cloud_species))).reshape(int(len(self.press)), \
                        int(len(self.cloud_species))),dtype='d',order='Fortran')
            self.r_g = np.array(np.zeros(int(len(self.press))* \
                        int(len(self.cloud_species))).reshape(int(len(self.press)), \
                        int(len(self.cloud_species))),dtype='d',order='Fortran')

    def interpolate_species_opa(self, temp):
        ''' Interpolate line opacities to given temperature structure. '''
        self.temp = temp
        if len(self.line_species) > 0:
            self.line_struc_kappas = fi.interpol_opa_ck(self.press,temp,self.line_TP_grid,self.line_grid_kappas)
        else:
            self.line_struc_kappas = np.zeros_like(self.line_struc_kappas)
            
    def mix_opa_tot(self, abundances, mmw, gravity,\
                        sigma_lnorm = None, fsed = None, Kzz = None, \
                        radius = None, gray_opacity = None):
        ''' Combine total line opacities, according to mass fractions (abundances),
            also add continuum opacities, i.e. clouds, CIA...'''
        
        self.scat = False
        self.mmw = mmw
        for i_spec in range(len(self.line_species)):
            self.line_abundances[:,i_spec] = abundances[self.line_species[i_spec]]
        self.continuum_opa = np.zeros_like(self.continuum_opa)
        self.continuum_opa_scat = np.zeros_like(self.continuum_opa_scat)

        # Calc. CIA opacity
        if self.H2H2CIA:
            self.continuum_opa = self.continuum_opa + fi.cia_interpol(self.freq,self.temp, \
                self.cia_h2h2_lambda,self.cia_h2h2_temp,self.cia_h2h2_alpha_grid, \
                    self.press,self.mmw,abundances['H2'],2.)
        if self.H2HeCIA:
            self.continuum_opa = self.continuum_opa + fi.cia_interpol(self.freq,self.temp, \
                self.cia_h2he_lambda,self.cia_h2he_temp,self.cia_h2he_alpha_grid, \
                self.press,self.mmw,np.sqrt(abundances['H2']*abundances['He']),np.sqrt(8.))

        # Add mock gray cloud opacity here
        if self.gray_opacity != None:
            self.continuum_opa = self.continuum_opa + self.gray_opacity

        # Add cloud opacity here, will modify self.continuum_opa
        if int(len(self.cloud_species)) > 0:
            self.calc_cloud_opacity(abundances, mmw, gravity, sigma_lnorm, fsed, Kzz, radius)

        # Interpolate line opacities, combine with continuum oacities
        self.line_struc_kappas = fi.mix_opas_ck(self.line_abundances, \
                                        self.line_struc_kappas,self.continuum_opa)
        if len(self.rayleigh_species) != 0:
            self.scat = True
            self.add_rayleigh(abundances)
        if (self.Pcloud != None):
            self.continuum_opa_scat[:,self.press>self.Pcloud*1e6] = 1e99

        # In the line-by-line case we can simply add the opacities of different species
        # in frequency space. All opacities are stored in the first species index slot
        if (self.mode == 'lbl') and (int(len(self.line_species)) > 1):
            self.line_struc_kappas[:,:,0,:] = np.sum(self.line_struc_kappas, axis = 2)

    def calc_cloud_opacity(self,abundances, mmw, gravity, sigma_lnorm, fsed = None, Kzz = None, radius = None):
        ''' Function to calculate cloud opacities for defined atmospheric structure. '''
        rho = self.press/nc.kB/self.temp*mmw*nc.amu
        for i_spec in range(int(len(self.cloud_species))):
            self.cloud_mass_fracs[:,i_spec] = abundances[self.cloud_species[i_spec]]
            if radius != None:
                self.r_g[:,i_spec] = radius[self.cloud_species[i_spec]]
        
        if radius != None:
            cloud_abs_opa_TOT,cloud_scat_opa_TOT,cloud_red_fac_aniso_TOT = \
              fs.calc_cloud_opas(rho,self.rho_cloud_particles,self.cloud_mass_fracs,self.r_g,sigma_lnorm, \
                                self.cloud_rad_bins,self.cloud_radii,self.cloud_lambdas, \
                                self.cloud_specs_abs_opa,self.cloud_specs_scat_opa, \
                                self.cloud_aniso)
        else:
            self.r_g = fs.get_rg_n(gravity,rho,self.rho_cloud_particles,self.temp,mmw,fsed,self.cloud_mass_fracs, \
                                    sigma_lnorm,Kzz)
            cloud_abs_opa_TOT,cloud_scat_opa_TOT,cloud_red_fac_aniso_TOT = \
              fs.calc_cloud_opas(rho,self.rho_cloud_particles,self.cloud_mass_fracs,self.r_g,sigma_lnorm, \
                                self.cloud_rad_bins,self.cloud_radii,self.cloud_lambdas, \
                                self.cloud_specs_abs_opa,self.cloud_specs_scat_opa, \
                                self.cloud_aniso)

        cloud_abs, cloud_scat, aniso, cloud_abs_tot_no_aniso = \
           fs.interp_integ_cloud_opas(cloud_abs_opa_TOT,cloud_scat_opa_TOT, \
            cloud_red_fac_aniso_TOT,self.cloud_lambdas,self.border_freqs)

        self.continuum_opa += cloud_abs
        self.continuum_opa_scat += cloud_abs_tot_no_aniso - cloud_abs
            
        return
    

    def add_rayleigh(self,abundances):
        ''' Add Rayleigh scattering cross-sections'''
        for spec in self.rayleigh_species:
            haze_multiply = 1.
            if (self.haze_factor != None):
                haze_multiply = self.haze_factor
            self.continuum_opa_scat = self.continuum_opa_scat + \
              haze_multiply*fs.add_rayleigh(spec,abundances[spec],self.lambda_angstroem, \
                                  self.mmw,self.temp,self.press)
        #print(self.continuum_opa_scat[:,0])

    def calc_opt_depth(self,gravity):
        ''' Calculate optical depth for the total opacity. '''
        if (self.mode == 'lbl') and (int(len(self.line_species)) > 1):
            self.total_tau[:,:,:1,:] = fs.calc_tau_g_tot_ck(gravity,self.press,self.line_struc_kappas[:,:,:1,:])
        else:
            self.total_tau = fs.calc_tau_g_tot_ck(gravity,self.press,self.line_struc_kappas)
        
    def calc_RT(self,contribution):
        ''' Calculate the flux'''
        if (self.mode == 'lbl') and (int(len(self.line_species)) > 1):
            self.flux, self.contr_em = fs.flux_ck(self.freq,self.total_tau[:,:,:1,:],self.temp, \
                self.mu,self.w_gauss_mu,self.w_gauss,contribution)
        else:
            self.flux, self.contr_em = fs.flux_ck(self.freq,self.total_tau,self.temp, \
                self.mu,self.w_gauss_mu,self.w_gauss,contribution)

    def calc_tr_rad(self,P0_bar,R_pl,gravity,mmw,contribution):
        ''' Calculate the transmission spectrum '''
        if (self.mode == 'lbl') and (int(len(self.line_species)) > 1):
            self.transm_rad = fs.calc_transm_spec(self.freq,self.line_struc_kappas[:,:,:1,:],self.temp, \
                                    self.press,gravity,mmw,P0_bar,R_pl,self.w_gauss,self.scat, \
                                                      self.continuum_opa_scat)
            if contribution:
                self.contr_tr = fs.calc_transm_spec_contr(self.freq,self.line_struc_kappas[:,:,:1,:], self.temp, \
                                    self.press,gravity,mmw,P0_bar,R_pl,self.w_gauss,self.transm_rad**2.,self.scat, \
                                    self.continuum_opa_scat)
            
        else:
            self.transm_rad = fs.calc_transm_spec(self.freq,self.line_struc_kappas,self.temp, \
                                    self.press,gravity,mmw,P0_bar,R_pl,self.w_gauss,self.scat, \
                                                      self.continuum_opa_scat)
            if contribution:
                self.contr_tr = fs.calc_transm_spec_contr(self.freq,self.line_struc_kappas,self.temp, \
                                    self.press,gravity,mmw,P0_bar,R_pl,self.w_gauss,self.transm_rad**2.,self.scat, \
                                    self.continuum_opa_scat)
        
    def calc_flux(self,temp,abunds,gravity,mmw,sigma_lnorm = None, \
                      fsed = None, Kzz = None, radius = None,contribution=False, \
                      gray_opacity = None):
        ''' Function to calc flux, called from outside '''
        self.gray_opacity = gray_opacity
        self.interpolate_species_opa(temp)
        self.mix_opa_tot(abunds,mmw,gravity,sigma_lnorm,fsed,Kzz,radius)
        self.calc_opt_depth(gravity)
        self.calc_RT(contribution)

    def calc_transm(self,temp,abunds,gravity,mmw,P0_bar,R_pl,sigma_lnorm = None, \
                        fsed = None, Kzz = None, radius = None,Pcloud=None, \
                        contribution=False,haze_factor=None, \
                        gray_opacity = None):
        ''' Function to calc transm. spectrum, called from outside '''
        self.Pcloud = Pcloud
        self.gray_opacity = gray_opacity
        self.interpolate_species_opa(temp)
        self.haze_factor = haze_factor
        self.mix_opa_tot(abunds,mmw,gravity,sigma_lnorm,fsed,Kzz,radius)
        self.calc_tr_rad(P0_bar,R_pl,gravity,mmw,contribution)
        

    def calc_flux_transm(self,temp,abunds,gravity,mmw,P0_bar,R_pl,sigma_lnorm = None, \
                             fsed = None, Kzz = None, radius = None,Pcloud=None, \
                             contribution=False,gray_opacity = None):
        ''' Function to calc flux and transmission spectrum, called from outside '''
        self.Pcloud = Pcloud
        self.gray_opacity = gray_opacity
        self.interpolate_species_opa(temp)
        self.mix_opa_tot(abunds,mmw,gravity,sigma_lnorm,fsed,Kzz,radius)
        self.calc_opt_depth(gravity)
        self.calc_RT(contribution)
        self.calc_tr_rad(P0_bar,R_pl,gravity,mmw,contribution)

########################################################
### Radtrans utility for temperature model computation
########################################################

### Box car conv. average
def box_car_conv(array,points):

    res = np.zeros_like(array)
    len_arr = int(len(array) + 0.01)
    for i in range(len(array)):
        if (i-int(points/2) >= 0) and (i+int(points/2) <= len_arr+1):
            smooth_val = array[i-int(points/2):i+int(points/2)]
            res[i] = np.sum(smooth_val)/len(smooth_val)
        elif (i+int(points/2) > len_arr+1):
            len_use = len_arr+1-i
            smooth_val = array[i-len_use:i+len_use]
            res[i] = np.sum(smooth_val)/len(smooth_val)
        elif i-int(points/2) < 0:
            smooth_val = array[:max(2*i,1)]
            res[i] = np.sum(smooth_val)/len(smooth_val)
    return res


### Global Guillot P-T formula with kappa/grav replaced by delta
def guillot_global(P,delta,gamma,T_int,T_equ):

    delta = np.abs(delta)
    gamma = np.abs(gamma)
    T_int = np.abs(T_int)
    T_equ = np.abs(T_equ)
    tau = P*1e6*delta
    T_irr = T_equ*np.sqrt(2.)
    T = (0.75*T_int**4.*(2./3.+tau) + \
      0.75*T_irr**4./4.*(2./3.+1./gamma/3.**0.5+ \
                         (gamma/3.**0.5-1./3.**0.5/gamma)*np.exp(-gamma*tau*3.**0.5)))**0.25
    return T

### Modified Guillot P-T formula
def guillot_modif(P,delta,gamma,T_int,T_equ,ptrans,alpha):
    return guillot_global(P,np.abs(delta),np.abs(gamma),np.abs(T_int),np.abs(T_equ))* \
      (1.-alpha*(1./(1.+np.exp((np.log(P/ptrans))))))

### Function to make temp
def make_press_temp(rad_trans_params):
    press_many = np.logspace(-6,5,210)
    index = (press_many <= 1e3) & (press_many >= 1e-6)
    press = press_many[index][::2]
    t = box_car_conv(guillot_modif(press_many, \
        1e1**rad_trans_params['log_delta'],1e1**rad_trans_params['log_gamma'], \
        rad_trans_params['t_int'],rad_trans_params['t_equ'], \
        1e1**rad_trans_params['log_p_trans'],rad_trans_params['alpha']),20)
    temp = t[index][::2]
    return press, temp