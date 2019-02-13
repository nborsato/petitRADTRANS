from __future__ import division, print_function

from . import fort_input as fi
from . import fort_spec as fs
from . import nat_cst as nc
from . import pyth_input as pyi

import numpy as np
import copy as cp
import os
import sys

class Radtrans:
    """ Class defining objects for carrying out spectral calculations for a
    given set of opacities

    Args:
        line_species (Optional):
            list of strings, denoting which line absorber species to include.
        rayleigh_species (Optional):
            list of strings, denoting which Rayleigh scattering species to 
            include.
        cloud_species (Optional): 
            list of strings, denoting which cloud opacity species to include.
        continuum_species (Optional):
            list of strings, denoting which continuum absorber species to
            include.
        H2H2CIA (Optional[bool]):
            Will be ``False`` by default.
            If ``True``, will add H2-H2 Collision induced
            absoprtion as continuum absorber (alternatively, put ``'H2-H2'``
            into continuum_species list).
        H2HeCIA (Optional[bool]):
            Will be ``False`` by default.
            If ``True``, will add H2-He Collision induced
            absoprtion as continuum absorber (alternatively, put ``'H2-He'``
            into continuum_species list).
        wlen_bords_micron (Optional):
            list containing left and right border of wavelength region to be
            considered, in micron. If nothing else is specified, it will be
            equal to ``[0.05, 300]``, hence using the full petitRADTRANS
            wavelength range (0.11 to 250 microns for ``'c-k'`` mode, 0.3 to 30
            microns for the ``'lbl'`` mode). The larger the range the longer the
            computation time.
        mode (Optinional[string]):
            if equal to ``'c-k'``: use low-resolution mode, at
            :math:`\lambda/\Delta \lambda = 1000`, with the correlated-k
            assumption. if equal to ``'lbl'``: use high-resolution mode, at
            :math:`\lambda/\Delta \lambda = 10^6`, with a line-by-line
            treatment.

    """
    
    def __init__(self, line_species=[], rayleigh_species=[], cloud_species=[], \
                     continuum_opacities = [], H2H2CIA=False, H2HeCIA=False, \
                     wlen_bords_micron=[0.05,300.], mode='c-k'):

        # Line-by-line or corr-k
        self.mode = mode

        # Line opacity species to be considered
        self.line_species = line_species

        # Rayleigh scattering species to be considered
        self.rayleigh_species = rayleigh_species

        # Cloud species to be considered
        self.cloud_species = cloud_species

        # Include continuum opacities?
        # Still allow for old way, when only CIA were continuum opacities
        self.H2H2CIA = H2H2CIA
        self.H2HeCIA = H2HeCIA
        # New species
        self.Hminus = False
        # Check what is supposed to be included.
        if len(continuum_opacities) > 0:
            for c in continuum_opacities:
                if c == 'H2-H2':
                    self.H2H2CIA = True
                elif c == 'H2-He':
                    self.H2HeCIA = True
                elif c == 'H-':
                    self.Hminus = True
            

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
        # Some necessary definitions, also prepare arrays for fluxes,
        # transmission radius...
        ###########################        

        self.lambda_angstroem = np.array(nc.c/self.freq/1e-8,dtype='d', \
                                             order='Fortran')
        self.flux = np.array(np.zeros(self.freq_len),dtype='d', \
                                 order='Fortran')
        self.transm_rad = np.array(np.zeros(self.freq_len),dtype='d', \
                                       order='Fortran')

        # Define frequency bins around grid for later interpolation
        # purposes when including
        # clouds...
        self.border_freqs = np.array(nc.c/self.calc_borders(nc.c/self.freq), \
                                         dtype='d',order='Fortran')
        self.border_lambda_angstroem = \
          np.array(self.calc_borders(self.lambda_angstroem))

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
        self.line_TP_grid = np.array(self.line_TP_grid.reshape( \
                    len(self.line_TP_grid[:,1]),2),dtype='d',order='Fortran')

        # Read actual opacities....
        # line_grid_kappas has the shape g_len,freq_len,len(line_species),
        # len(line_TP_grid[:,0])
        if len(self.line_species) > 0:

            tot_str = ''
            for sstring in self.line_species:
                tot_str = tot_str + sstring + ':'

            self.line_grid_kappas = fi.read_in_molecular_opacities( \
                    self.path,tot_str,freq_len_full,self.g_len, \
                    len(self.line_species),len(self.line_TP_grid[:,0]), \
                                                self.mode, arr_min, arr_max)
            if self.mode == 'c-k':
                self.line_grid_kappas = np.array( \
                    self.line_grid_kappas[:,index,:,:], \
                    dtype='d',order='Fortran')
            else:
                self.line_grid_kappas = \
                  np.array(self.line_grid_kappas,dtype='d',order='Fortran')
            
        # Read in g grid for correlated-k
        if self.mode == 'c-k':
            buffer = np.genfromtxt(self.path+'/opa_input_files/g_comb_grid.dat')
            self.g_gauss, self.w_gauss = buffer[:,0], buffer[:,1]
            self.g_gauss,self.w_gauss = np.array(self.g_gauss,dtype='d', \
                order='Fortran'),np.array(self.w_gauss, \
                dtype='d',order='Fortran')
        elif self.mode == 'lbl':
            self.g_gauss, self.w_gauss = np.ones(1), np.ones(1)
            self.g_gauss,self.w_gauss = np.array(self.g_gauss,dtype='d', \
                order='Fortran'),np.array(self.w_gauss, \
                dtype='d',order='Fortran')
        
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
        if self.H2H2CIA:
          print('  Read CIA opacities for H2-H2...')
          self.cia_h2h2_lambda, self.cia_h2h2_temp, \
            self.cia_h2h2_alpha_grid = fi.cia_read('H2H2',self.path)
          self.cia_h2h2_alpha_grid = np.array(self.cia_h2h2_alpha_grid, \
                                                dtype='d',order='Fortran')
        if self.H2HeCIA:
          print('  Read CIA opacities for H2-He...')
          self.cia_h2he_lambda, self.cia_h2he_temp, self.cia_h2he_alpha_grid = \
            fi.cia_read('H2He',self.path)
          self.cia_h2he_alpha_grid = np.array(self.cia_h2he_alpha_grid, \
                                                dtype='d',order='Fortran')
        if self.H2H2CIA or self.H2HeCIA:
          print(' Done.')
          print()

    def calc_borders(self,x):
        # Return bin borders for midpoints.
        xn = []
        xn.append(x[0]-(x[1]-x[0])/2.)
        for i in range(int(len(x))-1):
            xn.append(x[i]+(x[i+1]-x[i])/2.)
        xn.append(x[int(len(x))-1]+(x[int(len(x))-1]-x[int(len(x))-2])/2.)
        return np.array(xn)

    def read_cloud_opas(self):
        # Function to read cloud opacities
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
            '/opacities/continuum/clouds/MgSiO3_c/amorphous/mie/opa_0001.dat' \
                                                             )[:,0]))

        # Actual reading of opacities
        rho_cloud_particles, cloud_specs_abs_opa, cloud_specs_scat_opa, \
          cloud_aniso, cloud_lambdas, cloud_rad_bins, cloud_radii \
          = fi.read_in_cloud_opacities(self.path,tot_str_names,tot_str_modes, \
                            len(self.cloud_species),self.N_cloud_lambda_bins)

        self.rho_cloud_particles = \
          np.array(rho_cloud_particles,dtype='d',order='Fortran')
        self.cloud_specs_abs_opa = \
          np.array(cloud_specs_abs_opa,dtype='d',order='Fortran')
        self.cloud_specs_scat_opa = \
          np.array(cloud_specs_scat_opa,dtype='d',order='Fortran')
        self.cloud_aniso = np.array(cloud_aniso,dtype='d',order='Fortran')
        self.cloud_lambdas = np.array(cloud_lambdas,dtype='d',order='Fortran')
        self.cloud_rad_bins = np.array(cloud_rad_bins,dtype='d',order='Fortran')
        self.cloud_radii = np.array(cloud_radii,dtype='d',order='Fortran')
        
    # Preparing structures
    def setup_opa_structure(self,P):
        ''' Setup opacity arrays at atmospheric structure dimensions, 
        and set the atmospheric pressure array.

        Args:
            P:
                the atmospheric pressure (1-d numpy array, sorted in increasing
                order), in units of bar. Will be converted to cgs internally.  
        '''
        
        # bar to cgs
        self.press = P*1e6
        self.continuum_opa = np.array(np.zeros(self.freq_len*len(P)) \
                                          .reshape(self.freq_len, \
                                          len(P)),dtype='d',order='Fortran')
        self.continuum_opa_scat = np.array(np.zeros(self.freq_len*len(P)) \
                                            .reshape(self.freq_len, \
                                            len(P)),dtype='d',order='Fortran')
        self.contr_em = np.array(np.zeros(self.freq_len*len(P)) \
                                     .reshape(len(P),self.freq_len), \
                                     dtype='d',order='Fortran')
        self.contr_tr = np.array(np.zeros(self.freq_len*len(P)) \
                                     .reshape(len(P),self.freq_len), \
                                     dtype='d',order='Fortran')
        if len(self.line_species) > 0:
            self.line_struc_kappas = \
              np.array(np.zeros(self.g_len*self.freq_len*len(P)* \
                len(self.line_species)).reshape(self.g_len,self.freq_len, \
                len(self.line_species),len(P)),dtype='d',order='Fortran')
            self.total_tau = \
              np.array(np.zeros_like(self.line_struc_kappas), \
                           dtype='d',order='Fortran')
            self.line_abundances = \
              np.array(np.zeros(len(self.press)*len(self.line_species)) \
                           .reshape(len(self.press), \
                            len(self.line_species)),dtype='d',order='Fortran')
        else: # If there are no specified line species then we need at
              # least an array to contain the continuum opas
              # I'll (mis)use the line_struc_kappas array for that
            self.line_struc_kappas = \
              np.array(np.zeros(self.g_len*self.freq_len*len(P)) \
                           .reshape(self.g_len,self.freq_len, \
                            1,len(P)),dtype='d',order='Fortran')
            self.total_tau = \
              np.array(np.zeros_like(self.line_struc_kappas), \
                           dtype='d',order='Fortran')
            self.line_abundances = \
              np.array(np.zeros(len(self.press)) \
                           .reshape(len(self.press),1), \
                           dtype='d',order='Fortran')
        
        self.mmw = np.zeros_like(self.press)

        if len(self.cloud_species) > 0:
            self.cloud_mass_fracs = np.array(np.zeros(int(len(self.press))* \
                        int(len(self.cloud_species))) \
                            .reshape(int(len(self.press)), \
                            int(len(self.cloud_species))), \
                                                 dtype='d',order='Fortran')
            self.r_g = np.array(np.zeros(int(len(self.press))* \
                        int(len(self.cloud_species))) \
                                    .reshape(int(len(self.press)), \
                                            int(len(self.cloud_species))), \
                                            dtype='d',order='Fortran')

    def interpolate_species_opa(self, temp):
        # Interpolate line opacities to given temperature structure.
        self.temp = temp
        if len(self.line_species) > 0:
            self.line_struc_kappas = fi.interpol_opa_ck(self.press,temp, \
                                    self.line_TP_grid,self.line_grid_kappas)
        else:
            self.line_struc_kappas = np.zeros_like(self.line_struc_kappas)
            
    def mix_opa_tot(self, abundances, mmw, gravity, \
                        sigma_lnorm = None, fsed = None, Kzz = None, \
                        radius = None, gray_opacity = None, \
                        add_cloud_scat_as_abs = None):
        # Combine total line opacities,
        # according to mass fractions (abundances),
        # also add continuum opacities, i.e. clouds, CIA...
        
        self.scat = False
        self.mmw = mmw
        for i_spec in range(len(self.line_species)):
            self.line_abundances[:,i_spec] = abundances[self.line_species[i_spec]]
        self.continuum_opa = np.zeros_like(self.continuum_opa)
        self.continuum_opa_scat = np.zeros_like(self.continuum_opa_scat)

        # Calc. CIA opacity
        if self.H2H2CIA:
            self.continuum_opa = self.continuum_opa + \
              fi.cia_interpol(self.freq,self.temp, \
                self.cia_h2h2_lambda,self.cia_h2h2_temp, \
                self.cia_h2h2_alpha_grid, \
                self.press,self.mmw,abundances['H2'],2.)
        if self.H2HeCIA:
            self.continuum_opa = self.continuum_opa + \
              fi.cia_interpol(self.freq,self.temp, \
                self.cia_h2he_lambda,self.cia_h2he_temp, \
                self.cia_h2he_alpha_grid, \
                self.press,self.mmw, \
                np.sqrt(abundances['H2']*abundances['He']),np.sqrt(8.))

        # Calc. H- opacity
        if self.Hminus:
            self.continuum_opa = \
              self.continuum_opa + pyi.hminus_opacity(self.lambda_angstroem, \
                self.border_lambda_angstroem, \
                self.temp, self.press, mmw, abundances)

        # Add mock gray cloud opacity here
        if self.gray_opacity != None:
            self.continuum_opa = self.continuum_opa + self.gray_opacity

        # Add cloud opacity here, will modify self.continuum_opa
        if int(len(self.cloud_species)) > 0:
            self.calc_cloud_opacity(abundances, mmw, gravity, \
                                        sigma_lnorm, fsed, Kzz, radius, \
                                        add_cloud_scat_as_abs)

        if len(self.rayleigh_species) != 0:
            self.scat = True
            self.add_rayleigh(abundances)
        if (self.Pcloud != None):
            self.continuum_opa[:,self.press>self.Pcloud*1e6] += 1e99

        # Interpolate line opacities, combine with continuum oacities
        self.line_struc_kappas = fi.mix_opas_ck(self.line_abundances, \
                                    self.line_struc_kappas,self.continuum_opa)

        # In the line-by-line case we can simply
        # add the opacities of different species
        # in frequency space. All opacities are
        # stored in the first species index slot
        if (self.mode == 'lbl') and (int(len(self.line_species)) > 1):
            self.line_struc_kappas[:,:,0,:] = \
              np.sum(self.line_struc_kappas, axis = 2)

    def calc_cloud_opacity(self,abundances, mmw, gravity, sigma_lnorm, \
                               fsed = None, Kzz = None, \
                               radius = None, add_cloud_scat_as_abs = None):
        # Function to calculate cloud opacities
        # for defined atmospheric structure.
        rho = self.press/nc.kB/self.temp*mmw*nc.amu
        for i_spec in range(int(len(self.cloud_species))):
            self.cloud_mass_fracs[:,i_spec] = \
              abundances[self.cloud_species[i_spec]]
            if radius != None:
                self.r_g[:,i_spec] = radius[self.cloud_species[i_spec]]
        
        if radius != None:
            cloud_abs_opa_TOT,cloud_scat_opa_TOT,cloud_red_fac_aniso_TOT = \
              fs.calc_cloud_opas(rho,self.rho_cloud_particles, \
                                self.cloud_mass_fracs,self.r_g,sigma_lnorm, \
                                self.cloud_rad_bins,self.cloud_radii, \
                                self.cloud_lambdas, \
                                self.cloud_specs_abs_opa, \
                                self.cloud_specs_scat_opa, \
                                self.cloud_aniso)
        else:
            self.r_g = fs.get_rg_n(gravity,rho,self.rho_cloud_particles, \
                                       self.temp,mmw,fsed, \
                                       self.cloud_mass_fracs, \
                                       sigma_lnorm,Kzz)
            cloud_abs_opa_TOT,cloud_scat_opa_TOT,cloud_red_fac_aniso_TOT = \
              fs.calc_cloud_opas(rho,self.rho_cloud_particles, \
                                     self.cloud_mass_fracs, \
                                     self.r_g,sigma_lnorm, \
                                     self.cloud_rad_bins,self.cloud_radii, \
                                     self.cloud_lambdas, \
                                     self.cloud_specs_abs_opa, \
                                     self.cloud_specs_scat_opa, \
                                     self.cloud_aniso)

        cloud_abs, cloud_scat, aniso, cloud_abs_tot_no_aniso = \
           fs.interp_integ_cloud_opas(cloud_abs_opa_TOT,cloud_scat_opa_TOT, \
            cloud_red_fac_aniso_TOT,self.cloud_lambdas,self.border_freqs)

        self.continuum_opa_scat += cloud_abs_tot_no_aniso - cloud_abs

        if add_cloud_scat_as_abs != None:
            if add_cloud_scat_as_abs:
                self.continuum_opa += cloud_abs \
                  + 0.20*(cloud_abs_tot_no_aniso - cloud_abs)
            else:
                self.continuum_opa += cloud_abs
        else:
            self.continuum_opa += cloud_abs
            
        return
    

    def add_rayleigh(self,abundances):
        # Add Rayleigh scattering cross-sections
        for spec in self.rayleigh_species:
            haze_multiply = 1.
            if (self.haze_factor != None):
                haze_multiply = self.haze_factor
            self.continuum_opa_scat = self.continuum_opa_scat + \
              haze_multiply*fs.add_rayleigh(spec,abundances[spec], \
                                                self.lambda_angstroem, \
                                                self.mmw,self.temp,self.press)
                                                
    def calc_opt_depth(self,gravity):
        # Calculate optical depth for the total opacity.
        if (self.mode == 'lbl') and (int(len(self.line_species)) > 1):
            self.total_tau[:,:,:1,:] = fs.calc_tau_g_tot_ck(gravity, \
                                self.press,self.line_struc_kappas[:,:,:1,:])
        else:
            self.total_tau = fs.calc_tau_g_tot_ck(gravity,self.press, \
                                                      self.line_struc_kappas)
        
    def calc_RT(self,contribution):
        # Calculate the flux
        if (self.mode == 'lbl') and (int(len(self.line_species)) > 1):
            self.flux, self.contr_em = fs.flux_ck(self.freq, \
                                                  self.total_tau[:,:,:1,:], \
                                                  self.temp, \
                                                  self.mu, \
                                                  self.w_gauss_mu,
                                                  self.w_gauss, \
                                                  contribution)
        else:
            self.flux, self.contr_em = fs.flux_ck(self.freq, \
                             self.total_tau,self.temp, \
                             self.mu,self.w_gauss_mu, \
                             self.w_gauss,contribution)

    def calc_tr_rad(self,P0_bar,R_pl,gravity,mmw, \
                        contribution,variable_gravity):
        # Calculate the transmission spectrum
        if (self.mode == 'lbl') and (int(len(self.line_species)) > 1):
            self.transm_rad = fs.calc_transm_spec(self.freq, \
                                self.line_struc_kappas[:,:,:1,:],self.temp, \
                                self.press,gravity,mmw,P0_bar,R_pl, \
                                self.w_gauss,self.scat, \
                                self.continuum_opa_scat,variable_gravity)
            if contribution:
                self.contr_tr = fs.calc_transm_spec_contr(self.freq, \
                                self.line_struc_kappas[:,:,:1,:], self.temp, \
                                self.press,gravity,mmw,P0_bar,R_pl, \
                                self.w_gauss,self.transm_rad**2.,self.scat, \
                                self.continuum_opa_scat,variable_gravity)
            
        else:
            self.transm_rad = fs.calc_transm_spec(self.freq, \
                                    self.line_struc_kappas,self.temp, \
                                    self.press,gravity,mmw,P0_bar,R_pl, \
                                    self.w_gauss,self.scat, \
                                    self.continuum_opa_scat,variable_gravity)
            if contribution:
                self.contr_tr = fs.calc_transm_spec_contr(self.freq, \
                                    self.line_struc_kappas,self.temp, \
                                    self.press,gravity,mmw,P0_bar,R_pl, \
                                    self.w_gauss,self.transm_rad**2., \
                                    self.scat, \
                                    self.continuum_opa_scat,variable_gravity)
        
    def calc_flux(self,temp,abunds,gravity,mmw,sigma_lnorm = None, \
                      fsed = None, Kzz = None, radius = None, \
                      contribution=False, \
                      gray_opacity = None, add_cloud_scat_as_abs = None):
        ''' Method to calculate the atmosphere's emitted flux
        (emission spectrum).

            Args:
                temp:
                    the atmospheric temperature in K, at each atmospheric layer
                    (1-d numpy array, same length as pressure array).
                abunds:
                    dictionary of mass fractions for all atmospheric absorbers.
                    Dictionary keys are the species names.
                    Every mass fraction array
                    has same length as pressure array.
                gravity (float):
                    Surface gravity in cgs. Vertically constant for emission
                    spectra.
                mmw:
                    the atmospheric mean molecular weight in amu,
                    at each atmospheric layer
                    (1-d numpy array, same length as pressure array).
                sigma_lnorm (Optional[float]):
                    width of the log-normal cloud particle size distribution
                fsed (Optional[float]):
                    cloud settling parameter
                Kzz (Optional):
                    the atmospheric eddy diffusion coeffiecient in cgs untis
                    (i.e. :math:`\\rm cm^2/s`),
                    at each atmospheric layer
                    (1-d numpy array, same length as pressure array).
                radius (Optional):
                    dictionary of mean particle radii for all cloud species.
                    Dictionary keys are the cloud species names.
                    Every radius array has same length as pressure array.    
                contribution (Optional[bool]):
                    If ``True`` the emission contribution function will be
                    calculated. Default is ``False``.
                gray_opacity (Optional[float]):
                    Gray opacity value, to be added to the opacity at all
                    pressures and wavelengths (units :math:`\\rm cm^2/g`)
                add_cloud_scat_as_abs (Optional[bool]):
                    If ``True``, 20 % of the cloud scattering opacity will be
                    added to the absorption opacity, introduced to test for the
                    effect of neglecting scattering.
        '''
        self.gray_opacity = gray_opacity
        self.interpolate_species_opa(temp)
        self.mix_opa_tot(abunds,mmw,gravity,sigma_lnorm,fsed,Kzz,radius, \
                             add_cloud_scat_as_abs = add_cloud_scat_as_abs)
        self.calc_opt_depth(gravity)
        self.calc_RT(contribution)

    def calc_transm(self,temp,abunds,gravity,mmw,P0_bar,R_pl, \
                        sigma_lnorm = None, \
                        fsed = None, Kzz = None, radius = None, \
                        Pcloud = None, \
                        contribution = False, haze_factor = None, \
                        gray_opacity = None,variable_gravity=True):
        ''' Method to calculate the atmosphere's transmission radius
        (for the transmission spectrum).

            Args:
                temp:
                    the atmospheric temperature in K, at each atmospheric layer
                    (1-d numpy array, same length as pressure array).
                abunds:
                    dictionary of mass fractions for all atmospheric absorbers.
                    Dictionary keys are the species names.
                    Every mass fraction array
                    has same length as pressure array.
                gravity (float):
                    Surface gravity in cgs at reference radius and pressure.
                mmw:
                    the atmospheric mean molecular weight in amu,
                    at each atmospheric layer
                    (1-d numpy array, same length as pressure array).
                P0_bar (float):
                    Reference pressure P0 in bar where R(P=P0) = R_pl,
                    where R_pl is the reference radius (parameter of this
                    method), and g(P=P0) = gravity, where gravity is the
                    reference gravity (parameter of this method)
                R_pl (float):
                    Reference radius R_pl, in cm.
                sigma_lnorm (Optional[float]):
                    width of the log-normal cloud particle size distribution
                fsed (Optional[float]):
                    cloud settling parameter
                Kzz (Optional):
                    the atmospheric eddy diffusion coeffiecient in cgs untis
                    (i.e. :math:`\\rm cm^2/s`),
                    at each atmospheric layer
                    (1-d numpy array, same length as pressure array).
                radius (Optional):
                    dictionary of mean particle radii for all cloud species.
                    Dictionary keys are the cloud species names.
                    Every radius array has same length as pressure array.    
                contribution (Optional[bool]):
                    If ``True`` the transmission and emission
                    contribution function will be
                    calculated. Default is ``False``.
                gray_opacity (Optional[float]):
                    Gray opacity value, to be added to the opacity at all
                    pressures and wavelengths (units :math:`\\rm cm^2/g`)
                Pcloud (Optional[float]):
                    Pressure, in bar, where opaque cloud deck is added to the
                    scattering opacity.
                haze_factor (Optional[float]):
                    Scalar factor, increasing the gas Rayleigh scattering
                    cross-section.
                variable_gravity (Optional[bool]):
                    Standard is ``True``. If ``False`` the gravity will be
                    constant as a function of pressure, during the transmission
                    radius calculation.
                add_cloud_scat_as_abs (Optional[bool]):
                    If ``True``, 20 % of the cloud scattering opacity will be
                    added to the absorption opacity, introduced to test for the
                    effect of neglecting scattering.
        '''
        
        self.Pcloud = Pcloud
        self.gray_opacity = gray_opacity
        self.interpolate_species_opa(temp)
        self.haze_factor = haze_factor
        self.mix_opa_tot(abunds,mmw,gravity,sigma_lnorm,fsed,Kzz,radius)
        self.calc_tr_rad(P0_bar,R_pl,gravity,mmw,contribution,variable_gravity)
        

    def calc_flux_transm(self,temp,abunds,gravity,mmw,P0_bar,R_pl,\
                             sigma_lnorm = None, \
                             fsed = None, Kzz = None, radius = None, \
                             Pcloud=None, \
                             contribution=False,gray_opacity = None, \
                             add_cloud_scat_as_abs = None, \
                             variable_gravity=True):
        ''' Method to calculate the atmosphere's emission flux *and*
        transmission radius (for the transmission spectrum).

            Args:
                temp:
                    the atmospheric temperature in K, at each atmospheric layer
                    (1-d numpy array, same length as pressure array).
                abunds:
                    dictionary of mass fractions for all atmospheric absorbers.
                    Dictionary keys are the species names.
                    Every mass fraction array
                    has same length as pressure array.
                gravity (float):
                    Surface gravity in cgs at reference radius and pressure,
                    constant durng the emission spectrum calculation.
                mmw:
                    the atmospheric mean molecular weight in amu,
                    at each atmospheric layer
                    (1-d numpy array, same length as pressure array).
                P0_bar (float):
                    Reference pressure P0 in bar where R(P=P0) = R_pl,
                    where R_pl is the reference radius (parameter of this
                    method), and g(P=P0) = gravity, where gravity is the
                    reference gravity (parameter of this method)
                R_pl (float):
                    Reference radius R_pl, in cm.
                sigma_lnorm (Optional[float]):
                    width of the log-normal cloud particle size distribution
                fsed (Optional[float]):
                    cloud settling parameter
                Kzz (Optional):
                    the atmospheric eddy diffusion coeffiecient in cgs untis
                    (i.e. :math:`\\rm cm^2/s`),
                    at each atmospheric layer
                    (1-d numpy array, same length as pressure array).
                radius (Optional):
                    dictionary of mean particle radii for all cloud species.
                    Dictionary keys are the cloud species names.
                    Every radius array has same length as pressure array.    
                contribution (Optional[bool]):
                    If ``True`` the transmission contribution function will be
                    calculated. Default is ``False``.
                gray_opacity (Optional[float]):
                    Gray opacity value, to be added to the opacity at all
                    pressures and wavelengths (units :math:`\\rm cm^2/g`)
                Pcloud (Optional[float]):
                    Pressure, in bar, where opaque cloud deck is added to the
                    scattering opacity.
                haze_factor (Optional[float]):
                    Scalar factor, increasing the gas Rayleigh scattering
                    cross-section.
                variable_gravity (Optional[bool]):
                    Standard is ``True``. If ``False`` the gravity will be
                    constant as a function of pressure, during the transmission
                    radius calculation.
        '''
        self.Pcloud = Pcloud
        self.gray_opacity = gray_opacity
        self.interpolate_species_opa(temp)
        self.mix_opa_tot(abunds,mmw,gravity,sigma_lnorm,fsed,Kzz,radius, \
                             add_cloud_scat_as_abs = add_cloud_scat_as_abs)
        self.calc_opt_depth(gravity)
        self.calc_RT(contribution)
        self.calc_tr_rad(P0_bar,R_pl,gravity,mmw,contribution,variable_gravity)

    def get_opa(self,temp):
        # Function to calc flux, called from outside
        self.interpolate_species_opa(temp)
        return self.line_struc_kappas

