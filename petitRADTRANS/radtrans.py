from __future__ import division, print_function

from . import fort_input as fi
from . import fort_spec as fs
from . import nat_cst as nc
from . import pyth_input as pyi

import numpy as np
import copy as cp
import os
import sys
import h5py
import glob

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
            self.freq_len_full = cp.copy(self.freq_len)
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
            self.freq_len_full = self.freq_len
            # Read in the frequency range of the opcity data
            wlen = fi.read_wlen(arr_min, arr_max, self.freq_len, path_length)
            self.freq = nc.c/wlen

        if self.mode == 'c-k':
            # To cut opacity data to required range
            index = (nc.c/self.freq > wlen_bords_micron[0]*1e-4) & \
              (nc.c/self.freq < wlen_bords_micron[1]*1e-4)

            self.freq_full = self.freq
            self.freq = np.array(self.freq[index],dtype='d',order='F')
            self.freq_len = len(self.freq)
        else:
            index = None

        # Read in the angle (mu) grid for the emission spectral calculations.
        buffer = np.genfromtxt(self.path+'/opa_input_files/mu_points.dat')
        self.mu, self.w_gauss_mu = buffer[:,0], buffer[:,1]

        ###########################
        # Some necessary definitions, also prepare arrays for fluxes,
        # transmission radius...
        ###########################        

        self.lambda_angstroem = np.array(nc.c/self.freq/1e-8,dtype='d', \
                                             order='F')
        self.flux = np.array(np.zeros(self.freq_len),dtype='d', \
                                 order='F')
        self.transm_rad = np.array(np.zeros(self.freq_len),dtype='d', \
                                       order='F')

        # Define frequency bins around grid for later interpolation
        # purposes when including
        # clouds...
        self.border_freqs = np.array(nc.c/self.calc_borders(nc.c/self.freq), \
                                         dtype='d',order='F')
        self.border_lambda_angstroem = \
          np.array(self.calc_borders(self.lambda_angstroem))

        self.Pcloud = None
        self.haze_factor = None
        self.gray_opacity = None

        ############################
        ############################
        # START Reading in opacities
        ############################
        ############################

        ###########################
        # Read in line opacities
        ###########################

        self.read_line_opacities(index, arr_min, arr_max)        

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
                                                dtype='d',order='F')
        if self.H2HeCIA:
          print('  Read CIA opacities for H2-He...')
          self.cia_h2he_lambda, self.cia_h2he_temp, self.cia_h2he_alpha_grid = \
            fi.cia_read('H2He',self.path)
          self.cia_h2he_alpha_grid = np.array(self.cia_h2he_alpha_grid, \
                                                dtype='d',order='F')
        if self.H2H2CIA or self.H2HeCIA:
          print(' Done.')
          print()

        #############################
        #############################
        # END Reading in opacities
        #############################
        #############################


    def calc_borders(self,x):
        # Return bin borders for midpoints.
        xn = []
        xn.append(x[0]-(x[1]-x[0])/2.)
        for i in range(int(len(x))-1):
            xn.append(x[i]+(x[i+1]-x[i])/2.)
        xn.append(x[int(len(x))-1]+(x[int(len(x))-1]-x[int(len(x))-2])/2.)
        return np.array(xn)

    def read_line_opacities(self, index, arr_min, arr_max):
        # Reads in the line opacities for spectral calculation

        # First get the P-Ts position where the grid is defined.
        # This here is for the nominal, log-uniform 10 x 13 point
        # P-T grid.
        buffer = np.genfromtxt(self.path+'/opa_input_files/opa_PT_grid.dat')
        self.line_TP_grid = np.zeros_like(buffer)
        self.line_TP_grid[:,0] = buffer[:,1]
        self.line_TP_grid[:,1] = buffer[:,0]
        # Convert from bars to cgs
        self.line_TP_grid[:,1] = 1e6*self.line_TP_grid[:,1]
        self.line_TP_grid = np.array(self.line_TP_grid.reshape( \
                    len(self.line_TP_grid[:,1]),2),dtype='d',order='F')

        # Check if species has custom P-T grid and reads in this grid.
        # Grid must be sorted appropriately, but the pyi.get_custom_grid()
        # will do that for the user in case a randomly ordered PTpaths.ls
        # is specified by the user in the opacity folder of the relevant species.
        # Only condition: it needs to be rectangular.
        # Because it is easier, a custom grid is saved for every species,
        # also the ones that use the nominal P-T grid of petitRADTRANS
        
        self.custom_grid = {}
        
        if len(self.line_species) > 0:
            self.custom_line_TP_grid = {}
            self.custom_line_paths = {}
            self.custom_diffTs, self.custom_diffPs = {}, {}
            
            for i_spec in range(len(self.line_species)):

                # Check if it is an Exomol hdf5 file that needs to be read:
                Chubb = False
                if self.mode == 'c-k':
                    path_opa = self.path+'/opacities/lines/corr_k/'+self.line_species[i_spec]
                    if glob.glob(path_opa+'/*.h5') != []:
                        Chubb = True

                # If not Exomol k-table made by Katy Chubb
                if not Chubb:

                    # Check and sort custom grid for species, if defined.
                    custom_grid_data = \
                      pyi.get_custom_PT_grid(self.path, \
                                             self.mode, \
                                             self.line_species[i_spec])

                    # If no custom grid was specified (no PTpaths.ls found):
                    # take nominal grid. This assumes that the files indeed
                    # are following the nominal grid and naming convention.
                    # Otherwise it will take the info provided in PTpaths.ls
                    # which was filled into custom_grid_data.
                    if custom_grid_data == None:
                        self.custom_line_TP_grid[self.line_species[i_spec]] = \
                          self.line_TP_grid
                        self.custom_line_paths[self.line_species[i_spec]] = None
                        self.custom_diffTs[self.line_species[i_spec]], \
                          self.custom_diffPs[self.line_species[i_spec]] = 13, 10
                        self.custom_grid[self.line_species[i_spec]] = False
                    else:
                        self.custom_line_TP_grid[self.line_species[i_spec]] = \
                          custom_grid_data[0]
                        self.custom_line_paths[self.line_species[i_spec]] = \
                          custom_grid_data[1]
                        self.custom_diffTs[self.line_species[i_spec]], \
                          self.custom_diffPs[self.line_species[i_spec]] = \
                          custom_grid_data[2], \
                          custom_grid_data[3]
                        self.custom_grid[self.line_species[i_spec]] = True
                else:
                    # If the user wants to make use of an Exomol k-table.
                    # In this case the custom grid is defined by reading
                    # the grid coordinates from the Exomol hdf5 file.

                    path_opa = self.path+'/opacities/lines/corr_k/'+self.line_species[i_spec]
                    file_path_hdf5 = glob.glob(path_opa+'/*.h5')[0]
                    f = h5py.File(file_path_hdf5,'r')
                    lent = len(f['t'][:])
                    lenp = len(f['p'][:])
                    retVal = np.zeros(lent*lenp*2).reshape(lent*lenp,2)
                    for i_t in range(lent):
                        for i_p in range(lenp):
                                                     # convert from bar to cgs
                            retVal[i_t*lenp+i_p, 1] = f['p'][i_p]*1e6
                            retVal[i_t*lenp+i_p, 0] = f['t'][i_t]
                    self.custom_line_TP_grid[self.line_species[i_spec]] = retVal
                    self.custom_diffTs[self.line_species[i_spec]], \
                      self.custom_diffPs[self.line_species[i_spec]] = \
                          lent, \
                          lenp
                    self.custom_grid[self.line_species[i_spec]] = True
                    f.close()
                    
        # Read actual opacities....
        # The nominal petitRADTRANS opacity grid "line_grid_kappas"
        # has the shape g_len,freq_len,len(line_species),len(line_TP_grid[:,0])
        # line_grid_kappas_custom_PT's entries have the shape
        # g_len,freq_len,len(self.custom_line_TP_grid[self.line_species[i_spec]])
        # From now on also the nominal grid opacities are read into
        # line_grid_kappas_custom_PT, because this makes things easier.
        self.line_grid_kappas_custom_PT = {}
                
        if len(self.line_species) > 0:

            tot_str = ''
            for sstring in self.line_species:
                tot_str = tot_str + sstring + ':'

            custom_file_names = ''

            for i_spec in range(len(self.line_species)):

                # Read in opacities in the petitRADTRANS format, either
                # in pRT P-T grid spacing or custom P-T grid spacing.

                # Check if it is an Exomol hdf5 file that needs to be read:
                Chubb = False
                if self.mode == 'c-k':
                    path_opa = self.path+'/opacities/lines/corr_k/'+self.line_species[i_spec]
                    if glob.glob(path_opa+'/*.h5') != []:
                        Chubb = True
                
                if not Chubb:

                    if not self.custom_grid[self.line_species[i_spec]]:
                        len_TP = len(self.line_TP_grid[:,0])
                    else:
                        len_TP = len(self.custom_line_TP_grid[ \
                                self.line_species[i_spec]][:,0])

                    custom_file_names = ''
                    if self.custom_grid[self.line_species[i_spec]]:
                        for i_TP in range(len_TP):
                            custom_file_names = custom_file_names + \
                                self.custom_line_paths[self.line_species[i_spec]][i_TP] \
                                + ':'

                    self.line_grid_kappas_custom_PT[self.line_species[i_spec]] = \
                      fi.read_in_molecular_opacities( \
                        self.path, \
                        self.line_species[i_spec]+':', \
                        self.freq_len_full, \
                        self.g_len, \
                        1, \
                        len_TP, \
                        self.mode, \
                        arr_min, \
                        arr_max, \
                        self.custom_grid[self.line_species[i_spec]], \
                        custom_file_names)

                # Read in the Exomol k-table by Katy Chubb if requested by the user
                else:
                    print('  Read line opacities of '+self.line_species[i_spec]+'...')

                    path_opa = self.path+'/opacities/lines/corr_k/'+self.line_species[i_spec]
                    file_path_hdf5 = glob.glob(path_opa+'/*.h5')[0]
                    f = h5py.File(file_path_hdf5,'r')

                    lenf = len(f['bin_centers'][:])
                    freqs_chubb = nc.c*f['bin_centers'][:][::-1]
                    lent = len(f['t'][:])
                    lenp = len(f['p'][:])

                    # Some swapaxes magic is required because the tables are sorted
                    # differently when coming from the Exomol website.
                    k_table = np.array(f['kcoeff'])
                    k_table = np.swapaxes(k_table, 0, 1)
                    k_table2 = k_table.reshape(lenp*lent, lenf, 16)
                    k_table2 = np.swapaxes(k_table2, 0, 2)
                    k_table2 = k_table2[:,::-1,:]

                    # Initialize an empty array that has the same spectral entries as
                    # pRT nominally. Only fill those values where the Exomol tables
                    # have entries.
                    retVal = np.zeros(self.g_len* self.freq_len_full* \
                       len(self.custom_line_TP_grid[self.line_species[i_spec]])).reshape( \
                                          self.g_len, self.freq_len_full, 1, \
                                          len(self.custom_line_TP_grid[self.line_species[i_spec]]))
                    index_fill = (self.freq_full <= freqs_chubb[0]) & \
                      (self.freq_full >= freqs_chubb[-1])
                    retVal[:, index_fill, 0, :] = k_table2

                    # Divide by mass to go from cross-sections to opacities, the latter
                    # is what pRT requires.
                    exomol_mass = float(f['mol_mass'][0])
                    self.line_grid_kappas_custom_PT[self.line_species[i_spec]] = retVal/exomol_mass/nc.amu
                    print(' Done.')
                    
                    f.close()

                # Cut the wavelength range of the just-read species to the wavelength range
                # requested by the user
                if self.mode == 'c-k':
                    self.line_grid_kappas_custom_PT[self.line_species[i_spec]] = \
                      np.array(self.line_grid_kappas_custom_PT[ \
                        self.line_species[i_spec]][:,index,0,:], \
                                 dtype='d',order='F')
                else:
                    self.line_grid_kappas_custom_PT[self.line_species[i_spec]] = \
                    np.array(self.line_grid_kappas_custom_PT[ \
                        self.line_species[i_spec]][:,:,0,:], \
                                 dtype='d',order='F')

            print()
            
        # Read in g grid for correlated-k
        if self.mode == 'c-k':
            buffer = np.genfromtxt(self.path+'/opa_input_files/g_comb_grid.dat')
            self.g_gauss, self.w_gauss = buffer[:,0], buffer[:,1]
            self.g_gauss,self.w_gauss = np.array(self.g_gauss,dtype='d', \
                order='F'),np.array(self.w_gauss, \
                dtype='d',order='F')
        elif self.mode == 'lbl':
            self.g_gauss, self.w_gauss = np.ones(1), np.ones(1)
            self.g_gauss,self.w_gauss = np.array(self.g_gauss,dtype='d', \
                order='F'),np.array(self.w_gauss, \
                dtype='d',order='F')

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
          np.array(rho_cloud_particles,dtype='d',order='F')
        self.cloud_specs_abs_opa = \
          np.array(cloud_specs_abs_opa,dtype='d',order='F')
        self.cloud_specs_scat_opa = \
          np.array(cloud_specs_scat_opa,dtype='d',order='F')
        self.cloud_aniso = np.array(cloud_aniso,dtype='d',order='F')
        self.cloud_lambdas = np.array(cloud_lambdas,dtype='d',order='F')
        self.cloud_rad_bins = np.array(cloud_rad_bins,dtype='d',order='F')
        self.cloud_radii = np.array(cloud_radii,dtype='d',order='F')
        
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
                                          len(P)),dtype='d',order='F')
        self.continuum_opa_scat = np.array(np.zeros(self.freq_len*len(P)) \
                                            .reshape(self.freq_len, \
                                            len(P)),dtype='d',order='F')
        self.contr_em = np.array(np.zeros(self.freq_len*len(P)) \
                                     .reshape(len(P),self.freq_len), \
                                     dtype='d',order='F')
        self.contr_tr = np.array(np.zeros(self.freq_len*len(P)) \
                                     .reshape(len(P),self.freq_len), \
                                     dtype='d',order='F')
        if len(self.line_species) > 0:
            self.line_struc_kappas = \
              np.array(np.zeros(self.g_len*self.freq_len*len(P)* \
                len(self.line_species)).reshape(self.g_len,self.freq_len, \
                len(self.line_species),len(P)),dtype='d',order='F')
            self.total_tau = \
              np.array(np.zeros_like(self.line_struc_kappas), \
                           dtype='d',order='F')
            self.line_abundances = \
              np.array(np.zeros(len(self.press)*len(self.line_species)) \
                           .reshape(len(self.press), \
                            len(self.line_species)),dtype='d',order='F')
        else: # If there are no specified line species then we need at
              # least an array to contain the continuum opas
              # I'll (mis)use the line_struc_kappas array for that
            self.line_struc_kappas = \
              np.array(np.zeros(self.g_len*self.freq_len*len(P)) \
                           .reshape(self.g_len,self.freq_len, \
                            1,len(P)),dtype='d',order='F')
            self.total_tau = \
              np.array(np.zeros_like(self.line_struc_kappas), \
                           dtype='d',order='F')
            self.line_abundances = \
              np.array(np.zeros(len(self.press)) \
                           .reshape(len(self.press),1), \
                           dtype='d',order='F')
        
        self.mmw = np.zeros_like(self.press)

        if len(self.cloud_species) > 0:
            self.cloud_mass_fracs = np.array(np.zeros(int(len(self.press))* \
                        int(len(self.cloud_species))) \
                            .reshape(int(len(self.press)), \
                            int(len(self.cloud_species))), \
                                                 dtype='d',order='F')
            self.r_g = np.array(np.zeros(int(len(self.press))* \
                        int(len(self.cloud_species))) \
                                    .reshape(int(len(self.press)), \
                                            int(len(self.cloud_species))), \
                                            dtype='d',order='F')

    def interpolate_species_opa(self, temp):
        # Interpolate line opacities to given temperature structure.
        self.temp = temp
        if len(self.line_species) > 0:
            for i_spec in range(len(self.line_species)):
                self.line_struc_kappas[:,:,i_spec,:] = fi.interpol_opa_ck(self.press,temp, \
                                    self.custom_line_TP_grid[self.line_species[i_spec]], \
                                    self.custom_grid[self.line_species[i_spec]], \
                                    self.custom_diffTs[self.line_species[i_spec]], \
                                    self.custom_diffPs[self.line_species[i_spec]], \
                                    self.line_grid_kappas_custom_PT[self.line_species[i_spec]])
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
            self.scat = True
            self.calc_cloud_opacity(abundances, mmw, gravity, \
                                        sigma_lnorm, fsed, Kzz, radius, \
                                        add_cloud_scat_as_abs)

        # Calculate rayleigh scattering opacities
        if len(self.rayleigh_species) != 0:
            self.scat = True
            self.add_rayleigh(abundances)
        # Add gray cloud deck
        if (self.Pcloud != None):
            self.continuum_opa[:,self.press>self.Pcloud*1e6] += 1e99
        # Add power law opacity
        if (self.kappa_zero != None):
            self.scat = True
            wlen_micron = nc.c/self.freq/1e-4
            scattering_add = self.kappa_zero * \
                (wlen_micron/0.35)**self.gamma_scat
            self.continuum_opa_scat += \
                np.repeat(scattering_add[None], \
                int(len(self.press)), axis = 0).transpose()

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

        # aniso = (1-g)
        cloud_abs, cloud_abs_plus_scat_aniso, aniso, cloud_abs_plus_scat_no_aniso = \
           fs.interp_integ_cloud_opas(cloud_abs_opa_TOT,cloud_scat_opa_TOT, \
            cloud_red_fac_aniso_TOT,self.cloud_lambdas,self.border_freqs)

        self.continuum_opa_scat += cloud_abs_plus_scat_no_aniso - cloud_abs

        if add_cloud_scat_as_abs != None:
            if add_cloud_scat_as_abs:
                self.continuum_opa += cloud_abs \
                    + 0.20*(cloud_abs_plus_scat_no_aniso - cloud_abs)
                    
                  #+ (cloud_abs_plus_scat_aniso - cloud_abs)

                  #+ 1.*(cloud_abs_plus_scat_no_aniso - cloud_abs)* \
                  #  (aniso)
                  
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
                      gray_opacity = None, Pcloud = None, \
                      kappa_zero = None, \
                      gamma_scat = None, \
                      add_cloud_scat_as_abs = None):
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
                Pcloud (Optional[float]):
                    Pressure, in bar, where opaque cloud deck is added to the
                    absorption opacity.
                kappa_zero (Optional[float]):
                    Scarttering opacity at 0.35 micron, in cgs units (cm^2/g).
                gamma_scat (Optional[float]):
                    Has to be given if kappa_zero is definded, this is the
                    wavelength powerlaw index of the parametrized scattering
                    opacity.
                add_cloud_scat_as_abs (Optional[bool]):
                    If ``True``, 20 % of the cloud scattering opacity will be
                    added to the absorption opacity, introduced to test for the
                    effect of neglecting scattering.
        '''

        self.Pcloud = Pcloud
        self.kappa_zero = kappa_zero
        self.gamma_scat = gamma_scat
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
                        kappa_zero = None, \
                        gamma_scat = None, \
                        contribution = False, haze_factor = None, \
                        gray_opacity = None, variable_gravity=True):
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
                    absorption opacity.
                kappa_zero (Optional[float]):
                    Scarttering opacity at 0.35 micron, in cgs units (cm^2/g).
                gamma_scat (Optional[float]):
                    Has to be given if kappa_zero is definded, this is the
                    wavelength powerlaw index of the parametrized scattering
                    opacity.
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
        self.kappa_zero = kappa_zero
        self.gamma_scat = gamma_scat
        self.mix_opa_tot(abunds,mmw,gravity,sigma_lnorm,fsed,Kzz,radius)
        self.calc_tr_rad(P0_bar,R_pl,gravity,mmw,contribution,variable_gravity)
        

    def calc_flux_transm(self,temp,abunds,gravity,mmw,P0_bar,R_pl,\
                             sigma_lnorm = None, \
                             fsed = None, Kzz = None, radius = None, \
                             Pcloud = None, \
                             kappa_zero = None, \
                             gamma_scat = None, \
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
                    absorption opacity.
                kappa_zero (Optional[float]):
                    Scarttering opacity at 0.35 micron, in cgs units (cm^2/g).
                gamma_scat (Optional[float]):
                    Has to be given if kappa_zero is definded, this is the
                    wavelength powerlaw index of the parametrized scattering
                    opacity.
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
        self.kappa_zero = kappa_zero
        self.gamma_scat = gamma_scat
        self.interpolate_species_opa(temp)
        self.mix_opa_tot(abunds,mmw,gravity,sigma_lnorm,fsed,Kzz,radius, \
                             add_cloud_scat_as_abs = add_cloud_scat_as_abs)
        self.calc_opt_depth(gravity)
        self.calc_RT(contribution)
        self.calc_tr_rad(P0_bar,R_pl,gravity,mmw,contribution,variable_gravity)

    def get_opa(self,temp):
        ''' Method to calculate and return the line opacities (assuming an abundance
        of 100 % for the inidividual species) of the Radtrans object. This method
        updates the line_struc_kappas attribute within the Radtrans class. For the
        low resolution (`c-k`) mode, the wavelength-mean within every frequency bin
        is returned.

            Args:
                temp:
                    the atmospheric temperature in K, at each atmospheric layer
                    (1-d numpy array, same length as pressure array).

            Returns:
                * wavelength in cm (1-d numpy array)
                * dictionary of opacities, keys are the names of the line_species
                  dictionary, entries are 2-d numpy arrays, with the shape
                  being (number of frequencies, number of atmospheric layers).
                  Units are cm^2/g, assuming an absorber abundance of 100 % for all
                  respective species.

        '''
        
        # Function to calc flux, called from outside
        self.interpolate_species_opa(temp)
       
        return_opas = {}

        resh_wgauss = self.w_gauss.reshape(len(self.w_gauss), 1, 1)
    
        for i_spec in range(len(self.line_species)):
            return_opas[self.line_species[i_spec]] = np.sum( \
                self.line_struc_kappas[:, :, i_spec, :] * \
                resh_wgauss, axis = 0)
        
        return nc.c/self.freq, return_opas

