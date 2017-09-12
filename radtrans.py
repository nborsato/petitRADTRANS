import fort_input as fi
import fort_spec as fs
import numpy as np
import nat_cst as nc
import copy as cp
import read_bin as rb

class radtrans:
    """ Class carrying out spectral calcs for a given set of opacities """
    
    def __init__(self,line_species=[],rayleigh_species=[],H2H2CIA=False,H2HeCIA=False, \
                     wlen_bords_micron=[0.05,300.], mode='c-k'):

        # Line-by-line or corr-k
        self.mode = mode

        # Line opacity species to be considered
        self.line_species = line_species

        # Rayleigh scattering species to be considered
        self.rayleigh_species = rayleigh_species

        # Include CIA?
        self.H2H2CIA = H2H2CIA
        self.H2HeCIA = H2HeCIA

        # Get path to all input data (opacities, grids, etc.)
        f = open('path.txt')
        lines = f.readlines()
        self.path = lines[1][:-1]
        f.close()

        # Read in frequency grid
        if self.mode == 'c-k':
            self.freq_len, self.g_len = fi.get_freq_len(self.path)
            freq_len_full = cp.copy(self.freq_len)
            self.freq = fi.get_freq(self.path,self.freq_len)
        elif self.mode == 'lbl':
            x,y = rb.read_bin(self.path+'/opacities/lines/line_by_line/'+ \
                                'CO_all_iso/sigma_05_110.K_0.000001bar.dat')
            self.freq_len = int(len(x)+0.01)
            self.g_len = 1
            freq_len_full = cp.copy(self.freq_len)
            self.freq = nc.c/x
        
        index = (nc.c/self.freq > wlen_bords_micron[0]*1e-4) & \
          (nc.c/self.freq < wlen_bords_micron[1]*1e-4)

        self.freq = np.array(self.freq[index],dtype='d',order='Fortran')
        self.freq_len = len(self.freq)
        self.lambda_angstroem = np.array(nc.c/self.freq/1e-8,dtype='d',order='Fortran')
        self.flux = np.array(np.zeros(self.freq_len),dtype='d',order='Fortran')
        self.transm_rad = np.array(np.zeros(self.freq_len),dtype='d',order='Fortran')

        # Read in opacity grid
        buffer = np.genfromtxt(self.path+'/opa_input_files/opa_PT_grid.dat')
        self.line_TP_grid = np.zeros_like(buffer)
        self.line_TP_grid[:,0] = buffer[:,1]
        self.line_TP_grid[:,1] = buffer[:,0]
        # Convert from bars to cgs
        self.line_TP_grid[:,1] = 1e6*self.line_TP_grid[:,1]
        self.line_TP_grid = np.array(self.line_TP_grid.reshape(len(self.line_TP_grid[:,1]),2),dtype='d',order='Fortran')
        # Get opa grid
        tot_str = ''
        for sstring in line_species:
           tot_str = tot_str + sstring + ':'

        # line_grid_kappas has the shape g_len,freq_len,len(line_species),len(line_TP_grid[:,0])
        if len(self.line_species) > 0:
            self.line_grid_kappas = fi.get_opas_ck(self.path,tot_str,freq_len_full,self.g_len, \
                                                       len(self.line_species),len(self.line_TP_grid[:,0]),self.mode)
            self.line_grid_kappas = np.array(self.line_grid_kappas[:,index,:,:],dtype='d',order='Fortran')

        # Read in g grid
        if self.mode == 'c-k':
            buffer = np.genfromtxt(self.path+'/opa_input_files/g_comb_grid.dat')
            self.g_gauss, self.w_gauss = buffer[:,0], buffer[:,1]
            self.g_gauss,self.w_gauss = np.array(self.g_gauss,dtype='d',order='Fortran'),np.array(self.w_gauss,dtype='d',order='Fortran')
        elif self.mode == 'lbl':
            self.g_gauss, self.w_gauss = np.ones(1), np.ones(1)
            self.g_gauss,self.w_gauss = np.array(self.g_gauss,dtype='d',order='Fortran'),np.array(self.w_gauss,dtype='d',order='Fortran')
        
        # Read in RT mu grid
        buffer = np.genfromtxt(self.path+'/opa_input_files/mu_points.dat')
        self.mu, self.w_gauss_mu = buffer[:,0], buffer[:,1]

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

        # Define to not run into trouble later
        self.Pcloud = None
        self.haze_factor = None
        
    # Preparing structures
    def setup_opa_structure(self,P):
        ''' Setup opacity array for structure, as well as pressure array '''
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

    def interpolate_species_opa(self,temp):
        ''' Interpolate opacities to given temperature structure. '''
        self.temp = temp
        if len(self.line_species) > 0:
            self.line_struc_kappas = fi.interpol_opa_ck(self.press,temp,self.line_TP_grid,self.line_grid_kappas)
        else:
            self.line_struc_kappas = np.zeros_like(self.line_struc_kappas)
            
    def mix_opa_tot(self,abundances,mmw):
        ''' Combine total line opacities, according to mass fractions (abundances). '''
        self.scat = False
        self.mmw = mmw
        for i_spec in range(len(self.line_species)):
            self.line_abundances[:,i_spec] = abundances[self.line_species[i_spec]]
        self.continuum_opa = np.zeros_like(self.continuum_opa)
        self.continuum_opa_scat = np.zeros_like(self.continuum_opa_scat)
        if self.H2H2CIA:
            self.continuum_opa = self.continuum_opa + fi.cia_interpol(self.freq,self.temp, \
                self.cia_h2h2_lambda,self.cia_h2h2_temp,self.cia_h2h2_alpha_grid, \
                    self.press,self.mmw,abundances['H2'],2.)
        if self.H2HeCIA:
            self.continuum_opa = self.continuum_opa + fi.cia_interpol(self.freq,self.temp, \
                self.cia_h2he_lambda,self.cia_h2he_temp,self.cia_h2he_alpha_grid, \
                self.press,self.mmw,np.sqrt(abundances['H2']*abundances['He']),np.sqrt(8.))
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
        
    def calc_flux(self,temp,abunds,gravity,mmw,contribution=False):
        ''' Function to calc flux, called from outside '''
        self.interpolate_species_opa(temp)
        self.mix_opa_tot(abunds,mmw)
        self.calc_opt_depth(gravity)
        self.calc_RT(contribution)

    def calc_transm(self,temp,abunds,gravity,mmw,P0_bar,R_pl,Pcloud=None,contribution=False,haze_factor=None):
        ''' Function to calc transm. spectrum, called from outside '''
        self.Pcloud = Pcloud
        self.interpolate_species_opa(temp)
        self.haze_factor = haze_factor
        self.mix_opa_tot(abunds,mmw)
        self.calc_tr_rad(P0_bar,R_pl,gravity,mmw,contribution)
        

    def calc_flux_transm(self,temp,abunds,gravity,mmw,P0_bar,R_pl,Pcloud=None,contribution=False):
        ''' Function to calc flux, called from outside '''
        self.Pcloud = Pcloud
        self.interpolate_species_opa(temp)
        self.mix_opa_tot(abunds,mmw)
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
