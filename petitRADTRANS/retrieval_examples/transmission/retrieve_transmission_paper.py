import numpy as np
import sys
import emcee
import pickle as pickle
import time
from emcee.utils import MPIPool

from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc
import master_retrieval_model as rm
import rebin_give_width as rgw

sys.stdout.flush()

################################################################################
################################################################################
### Hyperparameter setup, where to save things to / read things from
################################################################################
################################################################################

retrieval_name = 'JWST_transmission_petitRADTRANSpaper'
absolute_path = '' # end with forward slash!
observation_files = {}
observation_files['NIRISS SOSS'] = 'NIRISS_SOSS.dat'
observation_files['NIRSpec G395M'] = 'NIRSpec_G395M.dat'
observation_files['MIRI LRS'] = 'MIRI_LRS.dat'

# For testing purposes only
plotting = False
if plotting:
    import pylab as plt

# Retrieval hyperparameters
stepsize = 1.75
n_walkers = 240
n_iter = 4200

cluster = False       # Submit to cluster
n_threads = 1         # Use mutliprocessing (local = 1)
write_threshold = 200 # number of iterations after which diagnostics are updated

# Wavelength range of observations, fixed parameters that will not be retrieved
WLEN = [0.8, 14.0]
LOG_G =  2.58
R_pl =   1.84*nc.r_jup_mean
R_star = 1.81*nc.r_sun

################################################################################
################################################################################
### READ IN OBSERVATION
################################################################################
################################################################################

# Read in data, convert all to cgs!
data_wlen = {}
data_flux_lambda = {}
data_flux_lambda_error = {}
data_wlen_bins = {}

for name in observation_files.keys():

    dat_obs = np.genfromtxt(observation_files[name])
    data_wlen[name] = dat_obs[:,0]*1e-4
    data_flux_lambda[name] = dat_obs[:,1]
    data_flux_lambda_error[name] = dat_obs[:,2]
    
    data_wlen_bins[name] = np.zeros_like(data_wlen[name])
    data_wlen_bins[name][:-1] = np.diff(data_wlen[name])
    data_wlen_bins[name][-1] = data_wlen_bins[name][-2]

################################################################################
################################################################################
### MODEL SETUP
################################################################################
################################################################################

### Create and setup radiative transfer object
# Create random P-T profile to create RT arrays of the Radtrans object.
temp_params = {}
temp_params['log_delta'] = -6.
temp_params['log_gamma'] = np.log10(0.4)
temp_params['t_int'] = 750.
temp_params['t_equ'] = 0.
temp_params['log_p_trans'] = -3.
temp_params['alpha'] = 0.
p, t = nc.make_press_temp(temp_params)

# Create the Ratrans object here
rt_object = Radtrans(line_species=['H2', 'CO_all_iso', 'H2O', \
                                  'CH4', 'NH3', 'CO2', 'H2S', \
                                  'Na', 'K'], \
                    rayleigh_species=['H2','He'], \
                    continuum_opacities = ['H2-H2','H2-He'], \
                    mode='c-k', \
                    wlen_bords_micron = WLEN)

# Create the RT arrays of appropriate lengths
rt_object.setup_opa_structure(p)

################################################################################
################################################################################
###   PRIOR SETUP (For all priors, neglect the constant that would be added 
###         because of the normalization of the normal distribution)
################################################################################
################################################################################

def b_range(x, b):
    if x > b:
        return -np.inf
    else:
        return 0.

def a_b_range(x, a, b):
    if x < a:
        return -np.inf
    elif x > b:
        return -np.inf
    else:
        return 0.

log_priors = {}
log_priors['log_delta']      = lambda x: -((x-(-5.5))/2.5)**2./2.                           
log_priors['log_gamma']      = lambda x: -((x-(-0.0))/2.)**2./2. 
log_priors['t_int']          = lambda x: a_b_range(x, 0., 1500.)
log_priors['t_equ']          = lambda x: a_b_range(x, 0., 4000.)
log_priors['log_p_trans']    = lambda x: -((x-(-3))/3.)**2./2.
log_priors['alpha']          = lambda x: -((x-0.25)/0.4)**2./2.
log_priors['log_g']          = lambda x: a_b_range(x, 2.0, 3.7) 
log_priors['log_P0']         = lambda x: a_b_range(x, -4, 2.)

# Priors for log mass fractions
log_priors['CO_all_iso']     = lambda x: a_b_range(x, -10., 0.)
log_priors['H2O']            = lambda x: a_b_range(x, -10., 0.)
log_priors['CH4']            = lambda x: a_b_range(x, -10., 0.)
log_priors['NH3']            = lambda x: a_b_range(x, -10., 0.)
log_priors['CO2']            = lambda x: a_b_range(x, -10., 0.)
log_priors['H2S']            = lambda x: a_b_range(x, -10., 0.)
log_priors['Na']             = lambda x: a_b_range(x, -10., 0.)
log_priors['K']              = lambda x: a_b_range(x, -10., 0.)

################################################################################
################################################################################
### DEFINE LOG PROBABILITY
################################################################################
################################################################################

# Declare diagnostics
function_calls = 0
computed_spectra = 0
NaN_spectra = 0
delta_wt = write_threshold

start_time = time.time()
file_object = open(absolute_path + 'diag_' + \
                       retrieval_name + '.dat', 'w').close()

def calc_log_prob(params):

    log_delta, log_gamma, t_int, t_equ, log_p_trans, alpha, \
      log_g, log_P0 = params[:-8]

    # Make dictionary for modified Guillot parameters
    temp_params = {}
    temp_params['log_delta'] = log_delta
    temp_params['log_gamma'] = log_gamma
    temp_params['t_int'] = t_int
    temp_params['t_equ'] = t_equ
    temp_params['log_p_trans'] = log_p_trans
    temp_params['alpha'] = alpha

    # Make dictionary for log 'metal' abundances
    ab_metals = {}
    ab_metals['CO_all_iso']     = params[-8:][0]
    ab_metals['H2O']            = params[-8:][1]
    ab_metals['CH4']            = params[-8:][2]
    ab_metals['NH3']            = params[-8:][3]
    ab_metals['CO2']            = params[-8:][4]
    ab_metals['H2S']            = params[-8:][5]
    ab_metals['Na']             = params[-8:][6]
    ab_metals['K']              = params[-8:][7]
    
    global function_calls
    global computed_spectra
    global NaN_spectra
    global write_threshold
    
    function_calls += 1
    
    # Prior calculation of all input parameters
    log_prior = 0.

    # Alpha should not be smaller than -1, this
    # would lead to negative temperatures!
    if alpha < -1:
        return -np.inf
    
    for key in temp_params.keys():
        log_prior += log_priors[key](temp_params[key])
        
    log_prior += log_priors['log_g'](log_g)
    log_prior += log_priors['log_P0'](log_P0)

    # Metal abundances: check that their
    # summed mass fraction is below 1.
    metal_sum = 0.
    for name in ab_metals.keys():
        log_prior += log_priors[name](ab_metals[name])
        metal_sum += 1e1**ab_metals[name]

    if metal_sum > 1.:
        log_prior += -np.inf

    # Return -inf if parameters fall outside prior distribution
    if (log_prior == -np.inf):
        return -np.inf
    
    # Calculate the log-likelihood
    log_likelihood = 0.

    # Calculate the forward model, this
    # returns the wavelengths in cm and the planet radius
    # in R_jup.
    wlen, flux_lambda = \
            rm.retrieval_model_plain(rt_object, temp_params, log_g, \
                                         log_P0, R_pl, ab_metals)
                                         
    # Just to make sure that a long chain does not die
    # Unexpectedly:
    # Return -inf if retrieval model returns NaN values
    if np.sum(np.isnan(flux_lambda)) > 0:
        print("NaN spectrum encountered")
        NaN_spectra += 1
        return -np.inf

    # Convert to observation for transmission case
    flux_sq = (flux_lambda*nc.r_jup_mean/R_star)**2 

    # Calculate log-likelihood
    for instrument in data_wlen.keys():

        # Rebin model to observation
        flux_rebinned = rgw.rebin_give_width(wlen, flux_sq, \
                        data_wlen[instrument], data_wlen_bins[instrument])

        if plotting:
            plt.errorbar(data_wlen[instrument], \
                             data_flux_lambda[instrument], \
                             data_flux_lambda_error[instrument], \
                             fmt = 'o', \
                             zorder = -20, \
                             color = 'red')

            plt.plot(data_wlen[instrument], \
                             flux_rebinned, \
                             's', \
                             zorder = -20, \
                             color = 'blue')

        # Calculate log-likelihood
        log_likelihood += \
               -np.sum(((flux_rebinned - data_flux_lambda[instrument])/ \
                    data_flux_lambda_error[instrument])**2.)/2.

    if plotting:
        plt.plot(wlen, flux_sq, color = 'black')
        plt.xscale('log')
        plt.show()
        
    computed_spectra += 1

    # Write diagnostics file
    if (function_calls >= write_threshold):
        
        write_threshold += delta_wt
        hours = (time.time() - start_time)/3600.0
        info_list = [function_calls, computed_spectra, NaN_spectra, \
                 log_prior + log_likelihood, hours] 

        file_object = open(absolute_path + 'diag_' + \
                               retrieval_name + '.dat', 'a')

        for i in np.arange(len(info_list)):
            if (i == len(info_list) - 1):
                file_object.write(str(info_list[i]).ljust(15) + "\n")
            else:
                file_object.write(str(info_list[i]).ljust(15) + " ")
        file_object.close()

    print(log_prior + log_likelihood)
    print("--> ", function_calls, " --> ", computed_spectra)
    return log_prior + log_likelihood

def lnprob(x):
    return calc_log_prob(x)

################################################################################
################################################################################
### MCMC
################################################################################
################################################################################

n_dim = len(log_priors)

# Set initial position vectors in parameter space
p0 = [np.array([np.random.normal(loc = -5.5, scale = 2.5, size=1)[0], \
                np.random.normal(loc = 0., scale = 2., size=1)[0], \
                0.+1500.*np.random.uniform(size=1)[0], \
                0.+4000.*np.random.uniform(size=1)[0], \
                np.random.normal(loc = -3., scale = 3., size=1)[0], \
                np.random.normal(loc = -0.25, scale = 0.4, size=1)[0], \
                LOG_G,
                -4.+6.*np.random.uniform(size=1)[0], \
                -10.+10.*np.random.uniform(size=1)[0], \
                -10.+10.*np.random.uniform(size=1)[0], \
                -10.+10.*np.random.uniform(size=1)[0], \
                -10.+10.*np.random.uniform(size=1)[0], \
                -10.+10.*np.random.uniform(size=1)[0], \
                -10.+10.*np.random.uniform(size=1)[0], \
                -10.+10.*np.random.uniform(size=1)[0], \
                -10.+10.*np.random.uniform(size=1)[0]] \
                ) for i in range(n_walkers)]

print('Run pre-burn-in!')
print()

# Initialize the MPI-based pool used for parallelization.
if cluster:
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob, \
                                        a = stepsize, pool = pool)
else:
    if n_threads > 1: 
        sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob, \
                                            a = stepsize, threads = n_threads)
    else:
        sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob, \
                                            a = stepsize)

# Run the actual MCMC
pre_burn_in_runs = int(np.min([399, n_iter/10])) + 3
pos, prob, state = sampler.run_mcmc(p0, pre_burn_in_runs)

# Get the best-fit position
highest_prob_index = np.unravel_index(sampler.lnprobability.argmax(), \
                                          sampler.lnprobability.shape)
best_position = sampler.chain[highest_prob_index]

f = open('best_position_pre_burn_in_' + retrieval_name + '.dat', 'w')
f.write(str(best_position))
f.close()

# Run actual chain
p0 = [np.array([best_position[0]+np.random.normal(size=1)[0]*0.8, \
                best_position[1]+np.random.normal(size=1)[0]*0.5, \
                best_position[2]+np.random.normal(size=1)[0]*70., \
                best_position[3]+np.random.normal(size=1)[0]*200., \
                best_position[4]+np.random.normal(size=1)[0]*0.5, \
                best_position[5]+np.random.normal(size=1)[0]*0.1, \
                LOG_G, \
                best_position[7]+np.random.normal(size=1)[0]*0.2, \
                best_position[8]+np.random.normal(size=1)[0]*0.3, \
                best_position[9]+np.random.normal(size=1)[0]*0.3, \
                best_position[10]+np.random.normal(size=1)[0]*0.3, \
                best_position[11]+np.random.normal(size=1)[0]*0.3, \
                best_position[12]+np.random.normal(size=1)[0]*0.3, \
                best_position[13]+np.random.normal(size=1)[0]*0.3, \
                best_position[14]+np.random.normal(size=1)[0]*0.3, \
                best_position[15]+np.random.normal(size=1)[0]*0.3] \
                   ) for i in range(n_walkers)]

print('Run chain!')
print()

if cluster:
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob, \
                                        a = stepsize, pool = pool)
else:
    if n_threads > 1: 
        sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob, \
                                            a = stepsize, threads = n_threads)
    else:
        sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob, \
                                            a = stepsize) 

pos, prob, state = sampler.run_mcmc(p0, n_iter) 

if cluster:
    pool.close()

################################################################################
################################################################################
### Save results
################################################################################
################################################################################

f = open('chain_pos_' + retrieval_name + '.pickle','wb')
pickle.dump(pos,f)
pickle.dump(prob,f)
pickle.dump(state,f)
samples = sampler.chain[:, :, :].reshape((-1, n_dim))
pickle.dump(samples,f)
f.close()

with open('chain_lnprob_' + retrieval_name + '.pickle', 'wb') as f:
    pickle.dump([sampler.lnprobability], \
                f, protocol=pickle.HIGHEST_PROTOCOL)

################################################################################
################################################################################
### Done.
################################################################################
################################################################################
