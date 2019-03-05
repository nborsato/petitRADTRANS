import numpy as np
import pickle

from petitRADTRANS import nat_cst as nc

def return_PT_envelopes(samples, \
                        envelope_file, \
                        N_samples = None, \
                        read = False, \
                        true_values = None):

    if N_samples == None:
        N_samples = len(samples)

    try:
        if not true_values.any() == None:
            temp_params_input = {}
            temp_params_input['log_delta'] = true_values[0]
            temp_params_input['log_gamma'] = true_values[1]
            temp_params_input['t_int'] = true_values[2]
            temp_params_input['t_equ'] = true_values[3]
            temp_params_input['log_p_trans'] = true_values[4]
            temp_params_input['alpha'] = true_values[5]

            # Compute inputinal input P-T profile 
            p, temp_input = nc.make_press_temp(temp_params_input)
    except:
        temp_params = {}
        temp_params['log_delta'] = -6.
        temp_params['log_gamma'] = np.log10(0.4)
        temp_params['t_int'] = 750.
        temp_params['t_equ'] = 0.
        temp_params['log_p_trans'] = -3.
        temp_params['alpha'] = 0.
        p, buffer = nc.make_press_temp(temp_params)

    if read:

        f = open(envelope_file, 'rb')
        temp_left_1 = pickle.load(f)
        temp_left_2 = pickle.load(f)
        temp_left_3 = pickle.load(f)
        temp_left_4 = pickle.load(f)
        temp_right_1 = pickle.load(f)
        temp_right_2 = pickle.load(f)
        temp_right_3 = pickle.load(f)
        temp_right_4 = pickle.load(f)
        f.close()

        #print(temp_left_1)

    else:

        temps = np.zeros(N_samples *len(p)).reshape(len(p), N_samples)

        i_str = 0
        i_start = 0
        for params in samples[int(len(samples))-N_samples:]:
            if i_start%100 == 0:
                print(str(i_start+1)+'\r',)
            i_start += 1
            temp_params = {}
            temp_params['log_delta'] = params[0]
            temp_params['log_gamma'] = params[1]
            temp_params['t_int'] = params[2]
            temp_params['t_equ'] = params[3]
            temp_params['log_p_trans'] = params[4]
            temp_params['alpha'] = params[5]
            p, temp = nc.make_press_temp(temp_params)

            if (np.max(temp) < 15000) and (np.min(temp) > 0.):
                temps[:,i_str] = temp
                i_str = i_str+1

                
        temp_left_1 = np.zeros_like(p)
        temp_left_2 = np.zeros_like(p)
        temp_left_3 = np.zeros_like(p)
        temp_left_4 = np.zeros_like(p)

        temp_right_1 = np.zeros_like(p)
        temp_right_2 = np.zeros_like(p)
        temp_right_3 = np.zeros_like(p)
        temp_right_4 = np.zeros_like(p)

        for i in range(len(p)):
            sort_temp = np.sort(temps[i,:i_str])
    
            temp_left_1[i] = sort_temp[int(0.05*i_str)]
            temp_left_2[i] = sort_temp[int(0.15*i_str)]
            temp_left_3[i] = sort_temp[int(0.25*i_str)]
            temp_left_4[i] = sort_temp[int(0.35*i_str)]
            
            temp_right_1[i] = sort_temp[int(0.95*i_str)]
            temp_right_2[i] = sort_temp[int(0.85*i_str)]
            temp_right_3[i] = sort_temp[int(0.75*i_str)]
            temp_right_4[i] = sort_temp[int(0.65*i_str)]
            
        f = open(envelope_file,'wb')
        pickle.dump(temp_left_1, f)
        pickle.dump(temp_left_2, f)
        pickle.dump(temp_left_3, f)
        pickle.dump(temp_left_4, f)
        pickle.dump(temp_right_1, f)
        pickle.dump(temp_right_2, f)
        pickle.dump(temp_right_3, f)
        pickle.dump(temp_right_4, f)
        f.close()

    #print(temp_left_1)

    try:
        if not true_values.any() == None:
            return [[p, temp_input], \
                     [5, 15, 25, 35], \
                     [temp_left_1, temp_left_2, temp_left_3, temp_left_4], \
                     [temp_right_1, temp_right_2, temp_right_3, temp_right_4]]
    except:
        return [[p, None], \
                     [5, 15, 25, 35], \
                     [temp_left_1, temp_left_2, temp_left_3, temp_left_4], \
                     [temp_right_1, temp_right_2, temp_right_3, temp_right_4]]
