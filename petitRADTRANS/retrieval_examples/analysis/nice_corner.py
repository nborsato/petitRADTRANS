import numpy as np
import seaborn as sns
import pylab as plt

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)
plt.rc('text', usetex=True)

def nice_corner(samples, \
                parameter_names, \
                output_file, \
                N_samples = None, \
                parameter_ranges = None, \
                parameter_plot_indices = None, \
                true_values = None, \
                max_val_ratio = None):

    if N_samples != None:
        S = N_samples
    else:
        S = len(samples)

    try:
        if parameter_plot_indices == None:
            parameter_plot_indices = np.linspace(0, len(parameter_names)-1, \
                                        len(parameter_names)-1).astype('int')
    except:
        pass
                                        
    if max_val_ratio == None:
        max_val_ratio = 5.

    data_list = []
    labels_list = []
    range_list = []
    
    for i in parameter_plot_indices:
        
        data_list.append(samples[len(samples)-S:,i])
        labels_list.append(parameter_names[i])

        try:
            if parameter_ranges[i] == None:
                range_mean = np.mean(samples[len(samples)-S:,i])
                range_std = np.std(samples[len(samples)-S:,i])
                range_take = (range_mean-4*range_std, range_mean+4*range_std)
                range_list.append(range_take)
            else:
                range_list.append(parameter_ranges[i])
        except:
            range_mean = np.mean(samples[len(samples)-S:,i])
            range_std = np.std(samples[len(samples)-S:,i])
            range_take = (range_mean-4*range_std, range_mean+4*range_std)
            range_list.append(range_take)
            
    try:
        truths_list = []
        for i in parameter_plot_indices:
            truths_list.append(true_values[i])
    except:
        pass

    dimensions = len(parameter_plot_indices)

    # Set up the matplotlib figure
    f, axes = plt.subplots(dimensions, dimensions, figsize=(13, 13), \
                           sharex=False, sharey=False)
    i_col = 0
    i_lin = 0

    for ax in axes.flat:

        #print(i_col, i_lin, dimensions, len(axes.flat))
        if i_col < i_lin:
            x, y = data_list[i_col], data_list[i_lin]
            gridsize = 50
            ax.hexbin(x, y, cmap='bone_r', gridsize=gridsize, \
                          vmax = int(max_val_ratio*S/gridsize**2.), \
                          rasterized=True)

            ax.set_xlim([range_list[i_col][0], range_list[i_col][1]])
            ax.set_ylim([range_list[i_lin][0], range_list[i_lin][1]])

            try:
                ax.axhline(truths_list[i_lin], color = 'red', \
                               linestyle = '--')
                ax.axvline(truths_list[i_col], color = 'red', \
                               linestyle = '--')
            except:
                pass
            
            if i_col > 0:
                ax.get_yaxis().set_visible(False)

        elif i_col == i_lin:
            
            med = np.median(data_list[i_col])
            up = str(np.round(np.percentile(data_list[i_col], 84) - med,2))
            do = str(np.round(med - np.percentile(data_list[i_col], 16),2))
            med = str(np.round(med,2))
            med = med.split('.')[0]+'.'+med.split('.')[1].ljust(2, '0')
            up = up.split('.')[0]+'.'+up.split('.')[1].ljust(2, '0')
            do = do.split('.')[0]+'.'+do.split('.')[1].ljust(2, '0')
            
            ax.set_title(med+r'$^{+'+up+'}_{-'+do+'}$', \
                             fontdict=font, fontsize = 10)
            sns.distplot(data_list[i_col], bins=20, kde=False, \
                             rug=False, ax=ax, color = 'gray')
                             
            ax.set_xlim([range_list[i_col][0], range_list[i_col][1]])
            ax.get_yaxis().set_visible(False)
            try:
                ax.axvline(truths_list[i_col], color = 'red', \
                                   linestyle = '--')
            except:
                pass
            
            ax.axvline(float(med)+float(up), color = 'black', \
                           linestyle = ':', linewidth=1.0)
            ax.axvline(float(med), color = 'black', \
                           linestyle = ':', linewidth=1.0)
            ax.axvline(float(med)-float(do), color = 'black', \
                           linestyle = ':', linewidth=1.0)

        else:
            ax.axis('off')

        labels = ax.yaxis.get_ticklabels(which='both')
        labels[-1].set_visible(False)

        i_col += 1
        if i_col%dimensions == 0:
            i_col = 0
            i_lin += 1

    for i_col in range(dimensions):
        plt.sca(axes[dimensions-1, i_col])
        range_use = np.linspace(range_list[i_col][0], \
                                    range_list[i_col][1], 5)[1:-1]
                                    
        labels_use = []
        for r in range_use:
            labels_use.append(str(np.round(r,1)))
        plt.xticks(range_use[::2], labels_use[::2])
        plt.xlabel(labels_list[i_col])

    for i_lin in range(dimensions):
        plt.sca(axes[i_lin,0])
        plt.ylabel(labels_list[i_lin])

    plt.subplots_adjust(wspace = 0, hspace = 0)

    plt.savefig(output_file)
    plt.show()
    plt.clf()
