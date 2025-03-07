""" Create a code that gets the already computed PSD and fits the FOOOF model to it. It
estimates the aperiodic part, including the knee, which will be used to calculate the tau."""

import os
import numpy as np
import re
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
from fooof import FOOOF
import pandas as pd
import matplotlib.ticker as ticker
import pickle

#subfunctions to use in the main code

def plot_laminar_psd(psd, area_label, session_label):
    """
    Plot the PSD of the laminar data. PSD is a xarray with all the info needed.
    :param psd: xarray with the PSD data
    :param area_label: str, the area label
    :param session_label: str, the session label
    ----
    :return: fig
    """
    # Specify the number of colors you want
    num_colors = len(psd.channels.values)

    # Get the colormap from the name
    colormap = plt.get_cmap('tab20c')

    # Generate a list of colors
    colors = [colormap(i / (num_colors - 1)) for i in range(num_colors)]
    fig, ax = plt.subplots(1, 1, figsize=(10, 12))

    for i_ch, channel_name in enumerate(psd.channels.values):
        ax.loglog(psd.freqs.values,
                  psd.mean(dim='trials')[i_ch].values,
                  label=channel_name, color=colors[i_ch]) #linear values but log log plot
    ax.legend(loc='best')
    ax.set_xlabel('Frequency (Hz)', weight='bold')
    ax.set_ylabel('Power', weight='bold',fontsize=18)
    ax.tick_params(axis= 'both', labelsize = 18)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    
    fig.suptitle(f'{session_label} - {area_label}', weight='bold')
    sns.despine(fig)
    return fig


def compute_and_plot_fooof(psd, fit_freq,area_label, session_label):
    """
    Compute the fooof model based on the single-channel psd. Compute the knee and tau. Plot
    the model fitting and the aperiodic part in two different figures, each subplot one channel.
    :param psd: xarray with the PSD data
    :param area_label: str, the area label
    :param session_label: str, the session label
    
    Changes made as of September 2, 2024
    
    first fit the psd of each channel in the given freq range with knee to check if teh channel has a knee or not
    If you not, then fit in fixed mode. 
    
    
    ----
    :return: fig, fig2, all_knees, all_exp, all_tau, all_e, all_r_sq (updated on Sept 2)
    
    - 2 figure maps , 
    all knee values(nan if there is no knee for a channel)
    all exponent values, 
    all tau estimates ; nan if no knee in a channel
    all error and R^2 fit values
                                                                      
    """
    # fit the FOOOF model, plot both the full spectrum fitting and just the the aperiodic part (w knee and tau)
    fig, ax = plt.subplots(8, 4, figsize=(20, 20))
    fig2, ax2 = plt.subplots(8, 4, figsize=(20, 20))
    lw1 = 1.0
   
    #initialise to add the variables for all channels later
    all_knees = np.zeros(psd.shape[1]) #psd has shape trials x chan x timepoints
    all_e = np.zeros(psd.shape[1]) 
    all_r_sq = np.zeros(psd.shape[1]) 
    all_exp = np.zeros(psd.shape[1]) 
    all_tau= np.zeros(psd.shape[1]) 


    for i_ch, ch_name in enumerate(psd.channels.values):
       
        #first fit with knee
        fm = FOOOF(aperiodic_mode='knee')  #removed peak_threshold = 0.1 criteria        
        fm.fit(psd.freqs.values,
               psd.mean(dim='trials')[i_ch].values, freq_range=fit_freq)
                
        # Extract fit metrics
        error = fm.error_
        r_squared = fm.r_squared_
        knee_value = fm.aperiodic_params_[1] #the params are offset, knee, exponent return array from the model
        exp = fm.aperiodic_params_[2]        
        knee_freq = pow(knee_value,(1/exp))
        
        if (knee_freq > fit_freq[0]) & (knee_freq< fit_freq[1]): #if knee actually exists then take the slope from the knee fit
        
            all_exp[i_ch]=exp
            tau = (1 / (knee_freq * 2 * np.pi)) * 1000
            all_tau[i_ch] = tau
            all_knees[i_ch] = knee_freq
            
            #keeping the error obtained from each channel fit
            all_e[i_ch]=error
            all_r_sq[i_ch]= r_squared
            
            #now plot the original data, periodic and aperiodic fit         
            fm.plot(plot_aperiodic=True, plt_log=True,
                    ax=ax.flatten()[i_ch],
                    add_legend=False, aperiodic_kwargs={'lw': lw1, 'zorder': 10},
                    peak_kwargs={'lw': lw1}, data_kwargs={'lw': lw1},
                    model_kwargs={'lw': lw1})
            
            ax.flatten()[i_ch].axvline(knee_freq,color= 'blue', linestyle = '--')           
            ax.flatten()[i_ch].set_title(f'{ch_name} - knee : {knee_freq:.2f} e: {exp:.2f}', fontsize=16)

            ax.flatten()[i_ch].set_ylabel('log(Power)')
            ax.flatten()[i_ch].set_xlabel('log(Freq)')
            
            #now plot in fig2 in loglog space but plotting actual power vs Frequency and the aperiodic part
            ax2.flatten()[i_ch].loglog(fm.freqs, 10**(fm._ap_fit), label=ch_name, color = 'blue')
            ax2.flatten()[i_ch].loglog(psd['freqs'].values , psd.mean(dim='trials')[i_ch].values, color='k')            
            ax2.flatten()[i_ch].axvline(knee_freq, color='blue', linestyle='--')
            
            ax2.flatten()[i_ch].set_title(f'{ch_name} - knee : {knee_freq:.2f} e: {exp:.2f}')
            ax2.flatten()[i_ch].set_xlabel('Freq',fontsize=16)
            ax2.flatten()[i_ch].set_ylabel('Power',fontsize=16)
            ax2.flatten()[i_ch].tick_params(axis='both', labelsize=18)
            ax2.flatten()[i_ch].xaxis.set_major_formatter(ticker.ScalarFormatter())
            ax2.flatten()[i_ch].yaxis.set_major_formatter(ticker.ScalarFormatter())
            
            
        else:
            #in case no knee exists just fit again in fixed mode and take the exponent
            fm = FOOOF(aperiodic_mode='fixed')
            fm.fit(psd.freqs.values,
                   psd.mean(dim='trials')[i_ch].values, freq_range=fit_freq)
            
            exp = fm.aperiodic_params_[1]
            all_exp[i_ch]= exp
            all_tau[i_ch]= np.nan            
            all_knees[i_ch]= np.nan
            
            #keeping the error obtained from updated each channel fit
            all_e[i_ch]=error
            all_r_sq[i_ch]= r_squared
            
             #now plot the original data, periodic and aperiodic fit         
            fm.plot(plot_aperiodic=True, plt_log=True,
                     ax=ax.flatten()[i_ch],
                     add_legend=False, aperiodic_kwargs={'lw': lw1, 'zorder': 10},
                     peak_kwargs={'lw': lw1}, data_kwargs={'lw': lw1},
                     model_kwargs={'lw': lw1})
            
            ax.flatten()[i_ch].set_title(f'{ch_name} - knee : {np.nan} e: {exp:.2f}')
            ax.flatten()[i_ch].set_ylabel('log(Power)')
            ax.flatten()[i_ch].set_xlabel('log(Freq)')
            
            #now plot in fig2 in loglog space but plotting actual power vs Frequency and the aperiodic part
            
            ax2.flatten()[i_ch].loglog(fm.freqs, 10**(fm._ap_fit), label=ch_name)
            ax2.flatten()[i_ch].loglog(psd['freqs'].values , psd.mean(dim='trials')[i_ch].values, color='k')
            
            ax2.flatten()[i_ch].set_title(f'{ch_name} - knee : {np.nan} e: {exp:.2f}',fontsize=16)
            ax2.flatten()[i_ch].set_ylabel('Power',fontsize=16)
            ax2.flatten()[i_ch].set_xlabel('Freq',fontsize=16)
            ax2.flatten()[i_ch].tick_params(axis='both', labelsize=18)
            ax2.flatten()[i_ch].xaxis.set_major_formatter(ticker.ScalarFormatter())
            ax2.flatten()[i_ch].yaxis.set_major_formatter(ticker.ScalarFormatter())
    # Remove empty subplots from index 20 to the end
    for i in range(i_ch+1, len(ax.flatten())):
        fig.delaxes(ax.flatten()[i])  # Delete the unused subplots
    
    fig.suptitle(f'FOOOF fit with knee_loglog - {session_label} - {area_label}',
                  weight='bold', y=0.99)
    
    # Remove empty subplots from index 20 to the end
    for i in range(i_ch+1, len(ax2.flatten())):
        fig2.delaxes(ax2.flatten()[i])  # 
    
    fig2.suptitle(f'Aperiodic fitting - {session_label} - {area_label}',
                  weight='bold', y=0.99)
    
    fig.tight_layout()
    fig2.tight_layout()
    return fig, fig2, all_knees, all_exp, all_tau, all_e, all_r_sq

# Main function    

if __name__ == '__main__':
    
    monkey = 'Tomy'
    psd_freq = [0,100] #frequencies to compute the psd
    epoch = ['preSC1', 'preGo', 'full_trial']
    fit_freq =[2,80]
    
    # list all the sessions that we want to preprocess
    if monkey == 'Mourad': 
       LAMINAR_SESS = ['Mo180330001','Mo180405001','Mo180405004','Mo180411001','Mo180412002',
                    'Mo180418002','Mo180419003','Mo180426004','Mo180503002', 'Mo180523002','Mo180524003', 
                    'Mo180525003','Mo180531002','Mo180614002','Mo180614006',
                    'Mo180615002','Mo180615005', 'Mo180619002','Mo180620004','Mo180622002',
                    'Mo180626003', 'Mo180627003','Mo180629005', 'Mo180703003','Mo180704003', 'Mo180705002'
                      'Mo180706002', 'Mo180710002','Mo180711004','Mo180712006']
    elif monkey == 'Tomy':
    
        LAMINAR_SESS =['t140924003','t140925001','t140926002','t140930001','t141001001','t141008001','t141010003','t150122001',
                    't150123001','t150204001','t150205004','t150303002','t150319003','t150320002','t150327002',
                    't150327003','t150415002','t150416002','t150423002','t150430002','t150520003','t150716001','t150702002']
        #LAMINAR_SESS = ['t150702002']
    for e in epoch:
        for session in LAMINAR_SESS:
            # check where are we running the code
            current_path = os.getcwd()
    
            if current_path.startswith('/Users'):
                server = '/Volumes/work'  # local w VPN
            elif current_path.startswith('/home/'):
                server = '/envau/work'  # niolon
            elif current_path.startswith('/envau'):
                server = '/envau/work'  # niolon
    
            # set the path
            path = f'{server}/comco/nandi.n/IntrinsicTimescales/Paper/data/LFP/{monkey}/PSD/Epo_{e}/{session}_{psd_freq}'
            path_plots = f'{server}/comco/nandi.n/IntrinsicTimescales/Paper/results/LFP/{monkey}/PSD_FOOOF/Epo_{e}/{session}_{psd_freq}_{fit_freq}'
            path_fooof = f'{server}/comco/nandi.n/IntrinsicTimescales/Paper/data/LFP/{monkey}/FOOOF/Epo_{e}/{session}_{fit_freq}'
           
            if not os.path.exists(path_plots):
                os.makedirs(path_plots)
                
            if not os.path.exists(path_fooof):
                os.makedirs(path_fooof)

            filenames = []
            for i in os.listdir(path):
                if (os.path.isfile(os.path.join(path,i)) and f'{session}' in i and '.nc' in i and 
                'PSD' in i and 'bipolar' in i):
                    filenames.append(i)
            # iterate in the filenames (i.e. sites)
            data_list = []
            for filename_site in filenames:
                # load the data array with the PSD of each site
                PSD_single_trial = xr.load_dataarray(os.path.join(path, filename_site))
    
                # get the area
                area = re.split('-', PSD_single_trial.channels.values[0])[0]
    
                # plot the data itself, average across trials, each channel one line.
                fig1 = plot_laminar_psd(PSD_single_trial, area_label=area, session_label=session)
                fig1.savefig(os.path.join(path_plots, f'{session}_{area}_PSD.png'), dpi=300)
    
                # # compute and plot the FOOOF model
                fig_fooof, fig_aperiodic,all_knees, all_exp, all_tau, all_e,\
                    all_r_sq = compute_and_plot_fooof(PSD_single_trial, fit_freq,area_label=area,
                                                            session_label=session)
                              
                fig_fooof.savefig(os.path.join(path_plots,f'{session}_{area}_FOOOF_fit.png'), dpi=300)   
                fig_aperiodic.savefig(os.path.join(path_plots,f'{session}_{area}_FOOOF_aperiodic_fit.png'), dpi=300)
                plt.close('all') 
                    
            # store the average psd with all the info about the laminar knee, tau, exponent , error fits
                PSD_avg = PSD_single_trial.mean(dim='trials')
                PSD_avg = PSD_avg.assign_coords(laminar_knees=('channels', all_knees), laminar_exp=('channels', all_exp),
                                                laminar_tau_ms=('channels', all_tau), fit_err= ('channels',all_e),
                                                fit_R2=('channels',all_r_sq))
                # save the PSD_avg
                PSD_avg.to_netcdf(os.path.join(path_fooof, f'{session}_{area}_PSD_avg_fooof.nc'))
                
                
            #store the layer info for each channel - that is the exponent, fitr error and the tau to use later          
                layers_data = PSD_single_trial.corrected_layers.values
                
                for l, layer in enumerate(layers_data):
                    data_list.append({
                    'session': session,
                    'area': area,
                    'layer': layer,
                    'exp' : all_exp[l],
                    'knee': all_knees[l],
                    'tau': all_tau[l],
                    'err' : all_e[l],
                    'R^2' : all_r_sq[l]
                    
                    })
                
            df = pd.DataFrame(data_list)
            # # Save the dataframe as a pickle file
            df.to_pickle(f'{server}/comco/nandi.n/IntrinsicTimescales/Paper/data/LFP/{monkey}/FOOOF/Epo_{e}/{session}_{fit_freq}/layerwise_fooof_results.pkl')
            print (f'Completed session {session}')
            
            
            
    

            
      
