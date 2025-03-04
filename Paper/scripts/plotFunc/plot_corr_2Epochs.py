#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 19:09:40 2024
This script is used to check the correlation of TAU values in the 2 epochs - currently using preSC1 and preGo epochs

Modified on Feb 21, 2025 to separately check the PMd / M1 sessions -- removing the criteria of disckarding with skewness to keep the same std threshold 
criteria as the whole session fits

@author: nilanjana
"""

#import files from pre_SC1 and pre_Go directories , then scatter them and check if they have a correlation
#%% sufunctions to use

def loadPickleFile (filename):  
    with open(filename,'rb') as file:
        fit_data =  pickle.load(file)
        autoc =  pickle.load(file)
    return fit_data, autoc

def loadTAU (fit_val) :
            
        if fit_val['bestfit'][0]=='mono':         
            mono_fit = fit_val[fit_val['fit_type'] == 'mono'] #find the mono fit        
            #print('Inside mono loop')
            t = mono_fit['TAU'].values
            se = np.round(mono_fit['TAU_sd'].values,2)
        else:
            t = np.nan
            se = np.nan

        return t,se
    

def plot_allSess(server,session_names,probe_num,area,ax,include_all_neurons,combine_m):
    
    '''
    Input Params
    server  - the current working directory
    session_names - list of all sessions for a given monkey and probe
    probe_num - list of probe numbers (1 or 2) corresponding to PMd/M1 in a given session folder
    area = str - PMd or M1
    ax - plot axis
    include_all_neurons - 0 or 1 - for 0 , all neurons are plotted irrespective of std error threshold- they are just masked in red
                            1 means neurons are plotted after removing based on std 
                            
    Output - figure map with corr plot in two monkeys in PMd/M1
    '''

    #defining plot colors    
    if area == 'PMd':        
        color_area = '#B1C086'
    else:
        color_area = '#EBCCFF'
    #definig empty lists to concatenate later
    tau_preGo = [] 
    se_preGo = []
    tau_preSC1 =[] 
    se_preSC1 = []
    for i, i_session in enumerate(session_names):        
        if i_session.startswith('Mo'):
             m_name = 'Mourad'
        else:
             m_name = 'Tomy'
       
        inFolder_TAU_preGo = f"{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/{m_name}/{i_session}/TAU_preGo_firstBinExc"
        inFolder_TAU_preSC1 = f"{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/{m_name}/{i_session}/TAU_preSC1_firstBinExc"

        if os.path.exists(inFolder_TAU_preGo):
            #search for filenames in the preGo session folder where the probe number should belong to the 'area' (PMd/M1)
            filenames_preGo= []
            
            for filename in os.listdir(inFolder_TAU_preGo):
                if (filename.endswith('.pkl')) and (f'probe{probe_num.iloc[i]}' in filename):                
                   filenames_preGo.append(filename)
            
            for j, filename in enumerate(filenames_preGo):
                path_preGo = os.path.join(inFolder_TAU_preGo,filename)
                path_preSC1 = os.path.join(inFolder_TAU_preSC1,filename)
                
                fit_preGo, autoc_preGo = loadPickleFile(path_preGo)
                fit_preSC1, autoc_preSC1 =  loadPickleFile(path_preSC1)
     #concatenated all tau, sd, from that session folder for both epochs for a single session     
                tau1,se1= loadTAU(fit_preGo)
                tau2,se2 = loadTAU(fit_preSC1)

                tau_preGo.append(tau1)
                tau_preSC1.append(tau2)
                
                se_preGo.append(se1)
                se_preSC1.append(se2)

        else:
            print(f'Session {i_session} does not exist for 2 epochs')
    if combine_m == 1: 
        title_m = 'combined'
    else:
        title_m = m_name
            
 #Now concatenate across all sessions      
    t_preGo = np.hstack(tau_preGo) #works like np.concatenate but also takes into account the nan values. with np.concatenate the array will give nan
    t_preSC1 =  np.hstack(tau_preSC1)

    std_preGo = np.hstack(se_preGo)
    std_preSC1 = np.hstack(se_preSC1)
 
# to remove the nan values if its present for a datapoint in either of the arrays 
   
    valid_indices = ~np.isnan(t_preGo) & ~np.isnan(t_preSC1)
    x_valid = t_preGo[valid_indices]  #x_valid is preGo 
    y_valid = t_preSC1[valid_indices] #y_valid is preSC1
    
    std_err_x = std_preGo[valid_indices]
    std_err_y = std_preSC1[valid_indices]
    
       # Set a threshold for standard error
    fiterr_threshold = 150 #same threshold for TAU fits for the whole session
    #sk_threshold = 5
   
    if not np.isnan(x_valid).all():       
        if include_all_neurons ==1 :
        
            # Create a combined mask based on the condition
            
            mask =  (std_err_x > fiterr_threshold) | (std_err_y > fiterr_threshold)
            #plot the scatter plot to show the correlation             
            ax.scatter(x_valid,y_valid,color = np.where(mask,'red',f'{color_area}'), label=f'{len(x_valid)}')    
            
        else:        
            goodNeurons_indices = (std_err_x < fiterr_threshold) & (std_err_y < fiterr_threshold)
            
            x_valid = x_valid[goodNeurons_indices]
            y_valid = y_valid[goodNeurons_indices]
            ax.scatter(x_valid,y_valid,color = color_area, label = f'{len(x_valid)}')
            
        #correlation of TAU in the two epochs    
        if len(x_valid)>2:
            c, p_value = pearsonr(x_valid, y_valid)               
            alpha = 0.05
            is_significant = p_value < alpha
              
            print("Correlation Coefficient:", c)
            print("P-value:", p_value)
            print("Is the correlation significant?", is_significant)
            
            s,pval = kruskal(x_valid,y_valid)
            significance = pval < alpha
            
            print(f"Comparing Median  between preSC1 and PreGo in monkey {title_m} in {area}")
            print("P-value:", pval)
            print("Is the difference significant?", significance)
            
            ax.set_xlabel(f'preGo -- Median_Tau = {np.median(x_valid):.2f}', fontsize = 14)
            ax.set_ylabel(f'preSC1 -- Median_Tau = {np.median(y_valid):.2f}', fontsize = 14)
            ax.set_title(f'{title_m}-{area} -- corr:{c:.2f} p:{p_value:.2f}', fontsize = 14)
                    
            # Fit a linear regression through origin
            slope_origin = np.mean(y_valid) / np.mean(x_valid) 
            y_range_origin = slope_origin * x_valid
            ax.plot(x_valid, y_range_origin, color='green', linestyle='--', label = '')
            ax.tick_params(axis = 'both', labelsize = 18)
            
            # Add correlation coefficient and p-value as text on the plot
            text_str = f'Corr: {c:.2f}\np: {p_value:.2f}'
            ax.text(0.05, 0.95, text_str, transform=ax.transAxes, fontsize=14,
                    verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))
                                    
        ax.legend(loc = 'upper right', fontsize = 14)          
       
        # In case you want to find the regression line that fits the data instead of passing only through origin
        
        #     slope, intercept, _, _, _ = linregress(x_valid, y_valid)
        #     slope_origin = np.mean(y_valid) / np.mean(x_valid)  # Slope for the line passing through the origin
            
        #     # Plot the regression line
        #     x_range = np.linspace(0, max(x_valid), 100)
        #     y_range = slope * x_range + intercept
        #     y_range_origin = slope_origin * x_range
          
        # Display the legend
        #plt.tight_layout()
        #sns.despine()
        
    else:
        print(f'No neurons left for comparison between 2 epochs in this {area} in monkey {m_name}')
        
        
       
            
#%%
if __name__=='__main__' :   
    #modules import
    import os, pickle
    import matplotlib.pyplot as plt 
    import numpy as np
    from scipy.stats import pearsonr, linregress, kruskal, skew
    import seaborn as sns
    import pandas as pd
    
    
    #path initialise
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
        
    monkey_names = ['Mourad', 'Tomy']
    probe_area = ['PMd', 'M1']    
    combine_monkey = [1] #combine monkey data or not

    
    for combine_m in combine_monkey:
        fig, axs = plt.subplots(2,2, figsize=(12,8), sharex=True, sharey=True)
        include_all_neurons = 0 # include all neurons based on error threshold (1) or not (0) for the correlation plot
        
        for p, area in enumerate(probe_area): #Pmd or M1
            df_sess_c =  []
            df_probeNum_c = []
            for m, monkey in enumerate(monkey_names):                
                #read the csv file for AP and area                      
                df = pd.read_excel(f'{server}/work/comco/nandi.n/IntrinsicTimescales/docs/{monkey}LaminarInfo.xlsx')                                 
                df_sess= df[df['Probe_Area'] == area]['Session'] #subset of sessions with specific probe area
                df_probeNum= df[df['Probe_Area'] == area]['Probe_Num'] #probe number for the corresponding area - PMd/M1
                ax = axs[m,p]            
                df_probeNum = df_probeNum.astype('int')
                if combine_m ==0:                                           
                    print(f'Plotting for single monkey {monkey} and probe area {area}')
                    # pass the df and plot                
                    plot_allSess(server,df_sess,df_probeNum,area,ax,include_all_neurons,combine_m)                                        
                    plt.tight_layout()
                   
                else:
                    #combine the df_sess, df_probeNum for later
                    df_sess_c.append(df_sess)
                    df_probeNum_c.append(df_probeNum)
                    
            if combine_m ==1:                              
                print(f'Plotting for combined monkey data in each probe area {area}')
                plot_allSess(server,pd.concat(df_sess_c),pd.concat(df_probeNum_c),area,ax,include_all_neurons,combine_m)
                # Remove the empty subplots
                plt.tight_layout()
                
        outFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/results/SUA/combinedMonkey' 
        if not os.path.exists(outFolder):
             os.makedirs(outFolder)           
        plt.savefig(os.path.join(outFolder,f'Corr_TAU_preSC1_preGo_{combine_m}'))
               
                   
            

    
    
