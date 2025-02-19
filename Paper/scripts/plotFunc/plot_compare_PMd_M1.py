#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 16:21:49 2024

@author: nilanjana

This script is used to compare the tau distribution between PMd and M1 either using combined data acros two monkeys
or single monkey (always included single tip sessions of Tomy while comparing PMd/M1 structures)

The main function calls 3 subfunctions - compare_PMd_M1_allLayers  , plotHistogram_PMd_M1 , compute_bootstrap
 
Inputs of the main function -

combine_monkey  = 0 or 1 
bootstrapData = 0 or 1

Before running this function, the AlighTAU_FR_Layerwise.py and then combineMonkeyData.py needs to be run since it
calls the data processed by the above two functions 


"""

#import modules
import os, pickle
import matplotlib.pyplot as plt
from scipy.stats import kruskal, bootstrap
import numpy as np
import seaborn as sns
import pandas as pd

#add subfunctions

def compare_PMd_M1_allLayers(inFolder,plot_combineData,probe_area): 
    ''' 
 Loads and combines data across all layers
    Parameters
    ----------
    inFolder : path to find the combined monkey data or single monkey data for PMd/M1
    plot_combineData : int - 0 or 1
    probe_area : list -[PMd,M1]

    Returns
    -------
    TAU_PMd : array containing tau values
    TAU_M1 : array containing tau values

    '''
   
    if plot_combineData ==1:        
        data = inFolder
        with open(data,'rb') as file:
            TAU = pickle.load(file)           
            
        TAU_PMd = (np.concatenate(TAU['PMd']['all_layers']))
        TAU_M1= (np.concatenate(TAU['M1']['all_layers']))
               
    else:
        
        for i, i_probe in enumerate(probe_area):
            data = f'{inFolder}/TAU_FR_Lat_{i_probe}.pkl'
                           
            with open(data,'rb') as file:
                TAU = pickle.load(file)
                t = np.concatenate(((np.concatenate(TAU['L2/3'])),(np.concatenate(TAU['L5'])),(np.concatenate(TAU['L6']))))           
                
                if i_probe == 'PMd':
                    TAU_PMd = (t)
                else:
                    TAU_M1 = (t)                    
                             
    return TAU_PMd , TAU_M1


def plotHistogram_PMd_M1(TAU, probe_area, outFolder,monkey_name):
    '''
    

    Parameters
    ----------
    TAU : list - [TAU_PMd, TAU_M1] returned by running subfunction above
    probe_area : list [PMd, M1] - DO NOT CHANGE THE ORDER
    outFolder : path to store the histogram figures
    monkey_name : list ['Mourad', 'Tomy']
    Returns
    -------
    None.

    '''
    
    # Create a new figure and subplots for histograms structure wise
    fig, ax = plt.subplots(figsize=(12,8))
    
    # Define colors for each area
    area_color = {'PMd': '#B1C086', 'M1': '#EBCCFF'} #change 
    data = pd.DataFrame({
        'tau' : np.concatenate(TAU),
        'area' : ['PMd'] * len(TAU[0]) + ['M1'] * len(TAU[1])
        }
        )
    
    # Create violin plot with both areas side by side
    sns.violinplot(x='area', y='tau', data=data, palette=area_color,hue = None, inner='box',ax=ax)
    # Calculate medians and add custom legend entries
    for i, area in enumerate(probe_area):
        # Calculate the median
        med = np.median(data[data['area'] == area]['tau'])
        
        # Add a custom legend entry for the median
        ax.plot([], [], color=area_color[area], linestyle='dashed', linewidth=2, 
                label=f'{area} Median: {med:.2f}')
        
    s, p_value = kruskal(TAU[0],TAU[1])
    if p_value < 0.05:
        
            # Define the positions for the bracket and star
        x1, x2 = 0, 1  # x-positions of the two violins
        y_max = max(data['tau'])  # Maximum y-value for positioning the bracket
        y_bracket = y_max * 1.1  # Position the bracket slightly above the violins
        y_star = y_bracket * 1.02  # Position the star slightly above the bracket
        
        # Draw the bracket
        ax.plot([x1, x1, x2, x2], [y_bracket, y_star, y_star, y_bracket], color='black', linewidth=1.5)
        
        # Add the star
        ax.text((x1 + x2) / 2, y_star, '*', fontsize=20, ha='center', va='bottom', color='black')

        
    #This part is commneted out. It can be used to plot the horizontal median line for each violin plot
    # Get x-tick positions
    # pos = {area: i for i, area in enumerate(probe_area)}
    # for area in probe_area:
    #     med = np.median(data[data['area'] == area]['tau'])
    #     x_pos = pos[area]  # Get the x position of the category 
    
    # # Add median lines
    # for area in probe_area:
    #     med = np.median(data[data['area'] == area]['tau'])
    #     x_pos = pos[area]  # Get the x position of the category
    #     ax.hlines(y=med, xmin=x_pos - 0.5, xmax=x_pos + 0.5, 
    #           color='k', linestyle='dashed', linewidth=4, label=f'{area} Median: {med:.2f}')
        
    
    ax.set_xlabel('area', fontsize=20)
    ax.set_ylabel('tau (ms)', fontsize=20)
    ax.tick_params(axis='both', labelsize=18)
    ax.legend(loc='center', fontsize=16)
    
    if monkey_name is not None:
        # Set a suptitle for the entire figure
        fig.suptitle(f'{monkey_name}',fontsize= 25)
    # Add legend for median lines in the upper right corner
   
           
    plt.tight_layout()
    plt.savefig(os.path.join(outFolder,f'{monkey_name}_PMd_M1_histogram'))
    
def compute_bootstrap(data,n_iteration,sample_k,outFolder,monkey_name):
    '''
    

    Parameters
    ----------
    data : list - [TAU_PMd, TAU_M1]
        
    sample_k : int - number of samples to pick in each bootstrap iterature
    outFolder : path to store the bootstrap figures
    monkey_name : list - [Mourad, Tomy]

    Returns
    -------
    None.

    '''
    ''
    probe_area = ['PMd', 'M1']
    mu_PMd = []
    mu_M1 = []
    
    for i in range(n_iteration):
        
        subset_PMd = 10**(np.random.choice(data[0],sample_k,replace=False)) #converting back from log 10 to actual values
        subset_M1 = 10**(np.random.choice(data[1], sample_k, replace = False))##converting back from log 10 to actual values
        
        subset_PMd = subset_PMd[subset_PMd>50]
        subset_M1 = subset_M1[subset_M1>50]
        
        mu_PMd.append(np.median(subset_PMd))
        mu_M1.append(np.median(subset_M1))
     
    TAU = [np.array(mu_PMd), np.array(mu_M1)]
    
  #----plot the distributions of the sampling medians of PMd and M1 and see if they overlap--------
    
    # Create a new figure and subplots for histograms structure wise
    fig_hist, axs_hist = plt.subplots(1, 1, figsize=(8, 5), sharex=True,sharey = True)
    area_color = ['darkblue','purple']
                
   # Compute 95% Confidence Intervals
    ci_lower_0, ci_upper_0 = np.percentile(TAU[0], [2.5, 97.5]) #PMd
    ci_lower_1, ci_upper_1 = np.percentile(TAU[1], [2.5, 97.5]) #M1
    
    ci_L = [ci_lower_0, ci_lower_1]
    ci_H = [ci_upper_0, ci_upper_1]

    # Create a histogram plot for the first row of subplots
    for j,i_probe in enumerate(probe_area): 
        counts,bin_edges  = np.histogram(TAU[j])
        proportions = counts/len(TAU[j])
        
        axs_hist.hist(bin_edges[:-1], bin_edges, weights=proportions, color=area_color[j], alpha = 0.5, edgecolor='none', label=f'95% CI: [{ci_L[j]:.2f}, {ci_H[j]:.2f}]')
        #sns.kdeplot(TAU[j], color=area_color[j], ax=axs_hist, linewidth=2)
        
        axs_hist.axvline(ci_L[j], color=area_color[j], linestyle='dashed', linewidth=2)
        axs_hist.axvline(ci_H[j], color=area_color[j], linestyle='dashed', linewidth=2)
        axs_hist.legend(loc='upper right', fontsize=20)
       
    axs_hist.tick_params(axis='y', labelsize=18)
    axs_hist.tick_params(axis='x', labelsize=18)
    axs_hist.set_ylim([0,0.5])
    axs_hist.set_xlim([120,240])
    
    axs_hist.set_xlabel('Median Values',fontsize= 20)
    axs_hist.set_ylabel('Frequency in each bin', fontsize = 20)
    fig_hist.suptitle(f'{monkey_name}: Bootstrapped Medians and 95% Confidence Intervals',fontsize= 25)
    #axs_hist.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig(os.path.join(outFolder,f'{monkey_name}_BOOTSTRAPPED_{n_iteration}_{sample_k}'))    

#add main function

if __name__ == '__main__':
    
    #server paths 
    #server - VPN or inside the server
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
        
    monkey = ['Mourad', 'Tomy']
    probe_area = ['PMd', 'M1']

    combine_monkey= 1
    bootstrapData  = 0
    n_iter = 10000
    k_sample = 100
    
    if combine_monkey == 1:   
            
        inFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/combinedMonkey/combined_ST.pkl'
        outFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/results/SUA/combinedMonkey'
        if not os.path.exists(outFolder):
            os.makedirs(outFolder)
            
        # # -------Subfunction to generate the PMd/M1 arrays to compare
       
        TAU_PMd, TAU_M1 = compare_PMd_M1_allLayers(inFolder,combine_monkey,probe_area)  
        
        #plot histogram - PMD/M1
        plotHistogram_PMd_M1([TAU_PMd,TAU_M1], probe_area, outFolder,'Combined')
        
        s, p_value = kruskal(TAU_PMd, TAU_M1) 
       
        if p_value <0.05:
            print (f'Significant difference of tau between PMd and M1 with p value of {p_value:.2f} and statistic {s:.2f} in combinedMonkey' )
        else:
            print ('No significant difference between tau of PMd and M1')
        print(f'Combined PMd : {len(TAU_PMd)} and M1: {len(TAU_M1)}')
         
    else:
        
        # -------First compare PMd/ M1 --------
        #run for loop for individual monkey
        for m, monkey_name in enumerate(monkey): 
           outFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/results/SUA/{monkey_name}'
           if not os.path.exists(outFolder):
                os.makedirs(outFolder)
            
           if monkey_name == 'Mourad':
               inFolder =  f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/{monkey_name}/combined_Layerwise'
               TAU_PMd, TAU_M1 = compare_PMd_M1_allLayers(inFolder,combine_monkey,probe_area)
               
           elif monkey_name == 'Tomy': #then i have to load the data combined with single tip neurons
               inFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/{monkey_name}/ST/combined_ST.pkl'               
               with open (inFolder, 'rb') as f:
                   t =  pickle.load(f)
            
               TAU_PMd = (np.concatenate(t['PMd']['all_layers']))
               TAU_M1 = (np.concatenate(t['M1']['all_layers']))
               
           
            #plot the histogram
           plotHistogram_PMd_M1([TAU_PMd,TAU_M1], probe_area, outFolder,monkey_name)            
           s, p_value = kruskal(TAU_PMd, TAU_M1) 
           
           if p_value <0.05:
                print (f'Significant difference of tau between PMd and M1 with p value of {p_value:.2f} and statistic {s:.2f} in {monkey_name}' )
           else:
                print (f'No significant difference between tau of PMd and M1 in {monkey_name}')
            
           print(f'{monkey_name} PMd : {len(TAU_PMd)} and M1: {len(TAU_M1)}')
               
            # #plot the bootstrap - create another subfunction
        
            # if bootstrapData ==1 :
            #     compute_bootstrap([TAU_PMd,TAU_M1],n_iter,k_sample,outFolder,monkey_name)
                
                
              
                
        
        
     
       
       
       
       
       
       
       
       
       
       
       
       
       
       
    
    