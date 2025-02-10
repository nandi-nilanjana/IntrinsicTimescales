#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 16:50:09 2024

@author: nilanjana
"""

#combine TomyAlmapSpikes2023 and the rest 


#load the csv data for Tomy 
#load the other M1 PMd data for Tomy and just combine all layers and give that as an output

#saving it as TAU_FR_LAT_LayerWise_M1 or TAU_FR_LAT_LayerWise_PMd to integrate inside the plot_compare_PMd_M1_TAU function
#also save in separate extension folder in Tomy


#best is load the layerwise data from all sessins of and just combine the M1 and PMd all layers extracted from TomyAlmap2023.tx

#%%
#import modules

import os
import pickle
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import kruskal, bootstrap

#%%


#========== Subfunction to concatenate the tau, lat and the firing rate of neurons across layers for a particular session and monkey



def TAU_allLayers(data_path_TAU):
   
    tau_array = []
   
    
    # Search for pickle files in the folder for that session
    
    file_names_TAU = [filename for filename in os.listdir(data_path_TAU) if filename.endswith(".pkl")]

    
    for i, file_name in enumerate(file_names_TAU):
        
        file_path_TAU = os.path.join(data_path_TAU, file_name) 
      
        
        
        with open(file_path_TAU,'rb') as file: #load TAU Data
            fit_val = pickle.load(file)
                
           
            
        if fit_val['bestfit'][0] == 'mono':
            
            mono_fit = fit_val[fit_val['fit_type'] == 'mono'] #find the mono fit
            
            if mono_fit['TAU_sd'].values <150:
            
                t_val = mono_fit['TAU'].iloc[0]
                
                tau_array.append(t_val)
                       
    return np.array(tau_array)

# Define a function to check if a value is empty
def is_empty(value):
    return value is None or (isinstance(value, (str, list)) and len(value) == 0)


#%%
def plotHistogram_PMd_M1(TAU, probe_area, outFolder,monkey_name):
    
    
    # Create a new figure and subplots for histograms structure wise
    fig_hist, axs_hist = plt.subplots(1, 2, figsize=(8, 3), sharex=True)
    area_color = ['darkblue','purple']
    # Define the desired x-axis ticks
    x_ticks = [1, 2, 3]  # Corresponding to 10, 100, 1000 in non-log scale
    
    
    # Create a histogram plot for the first row of subplots
    for j,i_probe in enumerate(probe_area): 
        c2 = len(TAU[j])
        
        med = np.median(TAU[j])
       
        
        axs_hist[j].hist(TAU[j], bins=20, color=area_color[j], alpha = 0.5, edgecolor='black')
        axs_hist[j].legend([f"{c2}"], loc='upper right',fontsize=20)
      #  Set y-axis limits to ensure it ranges from 0 to 1
        axs_hist[j].set_ylim(0, 75)
        
    
    
        axs_hist[j].set_xticks(x_ticks, [f'{10 ** tick:.0f}' for tick in x_ticks], fontsize= 16)
        axs_hist[j].tick_params(axis='y', labelsize=12)
    
        #axs_hist[j].set_title(f"{i_probe}")
    
        if j == 0:
            axs_hist[j].set_xlabel('tau(ms)', fontsize=20)
            axs_hist[j].set_ylabel('Neurons (n)',fontsize=20)
            axs_hist[j].tick_params(axis='y', labelsize=16) 
    
        axs_hist[j].axvline(med, color=area_color[j],linestyle='dashed', linewidth=2)
        axs_hist[j].annotate(f'Median: {10 ** med:.2f} ms', xy=(med, 0.8), xytext=(med - 1,45 ),
                                color=area_color[j])
        
    if monkey_name is not None:
        # Set a suptitle for the entire figure
        fig_hist.suptitle(f'{monkey_name}')
      


#%% main part to run the subfunctions above

#====== Part of the script to run the subfunction TAU_FR_Layerwise for all sessions 
#across layers for each monkey

#==== Need to give information about the monkey name, probe area, and the session names here first before running the script ====#




monkey_name = 'Tommy'
probes = ['PMd', 'M1']

#Read the information file to collect all the neurons and layer information for a particular probe area for a monkey

#excel sheet for PMd/ M1 info of TomyAlmapSpikes2023

data_info = f'/Volumes/work/comco/nandi.n/IntrinsicTimescales/docs/TomyAlmap2023_AreaInfo.csv'

df = pd.read_csv(data_info)

t_comb= {'PMd': [] ,
         'M1' : []}

for p, probe_area in enumerate(probes):
    df_subset =  df[df['Area']== probe_area]
    
    inFolder_TAU = f"/Volumes/work/comco/nandi.n/IntrinsicTimescales/data/Tommy/TomyAlmap2023/TAU_FM_excludeFirstBin_loess_0.1"
    tau = TAU_allLayers(inFolder_TAU)   # Passing probe_number as input to this function so that only the pkl files which belong to this probe give back the tau values
    
    #load the layerwise data from other sessions of Tommy
    
    inFolder_laminarSess = f'/Volumes/work/comco/nandi.n/IntrinsicTimescales/data/Tommy/combinedData_firstBinExcluded_loess_0.1_updatedFinal/TAU_FR_Lat_LayerWise_{probe_area}.pkl'
    
    with open (inFolder_laminarSess,'rb') as file:
        TAU_T = pickle.load(file)
        
        t_allLayers = np.concatenate(((np.array(TAU_T['L2/3'])).flatten(),(np.array(TAU_T['L5'])).flatten(),(np.array(TAU_T['L6'])).flatten()),0)
        
    t_combined =   np.log10(np.concatenate([tau,t_allLayers],axis=0))
    
    t_comb[probe_area] = t_combined
    

TAU = [t_comb['PMd'], t_comb['M1']]

#plotHistogram_PMd_M1(TAU, probes, None,monkey_name)
    
# s_PMd_M1, p_value = kruskal(TAU[0], TAU[1]) 
# significance = 0.05  

# if p_value <significance:
#     print (f'Significant difference of tau between PMd and M1 with p value of {p_value} and statistic {s_PMd_M1}')
# else:
#     print ('No significant difference between tau of PMd and M1')
 

mu_PMd = []
mu_M1 = []

n_iteration = 10000

for k in range(n_iteration):

    subset_PMd = np.random.choice(t_comb['PMd'],90,replace=True)
    subset_M1 = np.random.choice(t_comb['M1'], 90, replace = True) 
    
    mu_PMd.append(np.median(subset_PMd))
    mu_M1.append(np.median(subset_M1))
    
    
    TAU = [np.array(mu_PMd), np.array(mu_M1)]

#plot the distributions of the sampling medians of PMd and M1 and see if they overlap


# Create a new figure and subplots for histograms structure wise
fig_hist, axs_hist = plt.subplots(1, 1, figsize=(8, 3), sharex=True)
area_color = ['darkblue','purple']
# Define the desired x-axis ticks
x_ticks = [1, 2, 3]  # Corresponding to 10, 100, 1000 in non-log scale


# Create a histogram plot for the first row of subplots
for j,i_probe in enumerate(probe_area): 


    axs_hist.hist(TAU[j], bins=20, color=area_color[j], alpha = 0.5, edgecolor='black')

  # #  Set y-axis limits to ensure it ranges from 0 to 1
  #   axs_hist[j].set_ylim(0, 75)



  #   axs_hist[j].set_xticks(x_ticks, [f'{10 ** tick:.0f}' for tick in x_ticks], fontsize= 16)
  #   axs_hist[j].tick_params(axis='y', labelsize=12)

#axs_hist[j].set_title(f"{i_probe}")

    if j == 0:
        axs_hist.set_xlabel('tau(ms)', fontsize=20)
        axs_hist.set_ylabel('Neurons (n)',fontsize=20)
        axs_hist.tick_params(axis='y', labelsize=16)    




#%% combine with Mourad

monkey_name = 'Combined'
probes = ['PMd', 'M1']

#Read the information file to collect all the neurons and layer information for a particular probe area for a monkey

#excel sheet for PMd/ M1 info of TomyAlmapSpikes2023

data_info = f'/Volumes/work/comco/nandi.n/IntrinsicTimescales/docs/TomyAlmap2023_AreaInfo.csv'

df = pd.read_csv(data_info)

t_comb= {'PMd': [] ,
         'M1' : []}

for p, probe_area in enumerate(probes):
    df_subset =  df[df['Area']== probe_area]
    
    inFolder_TAU = f"/Volumes/work/comco/nandi.n/IntrinsicTimescales/data/Tommy/TomyAlmap2023/TAU_FM_excludeFirstBin_loess_0.1"
    tau = TAU_allLayers(inFolder_TAU)   # Passing probe_number as input to this function so that only the pkl files which belong to this probe give back the tau values
    
    #load the layerwise data from other sessions of Tommy
    
    inFolder_laminarSess = f'/Volumes/work/comco/nandi.n/IntrinsicTimescales/data/Tommy/combinedData_firstBinExcluded_loess_0.1_updatedFinal/TAU_FR_Lat_LayerWise_{probe_area}.pkl'
    inFolder_laminarSess_M = f'/Volumes/work/comco/nandi.n/IntrinsicTimescales/data/Mourad/combinedData_firstBinExcluded_loess_0.1_updatedFinal/TAU_FR_Lat_LayerWise_{probe_area}.pkl'
    
    with open (inFolder_laminarSess,'rb') as file:
        TAU_T = pickle.load(file)
        
    with open (inFolder_laminarSess_M,'rb') as file: #combine with Mourad
         TAU_M = pickle.load(file)
        
    t_allLayers = np.concatenate(((np.array(TAU_T['L2/3'])).flatten(),(np.array(TAU_T['L5'])).flatten(),(np.array(TAU_T['L6'])).flatten(),(np.array(TAU_M['L2/3'])).flatten(),(np.array(TAU_M['L5'])).flatten(),(np.array(TAU_M['L6'])).flatten()),0)
        
    t_combined =   np.log10(np.concatenate([tau,t_allLayers],axis=0))
    
    t_comb[probe_area] = t_combined
    

TAU = [t_comb['PMd'], t_comb['M1']]

plotHistogram_PMd_M1(TAU, probes, None,monkey_name)
    
s_PMd_M1, p_value = kruskal(TAU[0], TAU[1]) 
significance = 0.05  

if p_value <significance:
    print (f'Significant difference of tau between PMd and M1 with p value of {p_value} and statistic {s_PMd_M1}')
else:
    print ('No significant difference between tau of PMd and M1')


#%%
#do the AP gradient with timescales just for Tomy











