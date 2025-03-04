#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 12:28:11 2024

This script is used to compare the tau distribution between PMd and M1 layrwise either using combined data acros two monkeys
or single monkey 

The main function calls 2 subfunctions - compare_PMd_M1_layerwise  , plotHistogram_layerwise 
 
Inputs of the main function -

combine_monkey  = 0 or 1 


Before running this function, the AlighTAU_FR_Layerwise.py and then combineMonkeyData.py needs to be run since it
calls the data processed by the above two functions 

@author: nilanjana
"""

#import modules
import os, pickle
import matplotlib.pyplot as plt
from scipy.stats import kruskal, bootstrap
import numpy as np
import seaborn as sns
import pandas as pd
#import matplotlib.colors as mcolors

def compare_PMd_M1_layerWise(inFolder,plot_combineData,probe_area):
    '''
    

    Parameters
    ----------
    inFolder : path for combined monkey data/ single monkey data for PMd/M1
    plot_combineData : int - 0 or 1
    probe_area : list - [PMd, M1]

    Returns
    -------
    TAU_PMd : array
    TAU_M1 : array

    '''
    layers =['L2/3','L5','L6']
    
    for probe in probe_area:  
        log_tau = []
        
        if plot_combineData ==1:       
            data = inFolder
            with open(data,'rb') as file:
                combined_TAU = pickle.load(file) # dictionary containing combined data across 2 monkeys in the 3 layers and all layers
            # Iterate through each layer for the histogram
            for j, layer in enumerate(layers):           
                subset = np.concatenate(combined_TAU[probe][layer])
                log_tau.append((subset)) #removed np.log10                 
        else:
            data = f'{inFolder}/TAU_FR_Lat_{probe}.pkl'           
            with open(data,'rb') as file :
                tau = pickle.load(file)
            
            for j, layer in enumerate(layers):
                subset = np.concatenate(tau[f'{layer}'])
                log_tau.append((subset)) #removed np.log10 
            
        if probe == 'PMd':
            TAU_PMd = log_tau
        else :
            TAU_M1 = log_tau
       
    return TAU_PMd, TAU_M1


# Function to add star and bracket if p-value is significant
def add_significance(ax, x1, x2, y_max, p_value):
    if p_value < 0.05:  # Check if p-value is significant
        y_bracket = y_max * 1.1  # Position the bracket slightly above the violins
        y_star = y_bracket * 1.02  # Position the star slightly above the bracket
        
        # Draw the bracket
        ax.plot([x1, x1, x2, x2], [y_bracket, y_star, y_star, y_bracket], color='black', linewidth=1.5)        
        # Add the star
        ax.text((x1 + x2) / 2, y_star, '*', fontsize=20, ha='center', va='bottom', color='black')


def plotHistogram_layerWise(TAU, probe_area, outFolder,monkey_name):
    '''
    

    Parameters
    ----------
    TAU : list of length 3 - arrays containing tau at L2/3, L5 and L6
    probe_area : PMd/M1
    outFolder : path to store the hisotgram plot
    monkey_name : str - Mourad/Tomy

    Returns
    -------
    None.

    '''
    # Create a new figure and subplots for histograms layerwise
       
    layers =['L2/3','L5','L6']    
        # Define base colors
    # pmd_base_color = '#B1C086' #"#00008B"  # Dark Blue
    # m1_base_color =  '#EBCCFF'  #"#DA70D6"   # Orchid
    
    # # Define transparency levels for layers
    # alpha_values = [0.05, 0.6, 1.0]  # Adjust these for L2/3, L5, L6
    
    # # Convert hex to RGBA with transparency
    # pmd_layer_colors = [mcolors.to_rgba(pmd_base_color, alpha) for alpha in alpha_values]
    # m1_layer_colors = [mcolors.to_rgba(m1_base_color, alpha) for alpha in alpha_values]
    
    for p, area in enumerate(probe_area): 
        
        if area == 'PMd':
            layer_colors =['#F0F0D7', '#D0DDD0', '#AAB99A'] #['#E5E3C9','#B4CFB0','#94B49F']  #pmd_layer_colors #
        else:
            layer_colors = ['#FFDFEF', '#EABDE6', '#E9A8F2']#'#D69ADE']#m1_layer_colors
        
        # Create a new figure and subplots for histograms structure wise
        fig, ax = plt.subplots(figsize=(12,8))
        
        data = pd.DataFrame({
            'tau' : np.concatenate(TAU[p]),
            'layer' : ['L2/3'] * len(TAU[p][0]) + ['L5'] * len(TAU[p][1]) + ['L6']*len(TAU[p][2])
            })
                
             # Create violin plot with both areas side by side
        sns.violinplot(x='layer', y='tau', data=data, palette=layer_colors,hue = None, inner='box',ax=ax)
        
        #Calculate medians and add custom legend entries
        for i, layer in enumerate(layers):
            # Calculate the median
            med = np.median(data[data['layer'] == layer]['tau'])
            
            # Add a custom legend entry for the median
            ax.plot([], [], color=layer_colors[i], linestyle='dashed', linewidth=2, 
                    label=f'{area} Median: {med:.2f}')
        
        #Kruskal Test
        s, p_val= kruskal(TAU[p][0],TAU[p][1],TAU[p][2])
     
        if p_val <0.05:
            print (f'Significant difference of tau between layers in {area} in {monkey_name}')
            
            #PostHoc Analysis if there is difference in median
            
            
            # Perform Kruskal-Wallis test for each pair of groups
            s_1, p_L2_L5 = kruskal(TAU[p][0], TAU[p][1])
            s_2, p_L2_L6 = kruskal(TAU[p][0], TAU[p][2])
            s_3, p_L5_L6 = kruskal(TAU[p][1], TAU[p][2])
            
            # Print or analyze the results for each pair of groups
            print(f'Group L2/3 vs. Group L5:  p_val = {p_L2_L5:.2f} in {monkey_name}')
            print(f'Group L2 vs. Group L6:  p_val = {p_L2_L6:.2f} in {monkey_name}')
            print(f'Group L5 vs. Group L6:  p_val = {p_L5_L6:.2f} in {monkey_name}')
            
               # Get the maximum y-value for positioning the bracket
            y_max = max(data['tau'])
            # Add significance markers for each pair of groups
            add_significance(ax, x1=0, x2=1, y_max=y_max, p_value=p_L2_L5)  # L2/3 vs. L5
            add_significance(ax, x1=0, x2=2, y_max=y_max, p_value=p_L2_L6)  # L2/3 vs. L6
            add_significance(ax, x1=1, x2=2, y_max=y_max, p_value=p_L5_L6)  # L5 vs. L6
        
        #     # Define the positions for the bracket and star
        # x1, x2 = 0, 1  # x-positions of the two violinss
        # y_max = max(data['tau'])  # Maximum y-value for positioning the bracket
        # y_bracket = y_max * 1.1  # Position the bracket slightly above the violins
        # y_star = y_bracket * 1.02  # Position the star slightly above the bracket
        
        # # Draw the bracket
        # ax.plot([x1, x1, x2, x2], [y_bracket, y_star, y_star, y_bracket], color='black', linewidth=1.5)
        
        # # Add the star
        # ax.text((x1 + x2) / 2, y_star, '*', fontsize=20, ha='center', va='bottom', color='black')
        
        ax.set_xlabel('layer', fontsize=20)
        ax.set_ylabel('tau (ms)', fontsize=20)
        ax.tick_params(axis='both', labelsize=18)
        ax.legend(loc='upper right', fontsize=16)    
        
        if monkey_name is not None:
            # Set a suptitle for the entire figure
            fig.suptitle(f'{monkey_name}_{area}', fontsize=20)
    
        # Adjust layout to prevent overlap
        plt.tight_layout()
        #plt.savefig(os.path.join(outFolder,f'{monkey_name}_{probe_area}_layerWise'))
    
    
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
        
    # monkey = ['Combined','Mourad', 'Tomy']
    # probe_area = ['PMd', 'M1']  
    monkey = ['Combined', 'Mourad']
    probe_area = ['PMd', 'M1']    
    
    for m, monkey_name in enumerate(monkey):       
        if monkey_name == 'Combined':
            combine_monkey =1 
            outFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/results/SUA/combinedMonkey' 
            if not os.path.exists(outFolder):
                os.makedirs(outFolder)
                
            inFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/combinedMonkey/combined_Laminar.pkl'
            # # -------layerwise organize the lists for PMd/M1
            
            
        else:
            combine_monkey = 0            
            #run for loop for individual monkey               
            outFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/results/SUA/{monkey_name}' 
            if not os.path.exists(outFolder):
                os.makedirs(outFolder)
            
            inFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/{monkey_name}/combined_Layerwise'
            TAU_PMd, TAU_M1 = compare_PMd_M1_layerWise(inFolder,combine_monkey,probe_area)      
       
        TAU_PMd, TAU_M1 = compare_PMd_M1_layerWise(inFolder,combine_monkey,probe_area)       
        plotHistogram_layerWise(([TAU_PMd,TAU_M1]), probe_area, outFolder,monkey_name)
        
#printing the layer statistics
        for probe in probe_area:   
            if probe == 'PMd':                   
                TAU = TAU_PMd
            else:
                TAU = TAU_M1
         
            # #Kruskal Test
            # s, p_val= kruskal(TAU[0],TAU[1],TAU[2])
         
            # if p_val <0.05:
            #     print (f'Significant difference of tau between layers in {probe} in {monkey_name}')
                
            #     #PostHoc Analysis if there is difference in median
                
            #     # Perform Kruskal-Wallis test for each pair of groups
            #     s_1, p_L2_L5 = kruskal(TAU[0], TAU[1])
            #     s_2, p_L2_L6 = kruskal(TAU[0], TAU[2])
            #     s_3, p_L5_L6 = kruskal(TAU[1], TAU[2])
                
            #     # Print or analyze the results for each pair of groups
            #     print(f'Group L2/3 vs. Group L5:  p_val = {p_L2_L5:.2f} in {monkey_name}')
            #     print(f'Group L2 vs. Group L6:  p_val = {p_L2_L6:.2f} in {monkey_name}')
            #     print(f'Group L5 vs. Group L6:  p_val = {p_L5_L6:.2f} in {monkey_name}')        
                    
                    
                    
                    
                    
                    
            
    