#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 15:56:24 2023

Modifying this function to only combine the laminar data across two monkeys. Data is combined for the same area (probe_name) across
layers 2/3 , L5, L6 and also for all layers 

Input - 
dataFolder - main path to the folder where the data is stores which is generally /server/work/...

combineDataFolderExt - entension for the combined data folder for each monkey - since I saved combined data for TAU -  combined_Mourad_Tomy_loess_0.1_updatedFinal

monkey_name - the monkey names to combine the data aross them since the layerwise tau, fr and lat data are saved 
in each monkey folder 

probe_area - the probe area to combine across - PMd/M1
This function will only run when the combineData is saved . If there is no folder, then you need to run AlignTAU_FR_LayerWise first

Output

# dictionary structure with keys - PMd  and M1. Inside each key, you have 4 keys all_layers, L2/3, L5 and L6. 
The data inside each layer is a list (length of the list is the length combined across sessions and across the 2 monkeys)

Just use np.concatenate(data['PMd'/'M1']['layer_key']) to get an array combined across all the tau data from each list



@author: nilanjana
"""



#import modules
import pickle
import os
import pandas as pd
import numpy as np

#%%#combining data across monkeys 

def combineData(dataFolder,combineDataFolderExt,monkey_name,probe_area): #subfunction to combine data across two monkeys layerwise
   
    # Create dictionaries to store combined data for each probe area
    combined_TAU = {area: {} for area in probe_area}
    combined_FR = {area: {} for area in probe_area}
    combined_LAT = {area: {} for area in probe_area}
    
    for i, i_probe in enumerate(probe_area):
        
        for i_monkey in monkey_name:
                            
            data = f"{dataFolder}/{i_monkey}/{combineDataFolderExt}/TAU_FR_Lat_{i_probe}.pkl"
        
            with open(data, 'rb') as file:
                TAU = pickle.load(file)
                FR = pickle.load(file)
                LAT = pickle.load(file)
                
            #if combineLayers == 1:
            if "all_layers" not in combined_TAU[i_probe]:
                combined_TAU[i_probe]["all_layers"] = []
                combined_FR[i_probe]["all_layers"] = []
                combined_LAT[i_probe]["all_layers"] = []
                
    
            for key in TAU:
                if key not in combined_TAU[i_probe]:
                    combined_TAU[i_probe][key] = []
                    combined_FR[i_probe][key] = []
                    combined_LAT[i_probe][key] = []
                    
                
                combined_TAU[i_probe]["all_layers"].extend(TAU[key]) #combining across all layers 
                combined_FR[i_probe]["all_layers"].extend(FR[key])
                combined_LAT[i_probe]["all_layers"].extend(LAT[key])
                
                
                combined_TAU[i_probe][key].extend(TAU[key]) #combining across specific layer for each monkey
                combined_FR[i_probe][key].extend(FR[key])
                combined_LAT[i_probe][key].extend(LAT[key])
                
                
    return combined_TAU, combined_FR, combined_LAT


if __name__ == '__main__':
            
     #give monkey name 
    monkey_names  = ['Tomy', 'Mourad']
     
     #which area  M1/PMd
    probe_areas = ['PMd','M1']
    
    #server - VPN or inside the server
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
        
    inFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA'   
    combinedataFolder = 'combined_LayerWise'
    
    combined_TAU, combined_FR, combined_LAT = combineData(inFolder,combinedataFolder,monkey_names,probe_areas)
    
    outFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/combinedMonkey'
    
    if not os.path.exists(outFolder): 
        os.makedirs(outFolder)
    
    outputFileName =  f"{outFolder}/combined_Laminar.pkl"
        
    with open(outputFileName,"wb") as file :
        pickle.dump(combined_TAU,file)
        pickle.dump(combined_FR,file)
        pickle.dump(combined_LAT,file)
        
        
#%% 
'''
Combining the tau from single tip electrodes from Tomy to save the data just for Tomy for PMd,M1
Also combining the ST Tomy data with Mourad to save all values in the combinedMonkey folder

'''


#========== Subfunction to concatenate the tau, lat and the firing rate of neurons across layers 
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


if __name__=='__main__':
    
    imonkey = ['Tomy','combinedMonkey']
    #server - VPN or inside the server
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
    
    for monkey_name in imonkey:
        probes = ['PMd', 'M1']
        #excel sheet for PMd/ M1 info of TomyAlmapSpikes2023
        data_info = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/docs/TomyAlmap2023_AreaInfo.csv'
        
        df = pd.read_csv(data_info)    
     
        t_comb= {area: {} for area in probes}
     
        for p, probe_area in enumerate(probes):
            df_subset =  df[df['Area']== probe_area]
            
            inFolder_SingleTip = f"{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/Tomy/TomyAlmap2023/TAU_FM"
            tau = TAU_allLayers(inFolder_SingleTip)   
            
            if monkey_name == 'Tomy':
                outFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/{monkey_name}/ST'
                
                #load the layerwise data from other sessions of Tommy    
                inFolder_laminarSess = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/{monkey_name}/combined_Layerwise/TAU_FR_Lat_{probe_area}.pkl'            
                with open (inFolder_laminarSess,'rb') as file:
                    TAU = pickle.load(file)        
                    t_allLayers = np.concatenate(((np.concatenate(TAU['L2/3'])),(np.concatenate(TAU['L5'])),(np.concatenate(TAU['L6']))))
            
            elif monkey_name == 'combinedMonkey':
                outFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/{monkey_name}'
                inFolder_laminarSess = f'{server}/work/comco/nandi.n/IntrinsicTimescales/Paper/data/SUA/{monkey_name}/combined_Laminar.pkl'            
                with open (inFolder_laminarSess,'rb') as file:
                    TAU = pickle.load(file)        
                    t_allLayers = np.concatenate(((np.concatenate(TAU[probe_area]['L2/3'])),(np.concatenate(TAU[probe_area]['L5'])),(np.concatenate(TAU[probe_area]['L6']))))
               
            t_combined =   [tau,t_allLayers]          
            t_comb[probe_area]['all_layers'] = t_combined
            
        if not os.path.exists(outFolder): 
            os.makedirs(outFolder)
        
        outputFileName =  f"{outFolder}/combined_ST.pkl"
            
        with open(outputFileName,"wb") as file :
            pickle.dump(t_comb,file)
            
            
    
        
        
