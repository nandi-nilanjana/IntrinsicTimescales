#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 15:56:24 2023

Modifying this function to only combine the data across two monkeys. Data is combined for the same area (probe_name) across
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

#combining data across monkeys and creating figure 1  and figure 2 for GDR



#%%


def combineData(dataFolder,combineDataFolderExt,monkey_name,probe_area): #subfunction to combine data across two monkeys layerwise
   
    # Create dictionaries to store combined data for each probe area
    combined_TAU = {area: {} for area in probe_area}
    combined_FR = {area: {} for area in probe_area}
    combined_LAT = {area: {} for area in probe_area}
    
    for i, i_probe in enumerate(probe_area):
        
        for i_monkey in monkey_name:
                            
            data = f"{dataFolder}/{i_monkey}/{combineDataFolderExt}/TAU_FR_Lat_LayerWise_{i_probe}.pkl"
        
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
    
    
    #modules
    import pickle
    import os
    
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
        
    inFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/data'   
    combinedataFolder = 'combinedData_firstBinExcluded_loess_0.1_updatedFinal'
    
    combined_TAU, combined_FR, combined_LAT = combineData(inFolder,combinedataFolder,monkey_names,probe_areas)
    
    outFolder = f'{server}/work/comco/nandi.n/IntrinsicTimescales/data/combined_Mourad_Tomy_loess_0.1_updatedFinal'
    
    if not os.path.exists(outFolder): 
        os.makedirs(outFolder)
    
    outputFileName =  f"{outFolder}/combined_PMd_M1_excludeFirstBin.pkl"
        
    with open(outputFileName,"wb") as file :
        pickle.dump(combined_TAU,file)
        pickle.dump(combined_FR,file)
        pickle.dump(combined_LAT,file)
        
        
    
    
    



