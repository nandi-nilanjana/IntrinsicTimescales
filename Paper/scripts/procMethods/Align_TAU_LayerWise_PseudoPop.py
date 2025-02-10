#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:45:58 2024

@author: nilanjana

Script created based on AlignTAU_FR_LayerWise.py but here instead we want to only concatenate TAU and LAT for 
all the neural pseudo populations estimated by dividing each long session into 10min snippets 

Here apart from excluding the fits with sd > 150 , we also exclude fits with ISI intervals (<2000 in pre500 and post500)
time lag epochs. This was suggested after manually inspecting all the fits . 


"""


#========== Subfunction to concatenate the tau, lat and the firing rate of neurons across layers for a particular session and monkey

#import modules
import os
import pickle
import re
import pandas as pd
import numpy as np

def TAU_FR_Layerwise(data_path_TAU,probe_name,ranges,all_layer,ISI_counts):
    '''
    

    Parameters
    ----------
    data_path_TAU : path / directory containing tau values stored for each neuron inside a specific session    
    probe_name :  here an integer - probe 1 or 2 taken from the excel list selected as a subset for PMd/M1
    ranges : for each session the corresponding contacts for the layer ranges
    all_layer : list containing the 3 layer names
    ISI_counts : int - total number of ISI counts to be present in the pre500 and post500 time lag period
    Returns
    -------
    tau_array : list of arrays
    lat_array : list of arrays 
    
    '''
   
    tau_array = {f"{layers[i]}": [] for i in range(len(all_layer))}
    lat_array = {f"{layers[i]}": [] for i in range(len(all_layer))}
    
    # Search for pickle files in the folder for that session
    
    file_names_TAU = [filename for filename in os.listdir(data_path_TAU) if filename.endswith(".pkl")]
    
    
    # Load and process data for each session here
    
    for i, file_name in enumerate(file_names_TAU):
        
        file_path_TAU = os.path.join(data_path_TAU, file_name) 
        
        #first check if you have the corresponding file in the Firing rate directory
        if os.path.exists(file_path_TAU):                    
            if probe_name in file_name:  # Check if the probe name is in the file name belonging to both TAU and FR directory
            
                with open(file_path_TAU,'rb') as file: #load TAU Data
                    fit_val = pickle.load(file)
                    autoc = pickle.load(file)
             
                if fit_val['bestfit'][0] == 'mono':                    
                    mono_fit = fit_val[fit_val['fit_type'] == 'mono'] #find the mono fit                    
                    l_pre500 = autoc['mids']<=500
                    l_post500 = autoc['mids']>=500
                    
                    s1= autoc['counts'][l_pre500].sum()
                    s2 = autoc['counts'][l_post500].sum()
                    
                    if (mono_fit['TAU_sd'].values <150) and (s1>ISI_counts) and (s2>ISI_counts):
                    
                        #find the contact info
                        last_name = os.path.splitext(os.path.basename(file_name))[0]
                          
                        # Extract the information after "contact" using regular expression and string manipulation
                        contact_info = re.sub(r'.*contact(\d+).*', r'\1', last_name)
                        contact = int(contact_info)
                        
                        # Iterate through the ranges and assign values to arrays
                        for i, range_str in enumerate(ranges):
                            start, end = map(int, range_str.split('-'))
                            if start <= contact <= end:
                                t_val = mono_fit['TAU'].iloc[0]
                                lat_val = mono_fit['peak1_lat'].iloc[0]                                
                                tau_array[f"{layers[i]}"].append(t_val)
                                lat_array[f"{layers[i]}"].append(lat_val)
                                
   
    return tau_array,lat_array #returns the tau values for all the layers from this probe and session 


# Define a function to check if a value is empty
def is_empty(value):
    return value is None or (isinstance(value, (str, list)) and len(value) == 0)



if __name__ == '__main__' :
    
    '''
     Part of the script to run the subfunction TAU_FR_Layerwise for all sessions #across layers for each monkey
     Need to give information about the monkey name, probe area, and the session names here first before running the script 
    '''
    
    #give monkey name
    monkey = ['Mourad', 'Tomy']
        
    #which area to combine for    
    probe = ['PMd', 'M1']
        
    # Info required to create a dictionary of arrays for each range later in the code
    layers = ['L2/3', 'L5' , 'L6']
    
    #server - VPN or inside the server
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
        
    ISI_thres_count = 2000 #total number of ISI counts required before/after 500 time lag perioid to accept that fit
    
    for monkey_name in monkey: #Mourad/Tomy
        for probe_area in probe: #PMd/M1
                   
            if monkey_name == 'Tomy':                          
            
                session = ['t140918002','t140924003','t140924006','t140925001','t140926002','t140930001','t141001001',
                           't141008001','t141010003','t141119001','t150122001','t150123001','t150128001','t150129002',
                           't150204001','t150205004','t150210001','t150211003','t150212001','t150218001','t150219003',
                           't150303002','t150319003','t150324002','t150327002','t150327003','t150415002','t150416002',
                           't150423002','t150430002','t150520003','t150716001']
            
            elif monkey_name == 'Mourad':
            
                #Session Name - All sessions with qulaity 1 neurons
                
                session =['Mo180328001','Mo180330001','Mo180405001','Mo180405004','Mo180411001','Mo180412002','Mo180418002',
                          'Mo180419003','Mo180503002','Mo180523002','Mo180524003','Mo180525003','Mo180530004','Mo180531002',
                          'Mo180601001','Mo180614002','Mo180614006','Mo180615002','Mo180615005''Mo180426004', 'Mo180427003', 
                          'Mo180619002', 'Mo180622002', 'Mo180626003', 'Mo180627003', 'Mo180629005', 'Mo180703003',
                          'Mo180704003', 'Mo180705002', 'Mo180706002', 'Mo180710002','Mo180711004','Mo180712006']
                
                
                #quality 4 sessions are below
                # session = ['Mo180330001', 'Mo180405001', 'Mo180405004', 'Mo180410001', 
                #                  'Mo180410002', 'Mo180418002', 'Mo180419003', 'Mo180426004', 
                #                  'Mo180503002', 'Mo180524003', 'Mo180525003', 'Mo180530004', 'Mo180531005', 
                #                  'Mo180601001', 'Mo180614006', 'Mo180615002', 'Mo180615005',  'Mo180619002',
                #                  'Mo180620004', 'Mo180622002', 'Mo180626003', 'Mo180627003', 'Mo180629005', 
                #                  'Mo180703003', 'Mo180704003', 'Mo180705002', 'Mo180706002', 'Mo180710002', 
                #                  'Mo180711004,' 'Mo180712006']
            
                   
            
            
            #excel sheet with layer inormation for each monkey
              
            data_info = f'{server}/work/comco/nandi.n/IntrinsicTimescales/docs/{monkey_name}LaminarInfo.xlsx'
            
            df = pd.read_excel(data_info)
            df_subset =  df[df['Probe_Area']== probe_area ] #first take the subset of sessions and the probe numbers 
            #(somtimes 2 in each session) belonging to either PMd/M1
            
            # Create an empty list to store the combined data from all sessions
            combined_TAU = []
            combined_lat = []
            
            # Iterate through each session
            for i, sess in enumerate(session):
                #print (sess)
                # First see if the particular session exists or not for M1/PMd from the list
                df_sess = df_subset[df_subset['Session']==sess]
                df_sess = df_sess.reset_index(drop=True)
                 
                if not df_sess.empty:
             
            # Read the data and the laminar info
            # Data path
                   
                    inFolder_TAU = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey_name}/{sess}/TAU_FM_PseudoPopulation"                  
                    print(f"Processing Session {sess}")
                     
                    n_probes = len(df_sess) # Some sessions have M1/PMd in both the probes
                    
                    for probe_idx in range(n_probes):
                        probe_info = df_sess.iloc[probe_idx]
                        probe_number = int(probe_info['Probe_Num'])
                        l = probe_info['Layer_Info'] # Layer info in that probe for the particular session
                        
                        # Use regular expression to find all range patterns in the input string
                        range_pattern = r'(\d+-\d+)'
                        range_l = re.findall(range_pattern, l) # Gives the ranges for each layer in each session and probe
                         
                        # Segregate tau values from each session into different arrays
                        probe_name = f"probe{probe_number}"

# Passing probe_number as input to this function so that only the pkl files which belong to this probe give back the tau values
                        tau_l,lat_l= TAU_FR_Layerwise(inFolder_TAU,probe_name, range_l, layers,ISI_thres_count)   
                      
                        combined_TAU.append(tau_l) #Appen tau from that session and probe 
                        combined_lat.append(lat_l)
                        
                    else:
                        print(f"Skipping Session {sess}")    
            
            p = pd.DataFrame(combined_TAU) 
            p2 = pd.DataFrame(combined_lat)
             
            outFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey_name}/combinedData_TAU_PseudoPop"      
            if not os.path.exists(outFolder):
                os.makedirs(outFolder)                
                
            output_filename = f"TAU_Lat_LayerWise_{probe_area}.pkl"
            outfile =os.path.join(outFolder,output_filename)
            
            with open(outfile, "wb") as output_file:
                pickle.dump(p, output_file)
                pickle.dump(p2, output_file) 
               
                
            print('Successfully completed combining across all sessions')
