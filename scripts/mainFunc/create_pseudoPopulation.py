#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 17:28:49 2024

This script is used to create snippets of time series data. The duration is (~5min) recording time for each snippet f
from the entire recording session if the session lasted for more than 5 minutes. 

This is done to increase the number of tau values estimated, if we have long sessions. Eg, a 15 min session can be divided 
into 3 snippets and we can compute tau on each of these timeseries data give us 3 tau estimates instead of 1

This was done to check if we increase the number of tau values especially when we compare them for each monkey separately
for each layer


@author: nilanjana
"""

#monkey loop 

#loop the sessions

#inside each session choose the neuron 

#check the neuron was stable for more than 5 then divide, with each snippet having 5 min

#else take the entire time if length is approx 5 , 

#else skip the neuron 

def run_computeTAU(inFolder,outFolder,lag_limit,excludeFirstBin):
        

    # Create folder if it does not exist
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
        
    #import data
    # Search for mat files in the folder
    file_names = [filename for filename in os.listdir(inFolder) if filename.endswith(".mat")]
      
        # Loop through the file names
    for i, file_name in enumerate(file_names, start=1):
       
        print (f"Processing Session {file_name}")
        # Load time series data
        file_path = os.path.join(inFolder,file_name)
        
        data =  open_matlab_behaviour(file_path)       
        d = data['ts']
        # Extract the last string name from the file name
        last_name = os.path.splitext(os.path.basename(file_name))[0]
        
        # Extract the information after "unit" using regular expression and string manipulation
        extracted_info = re.sub(r'.*unit(\d+).*', r'\1', last_name)
        
        # Convert the extracted info to a numeric variable        
        if extracted_info.isdigit() :
            unit = int(extracted_info) 
        else :
            print ("Cannot Extract Unit Information")
            unit = 0; # assigning by default
        
        len_rec = d[-1] - d[0] #length of time for which the neuron was recorded
        chunk_size = 10 * 60 * 1000  # 10 minutes in ms (300,000 ms)
        n_iter = 0
        if len_rec > chunk_size:  
            
            start_idx = 0  # Start from the beginning  index          
            while (start_idx < len(d)) and ((d[-1] - d[start_idx]) >= chunk_size): 
                
            # Find the time at `start_idx + chunk_size`
                target_time = d[start_idx] + chunk_size
                
                # Find the index of the closest time point to `target_time`
                end_idx = (np.abs(d - target_time)).argmin()
                
                # Perform your next computation here for the current chunk (d[start_idx:end_idx+1])
                # Run main function
                d_subset = d[start_idx:end_idx]
                [out_fit, out_autoc] = mainFun_computeTAU(d_subset,unit,lag_limit,excludeFirstBin)
                
                # For the next iteration inside the while loop -- I need to put this here since I have a continue 
            #argument afterwards - For neurons where the SVD does not converge, instead of breaking the loop, the continue
        #will run the while loop again, but the new start_idx and n_iter needs to be updated 
        #If the n_iter and start_idx are placed afterwards then they do not get updated once the disruption occurs
                
                start_idx = end_idx + 1
                n_iter+=1
                
                if out_fit is None or out_autoc is None:
                    print(f"Skipping {n_iter} due to computation issue of SVD converge")
                    continue #skips this iteration and moves on to the next iteration
                    
                # Create filename for saving output
                output_filename = os.path.join(outFolder, f"{last_name}_SNIPPET_{n_iter}.pkl")
                 
                #Save output using pickle             
                with open(output_filename, "wb") as output_file:
                    pickle.dump(out_fit, output_file)
                    pickle.dump(out_autoc, output_file)   
                    
                # convert dataframe to array to save in mat 
                
                # arr1 = np.array(out_fit.to_records(index=False));
                # arr2 = np.array(out_autoc.to_records(index=False));
                
                # data_dic=  {"autoc" : arr2, "fit_val": arr1}
                # filename2 = os.path.join(outFolder, f"{last_name}.mat")
                # savemat(filename2, data_dic)
                
                
        else:
             print("Recording length is less than 10 minutes, skipping the neuron.")    
            
            
    if n_iter > 0:
        used_sess = 1
    else:
        used_sess = 0
        
        
    return used_sess  
               

if __name__=='__main__':
    
    import os        
    
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
        
    
    import sys,re, pickle
     #import pandas as pd

    if server == '/envau':
        sys.path.append('/envau/work/comco/nandi.n/IntrinsicTimescales/scripts/Processing_Methods')
    else:
        sys.path.append('/Users/nilanjana/Documents/LabWorkBench/IntrinsicTimescales/scripts/Processing_Methods')
    
    from computeTAU_FontanierMethod import mainFun_computeTAU
    from open_matlabData import open_matlab_behaviour     
    # from scipy.io import loadmat
    import numpy as np
    from scipy.io import savemat    
   
    monkey_names = ['Tomy' , 'Mourad']
    #monkey = 'Mourad'

    lag_limit = [1000]
    excludeFirstBin = 1 
    
    sess_m = []
    for monkey in monkey_names:
    
        if monkey == 'Tomy':
            session_names =['t140924003','t140924006','t140925001','t140926002','t140930001','t141001001','t141008001',
                            't141010003','t141119001','t150122001','t150123001','t150128001','t150129002','t150204001',
                            't150205004','t150210001','t150211003','t150212001','t150218001','t150219003','t150303002',
                            't150319003','t150324002','t150327002','t150327003','t150415002','t150416002','t150423002',
                            't150430002','t150520003','t150716001']
            #session_names = ['TomyAlmap2023']
        
        elif monkey == 'Mourad':
            session_names= ['Mo180328001','Mo180330001','Mo180405001','Mo180405004','Mo180411001','Mo180412002',
                            'Mo180411001','Mo180418002','Mo180419003','Mo180503002', 'Mo180523002','Mo180524003', 
                            'Mo180525003','Mo180530004','Mo180531002','Mo180601001','Mo180614002','Mo180614006',
                            'Mo180615002','Mo180615005','Mo180426004', 'Mo180427003', 'Mo180619002','Mo180622002',
                            'Mo180626003', 'Mo180627003','Mo180629005', 'Mo180703003','Mo180704003', 'Mo180705002'
                              'Mo180706002', 'Mo180710002','Mo180711004','Mo180712006',
                                  ]
           
            #new sessions updated after my sorting
            
            # ['Mo180426004', 'Mo180427003', 'Mo180619002','Mo180622002',
            # 'Mo180626003', 'Mo180627003','Mo180629005', 'Mo180703003','Mo180704003', 'Mo180705002'
            #  'Mo180706002', 'Mo180710002','Mo180711004','Mo180712006']
            
        
        for l in lag_limit:
            s = 0
            for i_session in session_names:
                session = i_session   
                inFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{session}/modifiedData"
                outFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{session}/TAU_FM_PseudoPopulation"
                sess_used = run_computeTAU(inFolder,outFolder,l,excludeFirstBin)
                
                if sess_used > 0 :
                    s+=1
            
        sess_m.append(s)
            
    print('Success computing all fits -- sessions used in {monkey_names[0]}: {sess_m[0]}')
    print('Success computing all fits -- sessions used in {monkey_names[1]}: {sess_m[1]}')
            
             



