#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:44:02 2024

@author: nilanjana
"""

"""
Always run the scripts on VPN since data is accessed and always stored inside the server 

This function can compute the mean firing rate of a neuron. For this it uses the subfunction compute_FR in Processing Methods folder

Inputs-

arg1 - path to access the spike data
arg2 - path to save the firing rate compute for each neuron
arg3 - akes input as either 1 or 0. 1 means the data over the session will be aligned trialwise

For now since the input data accessed from inFolder only has 3 info - the time series data from the whole session, start times of the Selection Onset and start times of Go signal
so if we want to align each trial will only be cut from start time of selection onset till start of Go signal
To change we need to save the mat files in the required structure

The open_matlab_behavior inside this function is currently reading 3 indices

first index - the entire spike data
second index - denotes the start times from where each trial will be cut
third index  - denotes the end times till will teh trial will be cut


Output - mean firing rate for each neuron in pickle file  even for aligned data. 
For aligned data, first FR is computed and then the mean of that is returned as an output and saved 



"""

def run_compute_FR(inFolder,outFolder,alignData,eventMarker):
    
    #modules import
    
    import sys, pickle, os
    sys.path.append('/Users/nilanjana/Documents/LabWorkBench/IntrinsicTimescales/scripts/Processing_Methods')
    from compute_FR import compute_firingRate
    from open_matlabData import open_matlab_behaviour
 
    
    # Create folder if it does not exist
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)


     #import data
     # Search for CSV files in the folder
    file_names = [filename for filename in os.listdir(inFolder) if filename.endswith(".mat")]
     
     # Loop through the file names
    for i, file_name in enumerate(file_names):
        
         
         print (f"Processing Session {file_name}")
         # Load time series data
         file_path = os.path.join(inFolder,file_name)
         
         data =  open_matlab_behaviour(file_path)
         l = list(data.keys())
         #modified on June 10 to extract firing rates 2s before preSC1 and 2 sec before preGo
         d = data[l[0]]
         
         # t_st= data[l[1]]
         # t_end= data[l[2]]
         
         
         if eventMarker  == 'preSC1' :
             #depends on what is stored in the mat data file 
             t_st= data[l[1]]-2000 #start time to cut the trial; keys -  1 - Sc1 and 2 - Go_time
             t_end= data[l[1]] #stop time to cut the trial 
             
         
         elif eventMarker == 'preGo':
             #depends on what is stored in the mat data file 
             t_st= data[l[2]]-2000 #start time to cut the trial; keys -  1 - Sc1 and 2 - Go_time
             t_end= data[l[2]] #stop time to cut the trial 
             
         elif eventMarker == 'preSC2':
             t_st= data[l[1]]-1000
             t_end= data[l[1]]
             
         elif eventMarker == 'full':
             t_st= data[l[1]]
             t_end= data[l[2]]
            
            
         
         
         # Extract the last string name from the file name
         last_name = os.path.splitext(os.path.basename(file_name))[0]
    
    
         # Run main function
         if alignData == 1:
         
             mu_FR = compute_firingRate(1,d,t_st,t_end) 
         
         else:
             
             mu_FR = compute_firingRate(0,d,'','') 
             
    
         # Create filename for saving output
         output_filename = os.path.join(outFolder, f"{last_name}.pkl")
          
         #Save output using pickle (or another preferred format)
         # Save output using pickle (or another preferred format)
         with open(output_filename, "wb") as output_file:
             pickle.dump(mu_FR, output_file)
                 
                 
                 
if __name__=='__main__':
   #monkey = 'Tomy'
    monkey = 'Mourad'

    import os 
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
       
    if monkey == 'Tomy':
        
    
        session_names =['t140924003','t140924006','t140925001','t140926002','t140930001','t141001001','t141008001','t141010003','t141119001','t150122001',
                        't150123001','t150128001','t150129002','t150204001','t150205004','t150210001','t150211003','t150212001','t150218001','t150219003',
                        't150303002','t150319003','t150324002','t150327002','t150327003','t150415002','t150416002','t150423002','t150430002','t150520003',
                        't150716001']
    if monkey == 'Mourad':
        
        session_names= ['Mo180328001','Mo180330001','Mo180405001','Mo180405004','Mo180411001','Mo180412002',
                        'Mo180411001','Mo180418002','Mo180419003','Mo180503002', 'Mo180523002','Mo180524003', 
                        'Mo180525003','Mo180530004','Mo180531002','Mo180601001','Mo180614002','Mo180614006',
                        'Mo180615002','Mo180615005','Mo180426004', 'Mo180427003', 'Mo180619002','Mo180622002',
                        'Mo180626003', 'Mo180627003','Mo180629005', 'Mo180703003','Mo180704003', 'Mo180705002',
                         'Mo180706002', 'Mo180710002','Mo180711004','Mo180712006']
           
    alignData = 0
    evtMarker = 'full'
    for i_session,sess in enumerate(session_names):
        
       print(f"Procession Session{sess}")
       inFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{sess}/modifiedData" 
       outFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{sess}/FR" 
       run_compute_FR(inFolder, outFolder,alignData, None)
       
    print('Success computing FR for all sessions')
         
 