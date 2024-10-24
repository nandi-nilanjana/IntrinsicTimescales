#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 11:31:13 2023

@author: nilanjana
"""
"""
This function has been modified from run_computeTAU_FM such that it can  incorporate the following changes :
 
    
 
    
alignData = 0 or  1 ; 1 means the trials will be aligned from the entire session based on a trial start tine and trial end time given by the user

t_st and t_end = these are start and end times to cut a trial. This information comes from the data stored in the mat files 
while storing in structs. The second and the third keys of data opened here using open_matlab_behaviour will always have the information

there is also an argument called trial_info which which carries the information about the trial type


This function uses the subfunction mainFun_computeTAU_TrialWise defined in folder Processing_Methods 
The entire function is the same as mainFun_computeTAU except that it takes arguments for aligning data, 
start and end times for aligning data and trial type information


At the end of this script you need to input the path variables for the experimental sessions

open_matlab_behaviour is a subfunction also called inside this function to read the mat data from the monkeys 
and convert them to arrays to be read in python



"""

"""
Before running the code always check the t_start and t_end point for cutting the trial and change the outFolder name
Also change whether you want to exclude some trials or not ; trial_info takes the value of trial numbers 1, 2 or 3
if you want to exclude blue , green or pink trials consecutively for your analysis . else put trial_info as 0 if you 
want to include all trials 
 

"""    

##modules import
import sys,re, pickle
sys.path.append('/Users/nilanjana/Documents/LabWorkBench/IntrinsicTimescales/scripts/Processing_Methods') 
from computeTAU_FontanierMethod_TrialWise import mainFun_computeTAU
from open_matlabData import open_matlab_behaviour
import os
import numpy as np
from scipy.io import savemat
    
def run_computeTAU_FM(inFolder,outFolder,lag_limit,alignData,eventMarker,excludeFirstBin):
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
        l = list(data.keys())
       
        if eventMarker == 'preSC2':
            trial_info = data[l[2]]
        else:
            trial_info = data[l[3]] #keys : 3 has trial type : 1 - blue ; 2 - green; 3 - pink 
            
        d = data[l[0]] #entire spike data from the whole session ; keys - 0 : whole session data
        
        if eventMarker  == 'preSC1' :
            #depends on what is stored in the mat data file 
            t_st= data[l[1]]-2000 #start time to cut the trial; keys -  1 - Sc1 and 2 - Go_time
            t_end= data[l[1]] #stop time to cut the trial 
            trial_exclude = 0 
        
        elif eventMarker == 'preGo':
            #depends on what is stored in the mat data file 
            t_st= data[l[2]]-2000 #start time to cut the trial; keys -  1 - Sc1 and 2 - Go_time
            t_end= data[l[2]] #stop time to cut the trial 
            trial_exclude = 3 
        elif eventMarker == 'preSC2':
            t_st= data[l[1]]-1000
            t_end= data[l[1]]
            trial_include = 1
            trial_exclude = trial_include #doing this here to to satisfy the code insidde the main function to only inlcude the blue or green trials when you extract the tau 
            
        
        
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
        
     
        if alignData == 1:
            # Run main function
            
            
            [out_fit, out_autoc] = mainFun_computeTAU(d,unit,lag_limit,alignData,trial_exclude,t_st,t_end,trial_info,excludeFirstBin)  
            
        else:
                
            # Run main function
            [out_fit, out_autoc] = mainFun_computeTAU(d,unit,lag_limit,0,'None','None','None','None',excludeFirstBin)  

        # Create filename for saving output
        output_filename = os.path.join(outFolder, f"{last_name}.pkl")
         
        #Save output using pickle (or another preferred format)
        # Save output using pickle (or another preferred format)
        with open(output_filename, "wb") as output_file:
            pickle.dump(out_fit, output_file)
            pickle.dump(out_autoc, output_file)   
            
    
        # convert dataframe to array to save in mat 
        
        arr1 = np.array(out_fit.to_records(index=False));
        arr2 = np.array(out_autoc.to_records(index=False));
        
        data_dic=  {"autoc" : arr2, "fit_val": arr1}
        filename2 = os.path.join(outFolder, f"{last_name}.mat")
        savemat(filename2, data_dic)
        
    
     
if __name__=='__main__':
    
    #path initialise

    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
    
    
    monkey = 'Mourad'             
    lag_limit =[1000]
    alignData = 1 #align the trials on eventMarker
    eventMarker = 'preGo' #can be preSC1 or preGo depending on the bhv data extracted from matlab 
    #eventMarker = 'preGo'
    
    excludeFirstBin = 1;
   
    if excludeFirstBin ==1:
        binExt = 'firstBinExcluded'
    elif excludeFirstBin ==0 :
        binExt = ''
    
    
    # if monkey == 'Tommy':
    
    #     session_names =['t140924003','t140924006','t140925001','t140926002','t140930001','t141001001','t141008001','t141010003','t141119001','t150122001',
    #                     't150123001','t150128001','t150129002','t150204001','t150205004','t150210001','t150211003','t150212001','t150218001','t150219003',
    #                     't150303002','t150319003','t150324002','t150327002','t150327003','t150415002','t150416002','t150423002','t150430002','t150520003',
    #                     't150716001']
    
    # elif monkey == 'Mourad':
    #     session_names= ['Mo180328001','Mo180330001','Mo180405001','Mo180405004','Mo180411001','Mo180412002','Mo180411001',
    #                     'Mo180418002','Mo180419003','Mo180503002', 'Mo180523002', 'Mo180524003', 'Mo180525003','Mo180530004',
    #                     'Mo180531002','Mo180601001','Mo180614002','Mo180614006','Mo180615002','Mo180615005', 'Mo180426004', 
    #                     'Mo180427003', 'Mo180619002','Mo180622002','Mo180626003','Mo180627003','Mo180629005', 'Mo180703003',
    #                     'Mo180704003', 'Mo180705002','Mo180706002', 'Mo180710002','Mo180711004','Mo180712006']
        
        #sessions I sorted later
        
        # session_names =['Mo180426004', 'Mo180427003', 'Mo180619002','Mo180622002','Mo180626003', 
        #                 'Mo180627003','Mo180629005', 'Mo180703003','Mo180704003', 'Mo180705002',
        #                 'Mo180706002', 'Mo180710002','Mo180711004','Mo180712006']
        
        
    
    #s ='t140918002' ; need to add SC1_Go in the mat file
    
    #session_names = ['Mo180412002']
    m_names = ['Mourad','Tomy']
    
    for monkey in m_names:
        for l in lag_limit:
            if monkey == 'Tommy':
            
                session_names =['t140924003','t140924006','t140925001','t140926002','t140930001','t141001001','t141008001','t141010003','t141119001','t150122001',
                                't150123001','t150128001','t150129002','t150204001','t150205004','t150210001','t150211003','t150212001','t150218001','t150219003',
                                't150303002','t150319003','t150324002','t150327002','t150327003','t150415002','t150416002','t150423002','t150430002','t150520003',
                                't150716001']
            
            elif monkey == 'Mourad':
                session_names= ['Mo180328001','Mo180330001','Mo180405001','Mo180405004','Mo180411001','Mo180412002','Mo180411001',
                                'Mo180418002','Mo180419003','Mo180503002', 'Mo180523002', 'Mo180524003', 'Mo180525003','Mo180530004',
                                'Mo180531002','Mo180601001','Mo180614002','Mo180614006','Mo180615002','Mo180615005', 'Mo180426004', 
                  'Mo180427003', 'Mo180619002','Mo180622002','Mo180626003','Mo180627003','Mo180629005', 'Mo180703003',
                                'Mo180704003', 'Mo180705002','Mo180706002', 'Mo180710002','Mo180711004','Mo180712006']
                    
            
            
            for i_session in session_names:
                session = i_session   
                inFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{session}/SC1_Go"
                outFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{session}/TAU_{eventMarker}_{binExt}_3.33"
                #bin windth is 3.33 ms 
                run_computeTAU_FM(inFolder,outFolder,l,alignData,eventMarker,excludeFirstBin)
            
            print('Success computing all fits')
     
   
        
 