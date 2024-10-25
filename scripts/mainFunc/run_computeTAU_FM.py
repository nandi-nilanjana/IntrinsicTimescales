#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 11:31:13 2023

@author: nilanjana
"""
"""

This function takes 3 arguments - folderpath for finding the spike data as inFolder,
path to store the processed data as outFolder and 
lag_limit as an integer variable. This varaible computes the lag over which the autocorrelogram is computed. 
Since the data is always stored in the server, you need to be always connected to VPN to access and store the data


This function uses the subfunction mainFun_computeTAU defined in folder Processing_Methods 

At the end of this script you need to input the path variables for the experimental sessions

open_matlab_behaviour is a subfunction also called inside this function to read the mat data from the monkeys 
and convert them to arrays to be read in python



"""    

def run_computeTAU_FM(inFolder,outFolder,lag_limit,excludeFirstBin):

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
        
     
    
        # Run main function
        [out_fit, out_autoc] = mainFun_computeTAU(d,unit,lag_limit,excludeFirstBin)
        
        if out_fit is None or out_autoc is None:
            print(f"Skipping {file_name} due to computation issue of SVD converge")
            continue #skips this iteration and moves on to the next neuron

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
    
    import os 
    
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
        
    import sys,re, pickle 
    if server == '/envau':
        sys.path.append('/envau/work/comco/nandi.n/IntrinsicTimescales/scripts/Processing_Methods')
    else:
        sys.path.append('/Users/nilanjana/Documents/LabWorkBench/IntrinsicTimescales/scripts/Processing_Methods')
    
    from computeTAU_FontanierMethod import mainFun_computeTAU
    from open_matlabData import open_matlab_behaviour     
    # from scipy.io import loadmat
    import numpy as np
    from scipy.io import savemat    
       
   
    monkey = 'Tomy'
    #monkey = 'Mourad'

    lag_limit = [1000]
    excludeFirstBin = 1 
    
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
            
        for i_session in session_names:
            session = i_session   
            inFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{session}/modifiedData"
            outFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{session}/TAU_FM_excludeFirstBin_loess_0.1"
            run_computeTAU_FM(inFolder,outFolder,l,excludeFirstBin)
        
        print('Success computing all fits')
         
 