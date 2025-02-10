#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 16:10:29 2024

@author: nilanjana
"""

def run_computeTAU_FM(inFolder,outFolder,monkey,lag_limit,excludeFirstBin):
    
##modules import
    import sys,re, pickle 
    sys.path.append('/Users/nilanjana/Documents/LabWorkBench/IntrinsicTimescales/scripts/Processing_Methods')
    import pandas as pd
    from computeTAU_FontanierMethod import mainFun_computeTAU
    from open_matlabData import open_matlab_behaviour
    import os
    from scipy.io import loadmat
    import numpy as np
    from scipy.io import savemat
    import csv


    
    # Create folder if it does not exist
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
        
    df = pd.read_csv(csvList)
    df_filename = df['last_name']
    
    saveData = []
    
    
        
        # Loop through the file names
    for i, file in enumerate(df_filename):
        
        sess_name = df['last_name'][i].split('_')[0]
        #data_path = f"/Volumes/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{sess_name}/modifiedData"
        data_path = '/Users/nilanjana/Documents/LabWorkBench/IntrinsicTimescales/data/Mo180412002/modifiedData'
        
        file_name =  (''.join([file,'.mat']))
        print (f"Processing Session {file_name}")
        print (f"{sess_name}")
        
        
   # Load time series data
        file_path = os.path.join(data_path,file)
        
        data =  open_matlab_behaviour(file_path)
        
        d = data['ts']
        # Extract the last string name from the file name
        last_name = file
        
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

        # Create filename for saving output
        output_filename = os.path.join(outFolder, f"{last_name}.pkl")
        
        
         
        #Save output using pickle (or another preferred format)
        # Save output using pickle (or another preferred format)
        # with open(output_filename, "wb") as output_file:
        #     pickle.dump(out_fit, output_file)
        #     pickle.dump(out_autoc, output_file)  
        saveData.append((last_name, out_fit['TAU'][0]))
        

            
    
        # convert dataframe to array to save in mat 
        
        # arr1 = np.array(out_fit.to_records(index=False));
        # arr2 = np.array(out_autoc.to_records(index=False));
        
        # data_dic=  {"autoc" : arr2, "fit_val": arr1}
        # filename2 = os.path.join(outFolder, f"{last_name}.mat")
        # savemat(filename2, data_dic)
    # Write the data to a CSV file
    output_file = os.path.join(outFolder,"combined_output.csv")
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Filename', 'tau'])  # Write header
        csvwriter.writerows(saveData)  # Write data
            
        
        
if __name__=='__main__':
   
    #monkey = 'Tommy'
    monkey = 'Mourad'

    lag_limit = [1000]
    excludeFirstBin = 1 
   
   
    
    for l in lag_limit:
            
            
            #csvList = f"//Users/nilanjana/Documents/LabWorkBench/IntrinsicTimescales/data/Mo180412002/R_computed.csv"
            csvList = '/Users/nilanjana/Documents/LabWorkBench/IntrinsicTimescales/data/Mo180412002/R_computed/combined_output_R.csv'
            outFolder = '/Users/nilanjana/Documents/LabWorkBench/IntrinsicTimescales/data/Mo180412002/Python_computed_lowess_0.05'
            #outFolder = f"/Volumes/work/comco/nandi.n/IntrinsicTimescales/test/{monkey}/TAU_FM_excludeFirstBin_lowess_0.03"
            run_computeTAU_FM(csvList,outFolder,monkey,l,excludeFirstBin)
        
            print('Success computing all fits')
                 