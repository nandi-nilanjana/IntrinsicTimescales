#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 11:47:30 2024

@author: nilanjana
"""

'''
1.This code is used to plot the probability distribution of ISIs in bins and also the exponential fit for each neuron for visualisation
For now the tau.sd is kept at 150. Any neuron with its TAU_mono fit std greater than 150 is removed. 


This code saves the fits in the folder TAU_fits_plots and also the csv file containing all the neurons 

First run the code with csvDataSave = 0 to save all the plots first. Go through each neuron
Sometimes the sd threshold is not enough to remove the bad tau fits. 

Then manually remove the neurons inside each session folder called TAU_exclude_fits

Then use this code again to extract all the fits and put them in CleanedSession folder. At this time use
csvDatasave =1 to save the neuron and the tau information inside a csv file which i can use later.

So till now i manually change the outFolder from AllSessions to CleanedSessions (after i manually remove the neurons
                                                                                 from each session folder to NotInUse)

INPUT :
    
    monkey name - Mourad / Tomy 
    csvSave = 0 or 1 (basically to save each neuron and tau information or not)
    run_pseudoPopulation = 0 or 1 (for pseudopopulation you have different folder paths )

'''



def runplotData(inFolder,outFolder,saveCSVData):
    
   
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
            
    
    
    file_names =  [filename for filename in os.listdir(inFolder) if filename.endswith(".pkl")]
    
    tau_save= []
    fl_nametoSave =[]
    for i, file_name in enumerate(file_names):
        
        print(f'Processing session {file_name}')
        last_name = os.path.splitext(os.path.basename(file_name))[0]
        
        data_path = os.path.join(inFolder, file_name)
        
        with open(data_path,'rb') as file:
            params = pickle.load(file)
            autoc = pickle.load(file)
   
            #fig, axs = plt.subplots(1, 2, figsize=(7, 3), sharex=True)
            fig, axs = plt.subplots(1, 1)
            x = autoc['mids']
            y = autoc['countsn']
            
            fig.suptitle(last_name)
            
            if len(np.shape(axs))>0:
                axs[0].plot(x,y)
                axs[0].set_xlabel('lag(ms)')
                axs[0].set_ylabel('data without smoothen')
                
                x2 = autoc['mids']
                y2 =  autoc['countsnsm']
                axs[1].plot(x2,y2)
            
            else :
                
                axs.plot(x,y,alpha = 0.5)
                axs.set_xlabel('lag(ms)')
               # axs.set_ylabel('data without smoothen')
                
                x2 = autoc['mids']
                y2 =  autoc['countsnsm']                
                axs.plot(x2,y2)
                
                l_pre500 = autoc['mids']<=500
                l_post500 = autoc['mids']>=500
                
                s1= autoc['counts'][l_pre500].sum()
                s2 = autoc['counts'][l_post500].sum()
                
                axs.text(5, y.max(), f'pre500_ISI_count = {s1}', fontsize=12, ha='center')
                axs.text(800, y.max(), f'post500_ISI_count = {s2}', fontsize=12, ha='center')
                
                
                
           
            # Generate the fit using the obtained parameters for a, b, and TAU if mono fit exists
            if params['bestfit'][0] == 'mono': #checks if the best fit is mono or not

                mono_fit = params[params['fit_type'] == 'mono'] #find the mono fit
        
                if mono_fit['TAU_sd'].values <150:
             
                           idx= autoc['mids'] >= params['peak1_lat'][0]  
                           data2 = autoc[idx] # filtering the data from first peak to end
                           data2 = data2.reset_index(drop=True)
                           x3 = data2['mids']
                           t = mono_fit['TAU'].values
                           a = mono_fit['A'].values
                           b = mono_fit['B'].values                  
                           t_sd =  mono_fit['TAU_sd'].values
                           
                           if saveCSVData ==1:
                               tau_save.append(np.round(t,2))
                               fl_nametoSave.append(last_name)
                               
                         # Calculate the fitted curve
                           
                           y_fit = a * np.exp(-x3 / t) + b
                           
                          # Calculate the upper and lower bounds using standard deviation of the estimates
                           y_upper = (a + mono_fit['A_sd'].values) * np.exp(-x3 / (t+mono_fit['TAU_sd'].values)) + (b+mono_fit['B_sd'].values)         
                           y_lower = (a - mono_fit['A_sd'].values) * np.exp(-x3 / (t-mono_fit['TAU_sd']).values) + (b-mono_fit['B_sd'].values)
                    
                           if len(np.shape(axs))>0:
                            # Plot the exp fit curve with confidence interval
                               axs[1].plot(x3, y_fit, label=f"TAU: {np.round(t[0],decimals=2)}", color = 'red') 
                               
                               axs[1].fill_between(x3, y_lower, y_upper, alpha=0.4,color = 'red')  
                           else:
                               
                               axs.plot(x3, y_fit, label=f"TAU: {np.round(t[0],decimals=2)}", color = 'red') 
                               
                               axs.fill_between(x3, y_lower, y_upper, alpha=0.4,color = 'red')  
            fig.tight_layout
            fig.savefig(os.path.join(outFolder,f'{last_name}.png'))
            plt.close(fig)
            
    return tau_save, fl_nametoSave
        
 

if __name__=='__main__':

    #import modules
    import pickle
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import csv
    import os
    
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
    
    monkey_names = ['Tomy', 'Mourad']
    #monkey = 'Mourad'
    run_pseudoPopulation = 1
    csvDataSave =0

    lag_limit = [1000]
    excludeFirstBin = 1 
    
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
                            'Mo180626003', 'Mo180627003','Mo180629005', 'Mo180703003','Mo180704003', 'Mo180705002',
                            'Mo180706002', 'Mo180710002','Mo180711004','Mo180712006']
        
        savetau1 = []
        savefile = [] 
        headers = ["Filename", "TAU"]
          
        
        for session in session_names:
            
            if run_pseudoPopulation == 0:   
                
                inFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{session}/TAU_FM_excludeFirstBin_loess_0.1"
                outFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/TAU_fits_plots/{monkey}/AllSessions_TAU_fits_loess_0.1_TAUsd_lessThan150"        
                csvOutFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/docs/{monkey}_CleanedSessions_TAU_fits_loess_0.1_TAUsd_lessThan150.csv"
                
            elif run_pseudoPopulation == 1 :
                
                inFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{session}/TAU_FM_PseudoPopulation"
                outFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/TAU_fits_plots/{monkey}/PseudoPop_TAU_fits_loess_0.1_TAUsd_lessThan150"        
                csvOutFolder = f"{server}/work/comco/nandi.n/IntrinsicTimescales/docs/{monkey}_PseudoPop_TAU_fits_loess_0.1_TAUsd_lessThan150.csv"
                
            tauPersession, filename_all= runplotData(inFolder,outFolder,csvDataSave)
            
            if csvDataSave ==1:
                 savetau1.append(np.concatenate(tauPersession))
                 savefile.append((filename_all))
    
            print('Success plotting all fits')  
            
        
        if csvDataSave ==1:
            
             data_to_save =  pd.DataFrame({'Neuron' : np.concatenate(savefile), 'TAU' : np.concatenate(savetau1)})
            
             #Save DataFrame to CSV file with specified file path
             data_to_save.to_csv(csvOutFolder, index=False)
             print('CSV File saved successfully with tau values')
        
       
        
       

            
            
            
            
            
            