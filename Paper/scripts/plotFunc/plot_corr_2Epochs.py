#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 19:09:40 2024
This script is used to check the correlation of TAU values in the 2 epochs - currently using preSC1 and preGo epochs

Modified on Sept 17, 2024 to separately check the PMd / M1 sessions
@author: nilanjana
"""

#import files from pre_SC1 and pre_Go directories , then scatter them and check if they have a correlation
#%%
#modules import
import os, pickle
import matplotlib.pyplot as plt 
import numpy as np
from scipy.stats import pearsonr, linregress, kruskal, skew
import seaborn as sns
import pandas as pd
#%% sufunction to use

def loadPickleFile (filename):
    
    with open(filename,'rb') as file:
        fit_data =  pickle.load(file)
        autoc =  pickle.load(file)
    return fit_data, autoc


def loadTAU (fit_val) :
    
        
        if fit_val['bestfit'][0]=='mono':
            
            mono_fit = fit_val[fit_val['fit_type'] == 'mono'] #find the mono fit        
            #print('Inside mono loop')
            t = mono_fit['TAU'].values
           # t = np.round(t,decimals=2)
            se = np.round(mono_fit['TAU_sd'].values,2)
            mid  = mono_fit['peak1_lat'].values
            A = mono_fit['A'].values
            B = mono_fit['B'].values
        
        else:
            t = np.nan
            se = np.nan
            mid = np.nan
            A = np.nan
            B = np.nan
           
        return t,se,mid,A,B
    
    
def test_skewness(t,a,b,autoc_data,lat):
    
    if ~ np.isnan(t):
    
        idx= autoc_data['mids'] >= lat[0]
        x = autoc_data[idx] # filtering the data from first peak to end
        y = a * np.exp(-x['mids'] / t) + b
        sk  = skew (y)
    else:
        sk = np.nan
  
    return sk   


def plot_allSess(server,session_names,probe_num,area,ax,include_all_neurons,combine_m):
    
    '''
    Input Params
    server  - the current working directory
    session_names - list of all sessions for a given monkey and probe
    probe_num - list of probe numbers (1 or 2) corresponding to PMd/M1 in a given session folder
    area = str - PMd or M1
    ax - plot axis
    include_all_neurons - 0 or 1 - for 0 , all neurons are plotted irrespective of std or skewness- they are just masked in red
                            1 means neurons are plotted after removing few based on std or skeeness
                            
    Output - figure map with corr plot in two monkeys in PMd/M1
    '''

    #defining plot colors    
    if area == 'PMd':        
        color_area = 'darkblue' 
    else:
        color_area = 'purple'
    #definig empty lists to concatenate later
    tau_preGo = [] 
    se_preGo = []
    tau_preSC1 =[] 
    se_preSC1 = []
    sk_preGo = []
    sk_preSC1 = []
        
    if combine_m == 0:      
        if session_names.iloc[0].startswith('Mo'):
            m_name = 'Mourad'
        else:
            m_name = 'Tommy'
            
        title_m =  m_name
    else:
        
        title_m = 'combined'
    
    for i, i_session in enumerate(session_names):
        
        if i_session.startswith('Mo'):
             m_name = 'Mourad'
        else:
            m_name = 'Tommy'
        
    
        inFolder_TAU_preSC1 = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{m_name}/{i_session}/TAU_preSC1_firstBinExcluded_3.33"
        inFolder_TAU_preGo = f"{server}/work/comco/nandi.n/IntrinsicTimescales/data/{m_name}/{i_session}/TAU_preGo_firstBinExcluded_3.33"

        if os.path.exists(inFolder_TAU_preGo):
            #search for filenames in the preGo session folder where the probe number should belong to the 'area' (PMd/M1)
            filenames_preGo= []
            #print(f'{i_session}')
            for filename in os.listdir(inFolder_TAU_preGo):
                if (filename.endswith('.pkl')) and (f'probe{probe_num.iloc[i]}' in filename):                
                   filenames_preGo.append(filename)
            
            for j, filename in enumerate(filenames_preGo):
                
                #print (f"Processing Session {filename}")
                path_preGo = os.path.join(inFolder_TAU_preGo,filename)
                path_preSC1 = os.path.join(inFolder_TAU_preSC1,filename)
                
                fit_preGo, autoc_preGo = loadPickleFile(path_preGo)
                fit_preSC1, autoc_preSC1 =  loadPickleFile(path_preSC1)
     #concatenated all tau, sd, and skewness from that session folder for both epochs for a single session     
                tau1,se1,lat1,A1,B1= loadTAU(fit_preGo)
                tau2,se2,lat2,A2,B2 = loadTAU(fit_preSC1)
                
                skpreGo= test_skewness(tau1, A1, B1, autoc_preGo, lat1)
                skpreSC1 =  test_skewness(tau2, A2, B2, autoc_preSC1, lat2)
                
                tau_preGo.append(tau1)
                tau_preSC1.append(tau2)
                
                se_preGo.append(se1)
                se_preSC1.append(se2)
                
                sk_preGo.append(skpreGo)
                sk_preSC1.append(skpreSC1)
        else:
            print(f'Session {i_session} does not exist for 2 epochs')
            
 #Now concatenate across all sessions      
    t_preGo = np.hstack(tau_preGo) #works like np.concatenate but also takes into account the nan values. with np.concatenate the array will give nan
    t_preSC1 =  np.hstack(tau_preSC1)

    std_preGo = np.hstack(se_preGo)
    std_preSC1 = np.hstack(se_preSC1)
 
# to remove the nan values if its present for a datapoint in either of the arrays and 
#aso if the skewness of the fitted values less than 5 in either of the arrays        
    valid_indices = ~np.isnan(t_preGo) & ~np.isnan(t_preSC1)
    x_valid = t_preGo[valid_indices]  #x_valid is preGo 
    y_valid = t_preSC1[valid_indices] #y_valid is preSC1
    
    std_err_x = std_preGo[valid_indices]
    std_err_y = std_preSC1[valid_indices]
    
    sk_x_valid = np.array(sk_preGo)[valid_indices]
    sk_y_valid = np.array(sk_preSC1)[valid_indices]
    
       # Set a threshold for standard error
    fiterr_threshold = 100
    sk_threshold = 5
   
    if not np.isnan(x_valid).all():
        
        if include_all_neurons ==1 :
        
            # Create a combined mask based on the condition
            mask = (std_err_x > fiterr_threshold) | (std_err_y > fiterr_threshold) | (sk_x_valid>sk_threshold) | (sk_y_valid>sk_threshold)
            #plot the scatter plot to show the correlation
             
            ax.scatter(x_valid,y_valid,color = np.where(mask,'red',f'{color_area}'), label=f'{len(x_valid)}')    
            #red for those datapoint whose fit error for estimating tau exceeds 100 either in preSc1 or preGo condition
            
        else:
        
            goodNeurons_indices  =  (std_err_x < fiterr_threshold) & (std_err_y < fiterr_threshold) & (sk_x_valid<sk_threshold) & (sk_y_valid<sk_threshold)
            
            x_good = x_valid[goodNeurons_indices]
            y_good = y_valid[goodNeurons_indices]
            
            x_valid = x_good
            y_valid = y_good
            ax.scatter(x_valid,y_valid,color = color_area, alpha = 0.8, label = f'{len(x_valid)}')
            
            
        if len(x_valid)>2:
            c, p_value = pearsonr(x_valid, y_valid)
               
            alpha = 0.05
            is_significant = p_value < alpha
              
            print("Correlation Coefficient:", c)
            print("P-value:", p_value)
            print("Is the correlation significant?", is_significant)
            
            s,pval = kruskal(x_valid,y_valid)
            significance = pval < alpha
            
            print(f"Comparing Median  between preSC1 and PreGo in monkey {title_m} in {area}")
            print("P-value:", pval)
            print("Is the difference significant?", significance)
            
            ax.set_xlabel(f'preGo -- Median_Tau = {np.median(x_valid):.2f}', fontsize = 18)
            ax.set_ylabel(f'preSC1 -- Median_Tau = {np.median(y_valid):.2f}', fontsize = 18)
            ax.set_title(f'{title_m}-{area} -- corr:{c:.2f} p:{p_value:.2f}', fontsize = 18)
                    
            # Fit a linear regression through origin
            slope_origin = np.mean(y_valid) / np.mean(x_valid) 
            y_range_origin = slope_origin * x_valid
            ax.plot(x_valid, y_range_origin, color='green', linestyle='--', label='Regression Line')
            ax.tick_params(axis = 'both', labelsize = 18)
        ax.legend(loc = 'upper right', fontsize = 20)  
       
        # In case you want to find the regression line that fits the data instead of passing only through origin
        
        #     slope, intercept, _, _, _ = linregress(x_valid, y_valid)
        #     slope_origin = np.mean(y_valid) / np.mean(x_valid)  # Slope for the line passing through the origin
            
        #     # Plot the regression line
        #     x_range = np.linspace(0, max(x_valid), 100)
        #     y_range = slope * x_range + intercept
        #     y_range_origin = slope_origin * x_range
          
        # Display the legend
       
        
    else:
        print(f'No neurons left for comparison between 2 epochs in this {area} in monkey {m_name}')
        
        
       
            
#%%

if __name__=='__main__' :    
    #path initialise
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
        
    monkey_names = ['Mourad', 'Tommy']
    probe_area = ['PMd', 'M1']
    
    combine_m = 1
    
    fig, axs = plt.subplots(2,2, figsize=(12,8), sharex=True, sharey=True)
    include_all_neurons = 0

    
    for p, area in enumerate(probe_area): #Pmd or M1
        df_sess_c =  []
        df_probeNum_c = []
        for m, monkey in enumerate(monkey_names):
            #read the csv file for AP and area           
           
            df = pd.read_excel(f'{server}/work/comco/nandi.n/IntrinsicTimescales/docs/{monkey}LaminarInfo.xlsx')                                 
            df_sess= df[df['Probe_Area'] == area]['Session'] #subset of sessions with specific probe area
            df_probeNum= df[df['Probe_Area'] == area]['Probe_Num'] #probe number for the corresponding area - PMd/M1
            ax = axs[m,p]            
            df_probeNum = df_probeNum.astype('int')
            
            if combine_m ==1:
                #combine the df_sess, df_probeNum
                df_sess_c.append(df_sess)
                df_probeNum_c.append(df_probeNum)
                #plot_allSess(server,np.concatenate(df_sess_c),np.concatenate(df_probeNum_c),area,ax,include_all_neurons,combine_m)
                 
            else:
                # pass the df and plot                
                plot_allSess(server,df_sess,df_probeNum,area,ax,include_all_neurons,combine_m)
                fig.suptitle('Comparison of TAU in 2 preSC1 and preGo')
                plt.tight_layout()
                sns.despine()
        
        if combine_m ==1:        
            plot_allSess(server,pd.concat(df_sess_c),pd.concat(df_probeNum_c),area,ax,include_all_neurons,combine_m)
            fig.suptitle('Comparison of TAU in 2 preSC1 and preGo')
            plt.tight_layout()
            sns.despine()     
            

    #save the figure 
    outFold = f'{server}/work/comco/nandi.n/IntrinsicTimescales/results/bothMonkeys/TAU_Compare'
    if not os.path.exists(outFold):
        os.makedirs(outFold)
        
    f_title = f'Comparison_TAU_preSC1_preGo_allNeurons_1_combinedM_{combine_m}'
    plt.savefig(os.path.join(outFold,f_title))
            
        
    

# if (monkey == 'Mourad') & (area == 'PMd'):
    
#     session_names = ''
    
# elif  (monkey == 'Mourad') & (area == 'M1'):
   
#     session_names = ''
    
# elif (monkey == 'Tomy') & (area == 'PMd'): 
#     session_names = ''

# elif (monkey == 'Tommy') & (area == 'M1'):
#     session_names =  ''
    




# plotData = 1
# include_all_neurons = 0



    

# tau_preGo = [] 
# se_preGo = []
# tau_preSC1 =[] 
# se_preSC1 = []
# sk_preGo = []
# sk_preSC1 = []



# for i_session in session_names:
#     session = i_session
    
#     inFolder_TAU_preSC1 = f"/Volumes/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{session}/TAU_preSC1_firstBinExcluded_3.33"
#     inFolder_TAU_preGo = f"/Volumes/work/comco/nandi.n/IntrinsicTimescales/data/{monkey}/{session}/TAU_preGo_firstBinExcluded_3.33"

#     #search for filenames in the folder where tau values estimated from the trial before preGo 
#     filenames_preGo = [filename for filename in os.listdir(inFolder_TAU_preGo) if filename.endswith(".pkl")]
    
#     for i, filename in enumerate(filenames_preGo):
        
#         print (f"Processing Session {filename}")
#         path_preGo = os.path.join(inFolder_TAU_preGo,filename)
#         path_preSC1 = os.path.join(inFolder_TAU_preSC1,filename)
        
#         fit_preGo, autoc_preGo = loadPickleFile(path_preGo)
#         fit_preSC1, autoc_preSC1 =  loadPickleFile(path_preSC1)
        
#         tau1,se1,lat1,A1,B1= loadTAU(fit_preGo)
#         tau2,se2,lat2,A2,B2 = loadTAU(fit_preSC1)
        
#         skpreGo= test_skewness(tau1, A1, B1, autoc_preGo, lat1)
#         skpreSC1 =  test_skewness(tau2, A2, B2, autoc_preSC1, lat2)
        
#         tau_preGo.append(tau1)
#         tau_preSC1.append(tau2)
        
#         se_preGo.append(se1)
#         se_preSC1.append(se2)
        
#         sk_preGo.append(skpreGo)
#         sk_preSC1.append(skpreSC1)     
        
        
        
        
# if plotData ==1:
#     t_preGo = np.hstack(tau_preGo) #works like np.concatenate but also takes into account the nan values. with np.concatenate the array will give nan
#     t_preSC1 =  np.hstack(tau_preSC1)

#     std_preGo = np.hstack(se_preGo)
#     std_preSC1 = np.hstack(se_preSC1)
    
#     valid_indices = ~np.isnan(t_preGo) & ~np.isnan(t_preSC1) # to remove the nan values if its present for a datapoint in either of the arrays and aso if teh skewness of the fitted values less than 5 in either of the arrays
    
#     x_valid = t_preGo[valid_indices]  #x_valid is preGo 
#     y_valid = t_preSC1[valid_indices] #y_valid is preSC1
    
#     std_err_x = std_preGo[valid_indices]
#     std_err_y = std_preSC1[valid_indices]
    
#     sk_x_valid = np.array(sk_preGo)[valid_indices]
#     sk_y_valid = np.array(sk_preSC1)[valid_indices]
    
#        # Set a threshold for standard error
#     fiterr_threshold = 100
#     sk_threshold = 5
   
#     if include_all_neurons ==1 :
    
#         # Create a combined mask based on the condition
#         mask = (std_err_x > fiterr_threshold) | (std_err_y > fiterr_threshold) | (sk_x_valid>sk_threshold) | (sk_y_valid>sk_threshold)
#         #plot the scatter plot to show the correlation
         
#         plt.scatter(x_valid,y_valid,color = np.where(mask,'red','blue'))    #red for those datapoint whose fit error for estimating tau exceeds 100 either in preSc1 or preGo condition
        
#     else:
    
#         goodNeurons_indices  =  (std_err_x < fiterr_threshold) & (std_err_y < fiterr_threshold) & (sk_x_valid<sk_threshold) & (sk_y_valid<sk_threshold)
        
#         x_good = x_valid[goodNeurons_indices]
#         y_good = y_valid[goodNeurons_indices]
        
#         x_valid = x_good
#         y_valid = y_good
#         plt.scatter(x_valid,y_valid,color = color_plot, alpha = 0.8)


#     c, p_value = pearsonr(x_valid, y_valid)
       
#     alpha = 0.05
#     is_significant = p_value < alpha
      
#     print("Correlation Coefficient:", c)
#     print("P-value:", p_value)
#     print("Is the correlation significant?", is_significant)
    
#     s,pval = kruskal(x_valid,y_valid)
#     significance = pval < alpha
    
#     print("Comparing Median  between preSc1 and PreGo")
#     print("P-value:", pval)
#     print("Is the difference significant?", significance)
    
     
#     plt.title(f"TAU from {monkey} -- n = {len(x_valid)}")
#     plt.xlabel(f"TAU preGo -- Median :{round(np.median(x_valid),2)} ms ")
#     plt.ylabel(f"TAUpreSC1 -- Median :{round(np.median(y_valid),2)} ms" )
    
#     # Annotate the plot with the correlation coefficient
#     #annotation_text = f'Correlation: {c:.2f}\nP-value: {p_value}'
#     #plt.annotate(annotation_text, xy=(0.5, 0.95), xycoords='axes fraction', ha='center', fontsize=10, color='grey')
     
    

    
    
#     plt.plot(x_range, y_range, color=color_plot , alpha = 0.5, label='Regression Line')
#     plt.plot(x_range, y_range_origin, color='green', linestyle='--', label='Regression Line (Through Origin)')
      
#     # Display the legend
#     plt.legend()    
#     plt.show()           
    
        
#     #if you wish to compare the medians between tau pre SC1 and pre -Go
    
#     statistic, pVal = kruskal(x_valid, y_valid) 
    
    

    
    
    
    
 #%%
# import pandas as pd
# df = pd.read_csv('/Volumes/work/comco/nandi.n/IntrinsicTimescales/docs/Tommy_AreaInfo.csv')  
    
# df_PMd = df[df['Session'] == 'PMd']['Session']   
# df_M1 = df[df['Session'] == 'M1']['Session']     
    
#     session_names =['t140924003','t140924006','t140925001','t140926002','t140930001','t141001001','t141008001','t141010003','t141119001','t150122001',
#                         't150123001','t150128001','t150129002','t150204001','t150205004','t150210001','t150211003','t150212001','t150218001','t150219003',
#                         't150303002','t150319003','t150324002','t150327002','t150327003','t150415002','t150416002','t150423002','t150430002','t150520003',
#                         't150716001']
#     if monkey == 'Mourad':
        
#         session_names= ['Mo180328001','Mo180330001','Mo180405001','Mo180405004','Mo180411001','Mo180412002','Mo180411001',
#                         'Mo180418002','Mo180419003','Mo180503002', 'Mo180523002', 'Mo180524003', 'Mo180525003','Mo180530004',
#                         'Mo180531002','Mo180601001','Mo180614002','Mo180614006','Mo180615002','Mo180615005', 'Mo180426004', 
#                         'Mo180427003', 'Mo180619002','Mo180622002','Mo180626003','Mo180627003','Mo180629005', 'Mo180703003',
#                         'Mo180704003', 'Mo180705002','Mo180706002', 'Mo180710002','Mo180711004','Mo180712006']
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
