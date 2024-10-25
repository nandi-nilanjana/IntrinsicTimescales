#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 10:37:44 2023

@author: nilanjana


Updated on 24 May,2024
 -  If the max ISI count occurs in the first bin, it is excluded and the next bin is taken as the peak, after which the exp fit is drawn. 

Subfunction1  - compute_autoc(arg1,arg2) : takes the timeseries data over the whole session as arg1 and 
lag over which the autocorrelogram is computed as arg2

First computes the ISI differences over order 100, creates a histogram with bin with 3.3 ms, where spike counts in each bin is normalised by the total spike count, 
First 10ms is always discarded and then the distribution is smoothened

output - dataframe containig the mid point value of each bin, actual spike counts in each bin, normalized spike counts 
and smoothened counts after normalization counts in each bin. For the rest of the analysis, the smoothened normalized spike counts are used

Subfunction 2 - dip_func(arg1,arg2) - takes the output (autoc data) from subfunc1 as arg1 and 
neuron number (integer var) (>1 if more than 1 neuron was sorted in a single channel) as arg2 

output - df_peak - returned by the function - type dictionary -  a dataframe with autoc val and lat value for peak1, dip and peak2

peak1_idx = index of the max value in autoc[‘autocnsm’]
If the index is 0 (the first value of the array) then it discards it and moves on to the next index. It runs a loop till it finds a local maxima and calls it peak1. 
If peak1 ~= 0 (the first value of the array) , then it is taken as the index of the max value, that is peak1.

curr_post_peak1 - entire autoc data after the first peak detected

curr_min - global minima of the entire data after the first peak detected
curr_max - global maxima of the entire data after the first peak detected 

curr_begin =  entire data between first peak and 100 ms post first peak.
loc_min = minima in curr_begin

Next the code checks if this loc_min is a local minima or not. 
dip - If local minima condition holds true then it is considered as the dip iff the loc_min is greater than equal to 75% of AC(max,loc_min)

Given dip exists ; now the code tries to search for a second peak 12 ms after the dip .
curr_end - data filtered as -   autoc[‘mids’]>= dip2+ 12ms
Then if there is some data left it tries to search for the second peak by searching for the maxima in curr_end  else it just saves NaN values for second peak value and lat.

So finally df_peak returns a dataframe with autoc val and lat value for peak1, dip and peak2

Subfunction 3 exp_fit_fun(arg1,arg2,arg3)
used to generate the A,B, and TAU parameters for the exponential function

Input:
- Takes the autoc data as input variable 
- unit refers to the neuronal unit - has to be an integer 
- st - initial starting values for A, B and TAU for fitting the data. 

Inside the function:
    
df_fit - dataframe to store the exponential fit parameters over the autoc data : it saves the following params:
- unit refers to the neuronal unit - has to be an integer 
- A, B and TAU optimised after the fit
- TAU.se - standard error for TAU
- fit.err - residual error for the fitted values 

For the exponential function of the form : y = a * (np.exp(-x / TAU)) + b 
		y = autoc data in each bin 
		x = value of the mid point of each bin
		a , b - constants ; TAU - decay constant 
        
Output - df_fit as output
        

Subfunction 4 -  exp_fit_autoc_post_hoc(arg1,arg2) 


Input - df_fit as arg1
unit - neuron number as arg2

- Takes all the fits returned by the exp_fit_func on the autoc data
- First filters out all fits that have negative TAU, A, B values and those TAU values whose covariances could not be estimated (that is TAU.se = infiniy)
- Then the filtered data is grouped by the fit_type and only values for each group are filtered which has the min value in fit.err
- Then it checks whether each fit type exists or not in has_fittype variables
- If any of the fits are missing, then it updates the autoc_fit dataframe by adding nan values for all the the variables like A, B, TAU .. to match the same length of the dataframe and only in fit_type updates the name of the fit.
- Calculates the rmse of TAU of mono fit and sum of rmse of TAU of dip and peak2
- Then the code uses np.select to select the best fit amongst all the 3 fit types -
- The resulting best_fit array will have values of "dip" where dum['rmse_dip'] is smaller, "dip" where 'rmse_dip' is not NaN and 'rmse_mono' is NaN, "mono" where 'rmse_mono' is not NaN, and "none" for all other cases. This is a way to categorize or label the data based on these conditions. 

- If there are no fits, then everything is returned as an empty array of nan values in bestfit 



"""


#currently the function is modified such that the mono fit is getting calculated even if the first peak starts from the first mid value; after removing first 10ms and  the width of the bins is 3.3 ms

#import libraries
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
from loess.loess_1d import loess_1d

#can be used to create the function compute_autoc -- Subfunction 1
#data - takes spiketimeseries data input as a row vector ; 

def compute_autoc(data,lag_limit) :
    # Specify the maximum lag order
    diff_order = 100  # Change this to the desired maximum lag order

    # Initialize an empty list to store lagged differences
    lagged_diffs = []

    # Loop through different lag orders and calculate lagged differences
    for lg in range(1, diff_order + 1):
        lagged_diff = data[lg:] - data[:-lg]
        lagged_diff = lagged_diff[(lagged_diff>0) & (lagged_diff<lag_limit)]
        lagged_diffs.append(lagged_diff)

    # concatenates the lagged differences between ISIs of diff_order order
    diff = np.concatenate(lagged_diffs)    

    #num_bins = 300
    bin_range = (0,lag_limit)
    
    bin_width = 3.33
    num_bins = int(np.ceil((lag_limit-0) / bin_width))

    # Create a histogram
    val, bin_edges = np.histogram(diff, num_bins ,range=bin_range)
    # Calculate the midpoints of bins
    bin_mids = (bin_edges[1:] + bin_edges[:-1]) / 2

    
    data2 = {
        'mids': bin_mids,
        'counts': val,


    }

    # Create a DataFrame
    curr_autoc = pd.DataFrame(data2)

    # Compute total counts and bin width
    curr_autoc_tot = curr_autoc['counts'].sum()
    wd = np.mean(np.diff(curr_autoc['mids']))


    # Compute the PDF values
    curr_autoc['countsn'] = (curr_autoc['counts'] / curr_autoc_tot) * (1/wd)

    # Smooth the PDF using loess
    # lowess = sm.nonparametric.lowess(curr_autoc['countsn'], curr_autoc['mids'], frac=0.05)
    # countsnsm = lowess[:, 1]
    
    #if the smoothening does not happen exit the function and move on to the next neuron
    try:
    
        x,countsnsm,w = loess_1d(curr_autoc['mids'].values,curr_autoc['countsn'].values,frac = 0.1)
    
        pass
    
    except Exception as e:
        
        print(f" loess1d encountered: {e}")
        return None  # Indicate failure to converge

   
    
   #  x =  curr_autoc['mids'].values
   #  # y = curr_autoc['countsn'].values

   #  # y_sm = localreg.localreg(x,y,degree =2 , frac = 0.1)    
    
    curr_autoc['countsnsm'] =  countsnsm
   #  curr_autoc['countsnsm'] = y_sm

    # # Filter the data to remove values with mids less than 10 -- just doing it according to Fontanier code

    #Filter out values based on mids > 10 and update the DataFrame
    filtered_indices = curr_autoc['mids'] > 10

    curr_autoc = curr_autoc[filtered_indices]
    curr_autoc=curr_autoc.reset_index(drop=True)

     # Pdf again after filtering the first 10ms ; updating both countsn and countsnsm columns

    curr_autoc_tot = curr_autoc['countsnsm'].sum()
    curr_autoc['countsnsm'] = (curr_autoc['countsnsm'] / curr_autoc_tot) * (1/wd) #smoothened
    curr_autoc['countsn'] = (curr_autoc['countsn'] / curr_autoc['countsn'].sum()) * (1/wd) #without smoothening
    
    
    
    
    
    # Find the maximum value in the countsnsm column
    max_value = curr_autoc['countsnsm'].max()
    #curr_autoc['countsnsm'] =  curr_autoc['countsnsm']/max_value
   
    
    
    return curr_autoc


## Detect peak and putative dip/additional peak - Sunfunction 2
# most of the variable names and peak detection is based on Fontanier's code in R

# Peak detection
 # Peak detection
def dip_func(autoc,unit,excludeFirstBinPeak):

#autoc = curr_autoc #later I will pass this as a variable in the function
   
    # Find the index of the maximum count
    peak1_idx = autoc['countsnsm'].idxmax()
   
    #here I have commented out the condition that if the peak1 is the first point, I need to search for the next global maxima to call it a peak1. 
    
   
    if np.isnan(peak1_idx):
    
        df_peak = {'peak1_lat' : 'NA', 'peak1_val' : 'NA'}
    
    else:
        
        if excludeFirstBinPeak ==1:
            
            if peak1_idx == 0:
                peak1_idx = 1
                while (
                    peak1_idx < len(autoc['countsnsm'])-1 and
                    not (
                        autoc['countsnsm'][peak1_idx - 1] < autoc['countsnsm'][peak1_idx] and
                        autoc['countsnsm'][peak1_idx + 1] < autoc['countsnsm'][peak1_idx]
                    )
                ):  
                    peak1_idx += 1  
                    #print(peak1_idx)

            #if there is a local maxima which is the first peak; first check if it is the end index or not;   

            if peak1_idx + 1 == len(autoc['countsnsm']) : #adjusting by adding 1 since the last index is 1 less than the length of dataframe

                print('No peak found for autocorrelation analysis')
                
                df_peak = {'peak1_lat' : 'NA', 
                            'peak1_val' : 'NA',}
            else:
                
                #creating dataframe to save the mid point 1of the bin that detected the first peak(global/local) and the corresponding value
                df_peak = {'peak1_lat' : autoc.at[peak1_idx,'mids'], 
                                                'peak1_val' : autoc.at[peak1_idx,'countsnsm']}
                
                
            
        else:
            
                
            #creating dataframe to save the mid point 1of the bin that detected the first peak(global/local) and the corresponding value
            df_peak = {'peak1_lat' : autoc.at[peak1_idx,'mids'], 
                                            'peak1_val' : autoc.at[peak1_idx,'countsnsm']}
        
    #if there is global peak then there is no point searching for dip and second peak;
    #so creating the if else statemnet such that it runs the entire part below only if there is a global peak
    #if not by default all other dip and peak2 values gets defined as nan
    
        
    if not(df_peak['peak1_lat']=='NA'):
            # Putative dip detection from first peak
            filt_idx =autoc['mids'] >= df_peak['peak1_lat']
            curr_post_peak1 = autoc[filt_idx] #entire data after the first peak
    
            ## Global (post peak1) min and max
           # curr_min =np.min(curr_post_peak1['countsnsm'])
            curr_max = np.max(curr_post_peak1['countsnsm'])
    
            ## Minimum of the AC from first peak until 100ms after it
    
            # filter the data between post peak1 and 100 ms after it
            filt_data = curr_post_peak1[curr_post_peak1['mids']==(df_peak['peak1_lat'])+100]
            
            #it checks whether there is any data 100 ms after the first peak - leads to no dip
            if filt_data.empty:
                
                dip = False
                
            else:
                
                idx  = filt_data.index[0]
                curr_begin = curr_post_peak1.loc[:idx]
                curr_begin=curr_begin.reset_index(drop=True) #update the indices from 0 to end
                
            #here it checks if data even exists betweeen post peak 1 and 100 ms after it - leads to no dip
            
                if curr_begin.empty:
                    dip = False
                else: 
                    min_idx = curr_begin['countsnsm'].idxmin() #index of the local minima
                    loc_min = curr_begin['countsnsm'][min_idx]
        
              # check if is  a local minimum, (that is AC(xi-1)>AC(xi)<AC(xi+1)) 
        
        # if index is at the end of curr_begin it will throw an error - so this will also lead to no dip
            
                    if ( min_idx < len(curr_begin['countsnsm']) and
                         min_idx+1 < len(curr_begin['countsnsm']) and
                         (
                            curr_begin['countsnsm'][ min_idx - 1] > curr_begin['countsnsm'][min_idx] and
                            curr_begin['countsnsm'][min_idx + 1] > curr_begin['countsnsm'][min_idx]
                        )
                    ):  
                        dip = curr_max - loc_min >= 3/4 * (curr_max-loc_min) #only considered as the second dip if the local minima is greater than equal to 75% of AC(max,min)
                    else:
                        dip = False
    
            #second peak detection if there is a dip 
    
            if (dip== True):  
    
                df_peak['dip_lat'] = curr_begin['mids'][min_idx] # add updated df_peak with the value of local minima and the mid point of the bin often referred as lat by Fontanier
                df_peak['dip_countsnsm'] = loc_min
                
                #filter the autoc data 12ms after the dip 
                idx  = autoc['mids'] >= df_peak['dip_lat'] + 12 #according to Fontanier its 12ms
                curr_end = autoc[idx]
                curr_end = curr_end.reset_index(drop=True)
                
    
                # here also need to include criteria such that if there is no second peak after the dip then you write lat and countsnsm values as NA for peak2 
                if curr_end.empty:
                    df_peak['peak2_lat'] = 'NA'
                    df_peak['peak2_countsnsm'] = 'NA'
                else:
                    ind_max = curr_begin['countsnsm'].idxmax() #index with the maximum value/second peak in curr_begin
                    df_peak['peak2_lat'] = curr_end['mids'][ind_max]
                    df_peak['peak2_countsnsm'] = curr_end['countsnsm'][ind_max]
    
    
            else:
                
                df_peak['dip_lat'] = 'NA'
                df_peak['dip_countsnsm'] = 'NA'
                df_peak['peak2_lat'] = 'NA'
                df_peak['peak2_countsnsm'] = 'NA'
    else:
                
            df_peak['dip_lat'] = 'NA'
            df_peak['dip_countsnsm'] = 'NA'
            df_peak['peak2_lat'] = 'NA'
            df_peak['peak2_countsnsm'] = 'NA'
            

    return df_peak



### Define starting values for the fits --- Subfunction 3
  # 50 Random selection in selected range, we keep the same starting values for all the units
  # A [0 ; 2*range(Autoc)]
  # B [0 ; 2*min(Autoc)]
  # TAU [0 ; 1000]
def start_val(autoc):
    TAU = np.random.uniform(low=0, high=1000, size=50)
    A = np.random.uniform(low=0, high=2 * (max(autoc['countsnsm']) - min(autoc['countsnsm'])), size=50)
    B = np.random.uniform(low=0, high=2 * min(autoc['countsnsm']), size=50)
    
    st = pd.DataFrame({'A': A, 'B': B,'TAU': TAU, })
    
    return st


#defining the expfit function -- subfunction 4
def expfit(x, a, b, TAU):
    return a * (np.exp(-x / TAU)) + b


#subfunction used to fit the autoc data using expfit function -- sunfunction 5
def exp_fit_fun(autoc, unit, st):
    df_fit = pd.DataFrame()
    try:
        #pdb.set_trace()
        
        popt, pcov = curve_fit(expfit, autoc['mids'].values, autoc['countsnsm'].values, p0=np.array(st)) ##changed just now
        #estimates
        A = popt[0]
        B = popt[1]
        TAU = popt[2]
        
        #standard deviation of the estimates
        A_sd = np.sqrt(np.diag(pcov))[0]
        B_sd = np.sqrt(np.diag(pcov))[1]
        TAU_sd = np.sqrt(np.diag(pcov))[2]
        
       
        fitval = expfit(autoc['mids'],A, B, TAU) #fitted counstssm values based on optimised A,B and TAU params
        fit_err = np.sqrt(np.sum((fitval - autoc['countsnsm']) ** 2))
        
        df_fit = pd.DataFrame({
            'unit': [unit],
            'TAU': [TAU],
            'TAU_sd': [TAU_sd],
            'A': [A],
            'A_sd' : [A_sd],
            'B_sd' : [B_sd],
            'B': [B],
            'fit_err': [fit_err]
        }
        )
    except Exception as e:
        df_fit = pd.DataFrame({
            'unit': [unit],
            'TAU': [np.nan],
            'TAU_sd': [np.nan],
            'A': [np.nan],
            'A_sd': [np.nan],
            'B': [np.nan],
            'B_sd': [np.nan],
            'fit_err': [np.nan]
        })
    
    return df_fit


#curr_autoc as the input variable for autoc -- Subfunction 6
#function to obtain the fits on the autoc data; from peak1 to end and if existent from second peak to end; 
#it runs the subfunctions 3-5 within this

def exp_fit_autoc(autoc,unit,excludeFirstBin):  

    d = dip_func(autoc,unit,excludeFirstBin) 
    
    # Fit from the 1st peak till the end if there is a global maxima apart from the first value in compute autoc data
    
    if not(d['peak1_val']=='NA'):
        filt_idx =autoc['mids'] >= d['peak1_lat']  # Fit from the 1st peak till the end
        curr_mono = autoc[filt_idx] #entire data after the first peak
        curr_mono = curr_mono.reset_index(drop = True)
        
        # Repeat the fit 50 times with random starting values
        st  = start_val(curr_mono)
        curr_mono_fit =  pd.DataFrame()
        
       
        for i, row in st.iterrows():
           
            A = row['A']
            B = row['B']
            TAU = row['TAU']
  
            r = exp_fit_fun(curr_mono, unit, [A,B,TAU])
            curr_mono_fit = pd.concat([curr_mono_fit,r], ignore_index=True)
    
        curr_mono_fit['dof']= curr_mono.shape[0] #same as len(curr_mono)
        curr_mono_fit['fit_type'] = 'mono'
    
    else:
    
        curr_mono_fit = pd.DataFrame({
            'unit': [unit],
            'TAU': [np.nan],
            'TAU_sd': [np.nan],
            'A_sd': [np.nan],
            'B_sd': [np.nan],
            'A': [np.nan],
            'B': [np.nan],
            'fit_err': [np.nan],
            'dof': [np.nan],
            'fit_type': ["mono"],
        })
        
    curr_fit = curr_mono_fit
   # pdb.set_trace()
  
 # Fit from 1st peak to dip if there is a dip after the first peak

    if not(d['dip_countsnsm'] == 'NA'): 
        curr_first_peak_to_dip = autoc[
            (autoc['mids'] >= d['peak1_lat']) &
            (autoc['mids'] <= d['dip_lat'])
        ]
        curr_first_peak_to_dip = curr_first_peak_to_dip.reset_index(drop= True)
        
        st2 = start_val(curr_first_peak_to_dip)
        curr_first_peak_to_dip_fit = pd.DataFrame()
        
        for index, row in st2.iterrows():
            
            r2 = exp_fit_fun(curr_first_peak_to_dip, unit, [row['A'], row['B'], row['TAU']])
            
            curr_first_peak_to_dip_fit = pd.concat([curr_first_peak_to_dip_fit, r2], ignore_index=True)
        
        curr_first_peak_to_dip_fit['dof'] = len(curr_first_peak_to_dip)
        curr_first_peak_to_dip_fit['fit_type'] = "dip1"
    else:
        curr_first_peak_to_dip_fit = pd.DataFrame({
            'unit': [unit],
            'TAU': [np.nan],
            'TAU_sd': [np.nan],
            'A_sd': [np.nan],
            'B_sd': [np.nan],
            'A': [np.nan],
            'B': [np.nan],
            'fit_err': [np.nan],
            'dof': [np.nan],
            'fit_type': ["dip1"],
        })

        
    curr_fit = pd.concat([curr_fit,curr_first_peak_to_dip_fit],ignore_index= True)
        
 ## Fit second peak to the end if there is a second peak

    if not(d['peak2_countsnsm'] == 'NA'):
        curr_second_peak_to_end = autoc[autoc['mids'] >= d['peak2_lat']]
        st3 = start_val(curr_second_peak_to_end)
        curr_second_peak_to_end_fit = pd.DataFrame()
        
        for index, row in st3.iterrows():

            r3 = exp_fit_fun(curr_second_peak_to_end, unit, [row['A'], row['B'], row['TAU']])
            curr_second_peak_to_end_fit = pd.concat([curr_second_peak_to_end_fit, r3], ignore_index=True)

        curr_second_peak_to_end_fit['dof'] = len(curr_second_peak_to_end)
        curr_second_peak_to_end_fit['fit_type'] = "dip2"
    else:
        
        curr_second_peak_to_end_fit = pd.DataFrame({
            'unit': [unit],
            'TAU': [np.nan],
            'TAU_sd': [np.nan],
            'A_sd': [np.nan],
            'B_sd': [np.nan],
            'A': [np.nan],
            'B': [np.nan],
            'fit_err': [np.nan],
            'dof': [np.nan],
            'fit_type': ["dip2"],
        })
        
    curr_fit = pd.concat([curr_fit,curr_second_peak_to_end_fit],ignore_index= True)
    curr_fit['peak1_lat'] = d['peak1_lat']
    curr_fit['peak1_autocnsm'] = d['peak1_val']
    curr_fit['dip_lat'] = d['dip_lat']
    curr_fit['dip_autocnsm']= d['dip_countsnsm']
    curr_fit['peak2_lat'] =  d['peak2_lat']
    curr_fit['peak2_autocnsm'] = d['peak2_countsnsm']
    return curr_fit
    

## -- Subfunction 7 
## function to choose the best fit - fit_data takes the input as all the fits returned from subfunction 6 
def exp_fit_autoc_post_hoc(fit_data,unit): 

    #choosing the fits such that TAU<1000, A, B >0 and the covariance of parameters were estimated before.  so TAU_sd should not be infinity
    autoc_fit = fit_data[(fit_data['TAU'] < 1000) & (fit_data['TAU'] > 0) & (fit_data['A'] > 0) & (fit_data['B'] > 0)& ~np.isinf(fit_data['TAU_sd'])]

    if not autoc_fit.empty:
        # Keep the fit with the minimal error for each unit and fit type
        autoc_fit = autoc_fit.groupby(['unit', 'fit_type']).apply(lambda group: group.loc[group['fit_err'].idxmin()]).reset_index(drop=True)
    
        # Check if each fit type exists for each unit
        has_mono = (autoc_fit['fit_type'] == "mono").sum() == 1
        has_dip1 = (autoc_fit['fit_type'] == "dip1").sum() == 1
        has_dip2 = (autoc_fit['fit_type'] == "dip2").sum() == 1
    
        iunit = autoc_fit['unit'].unique()
    
        # Append missing fits
        if not has_mono:
            missing_mono = pd.DataFrame({
                'unit': iunit,
                'TAU': [np.nan],
                'TAU_sd': [np.nan],
                'A_sd': [np.nan],
                'B_sd': [np.nan],
                'A': [np.nan],
                'B': [np.nan],
                'fit_err': [np.nan],     
                'dof': [np.nan],
                'fit_type': ["mono"],
                'peak1_lat': [np.nan],
                'peak1_autocnsm': [np.nan],
                'dip_lat': [np.nan],
                'dip_autocnsm': [np.nan],
                'peak2_lat': [np.nan],
                'peak2_autocnsm': [np.nan]
            })
            autoc_fit = pd.concat([autoc_fit, missing_mono], ignore_index=True)
    
        if not has_dip1:
            missing_dip1 = pd.DataFrame({
                'unit': iunit,
                'TAU': [np.nan],
                'TAU_sd': [np.nan],
                'A_sd': [np.nan],
                'B_sd': [np.nan],
                'A': [np.nan],
                'B': [np.nan],
                'fit_err': [np.nan],     
                'dof': [np.nan],
                'fit_type': ["dip1"],
                'peak1_lat': [np.nan],
                'peak1_autocnsm': [np.nan],
                'dip_lat': [np.nan],
                'dip_autocnsm': [np.nan],
                'peak2_lat': [np.nan],
                'peak2_autocnsm': [np.nan]
            })
            autoc_fit = pd.concat([autoc_fit, missing_dip1], ignore_index=True)
    
        if not has_dip2:
            missing_dip2 = pd.DataFrame({
                'unit': iunit,
                'TAU': [np.nan],
                'TAU_sd': [np.nan],
                'A_sd': [np.nan],
                'B_sd': [np.nan],
                'A': [np.nan],
                'B': [np.nan],
                'fit_err': [np.nan],     
                'dof': [np.nan],
                'fit_type': ["dip2"],
                'peak1_lat': [np.nan],
                'peak1_autocnsm': [np.nan],
                'dip_lat': [np.nan],
                'dip_autocnsm': [np.nan],
                'peak2_lat': [np.nan],
                'peak2_autocnsm': [np.nan]
            })
            
            autoc_fit = pd.concat([autoc_fit, missing_dip2], ignore_index=True)
            
          
        # Pivot wider to get fit_err and dof for each fit type
        tbl=autoc_fit.pivot_table(index='fit_type', values=['fit_err','dof','unit'])
    
        dum = pd.DataFrame()
        # Calculate RMSE and best fit type
        dum['rmse_mono'] = [tbl.loc[('mono', 'fit_err')] / np.sqrt(tbl.loc[('mono', 'dof')])]
        dum['rmse_dip']= [np.sqrt(np.sum(np.array([tbl.loc[('dip1', 'fit_err')],tbl.loc[('dip2', 'fit_err')]])**2)/(tbl.loc[('dip1', 'dof')]+tbl.loc[('dip2', 'dof')]))]
    
       
        best_fit = np.select(
            [
                dum['rmse_dip'] < dum['rmse_mono'],
                ~dum['rmse_dip'].isna() & dum['rmse_mono'].isna(),
                ~dum['rmse_mono'].isna()
            ],
            ["dip", "dip", "mono"],
            "none"
        )
    
        # Select necessary columns
        dum['bestfit'] = [best_fit[0]]
        autoc_fit['bestfit'] = pd.Series([dum['bestfit'][0]] * len(autoc_fit))
        
    else:
        
        autoc_fit = pd.DataFrame({
            'unit': [unit],
            'TAU': [np.nan],
            'TAU_sd': [np.nan],
            'A_sd': [np.nan],
            'B_sd': [np.nan],
            'A': [np.nan],
            'B': [np.nan],
            'fit_err': [np.nan],     
            'dof': [np.nan],
            'fit_type': ["No valid fit"],
            'peak1_lat': [np.nan],
            'peak1_autocnsm': [np.nan],
            'dip_lat': [np.nan],
            'dip_autocnsm': [np.nan],
            'peak2_lat': [np.nan],
            'peak2_autocnsm': [np.nan],
            'bestfit' : ["No valid fit"] 
            })
            
   
    return autoc_fit
   

def mainFun_computeTAU(data,unit,lag_limit,excludeFirstBin):
    print('Inside main function to compute TAU')
    curr_autoc = compute_autoc(data,lag_limit)
    
    
    if curr_autoc is not None:
        
        print('Success compute_autoc')
    else:
        print('compute_autoc did not converge, moving on to next neuron')
        return None, None
    
    curr_fit = exp_fit_autoc(curr_autoc,unit,excludeFirstBin)
    print('Success estimating all fits')
    best_fit = exp_fit_autoc_post_hoc(curr_fit,unit)
    print('Success choosing best fit')
    
    return best_fit, curr_autoc


        
        
        
        