#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 22:34:46 2023

@author: nilanjana
"""

"""
Inputs

arg1 - takes input as either 1 or 0. 1 means the data over the session will be aligned trialwise
arg2 - entire spike time series data as a vector
arg3 - if arg1 is 1 then the start time for the alignment for each trial
arg4 - if arg1 is 1 then the end time for the alignment for each trial


arg1 or alignData is 1 when i want to save the firing rates between Sel and Go onset in each trial. 
Use this when i have extracted the TAU trialwise, that is, comparison between preSC1 and preGo 

arg1 is 0 when i am saving the FR for the entire session. This is saved like this when i want to this FR and 
the TAU estimated from the entire session together.


"""



#function to compute firing rate in specific time periods and save them. clean this code later

import numpy as np



def compute_firingRate(alignData,spike_data,t_st,t_end): #give input of onset timings of the two periods in which firing rate is to bw calculated 
    
    
    # Initialize an empty list to store firing rate of neuron for each trial
    firing_rate = []
    
    if alignData ==1:
        
        # Iterate through trials and align spike data
        for trial_start, trial_end in zip(t_st, t_end):
            delta_t = (trial_end - trial_start)#converting it to seconds
           
            # Extract spike data within the trial's start and end times
            aligned_spikes = spike_data[(spike_data >= trial_start) & (spike_data < trial_end)]
            
            # Align the extracted spike data by subtracting the trial start time
            aligned_spikes = aligned_spikes - trial_start
            n =len(aligned_spikes) 
            #print(n)
            #print('One trial done')
            fr_perTrial =( n/delta_t) *1000   #firing rate per Trial spikes/sec
            
            # Append the aligned spike data to the list
            firing_rate.append(fr_perTrial)
         
        
        mu_FR = np.mean(firing_rate)
        
    else:
        n = len(spike_data) #total number of spikes in the recording duration
        deltaT = (spike_data[len(spike_data)-1]- spike_data[0])#converting it back to spikes/sec
        mu_FR = (n/deltaT) *1000 #converting it back to spikes/sec
        
    return (mu_FR)




