#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 11:09:50 2024

@author: nilanjana

skipped Mo180426004 due to extraction failure in MATLAB

This wont work !!! Change it  !!

"""


#import modules
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
import os
from scipy.io import loadmat 
from sklearn.model_selection import KFold
import pickle
import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt



# def decode_directions_with_svm(data, labels, n_classes, n_bins, n_splits, bin_size):
#     """
#     Perform decoding of movement directions using an SVM with k-fold cross-validation.
    
#     Args:
#         data (np.array): Pseudopopulation data of shape (n_classes * min_trials, num_neurons, timepoints).
#         labels (np.array): Class labels for each trial.
#         n_classes (int): Number of classes.
#         n_bins (int): Number of time bins.
#         n_splits (int): Number of splits for cross-validation.
#         bin_size (int): Size of each time bin.
        
#     Returns:
#         list: List of accuracies for each fold.
#     """
#     # Compute firing rates per trial in each time bin
#     fr_data = compute_firing_rate(data, bin_size=bin_size)  # Shape: (total_trials, num_neurons, n_bins)
    
#     # Reshape for SVM: (total_trials, num_neurons * n_bins)
#     fr_data = fr_data.reshape(fr_data.shape[0], -1)
    
#     # Initialize k-fold cross-validation
#     kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
#     accuracies = []
    
#     # Perform k-fold cross-validation
#     for train_index, test_index in kf.split(fr_data):
#         train_data, test_data = fr_data[train_index], fr_data[test_index]
#         train_labels, test_labels = labels[train_index], labels[test_index]
        
#         # Initialize the SVM classifier
#         svm_classifier = SVC(kernel='linear')
        
#         # Train the SVM on the training data
#         svm_classifier.fit(train_data, train_labels)
        
#         # Predict on the test set
#         test_predictions = svm_classifier.predict(test_data)
        
#         # Compute accuracy
#         accuracy = accuracy_score(test_labels, test_predictions)
#         accuracies.append(accuracy)
    
#     return accuracies

def neurons_list_minC(data,ttype,minC):
    '''
    arg1 - dict -- full dictionary data
    arg2: int -- trial type - 1/2/3
    minC: int -- minimum number of trials in each direction for each  neuron required

    Returns
    -------
    indices of neurons which has minimum number of trials in all unique directions
    

    '''
    
    # Create a dictionary to store the results
   
    dir_key = list(data[1].keys())
    neuron_results = {}
      
    # Iterate over each neuron
    for n in range(len(data[ttype][dir_key[0]])):  # Assuming all directions have the same number of neurons
       
        neuron_results[n] = {}
        
        # Check each direction for the number of trials
        for direction in dir_key:
            num_trials = len(data[ttype][direction][n])  # Get the number of trials for neuron `n` in `direction`
            
            # Check if the number of trials is greater than or equal to `min_C`
            if num_trials >= minC:
                neuron_results[n][direction] = 1  # Passes the condition
            else:
                neuron_results[n][direction] = 0  # Fails the condition
    
    # Convert the results dictionary to a DataFrame
    df_results = pd.DataFrame.from_dict(neuron_results, orient='index')
   
    minC_allDir = df_results[(df_results == 1).all(axis=1)]
    # Get the indices of these neurons
    idx_neurons = minC_allDir.index
    
    return idx_neurons


if __name__ == "__main__":
    
    monkey = 'Mourad'
    probe = 'PMd'
    
    bin_sz= 100
    sl_wd = 0
    
    #need to give the path to find each neuron information
    current_path = os.getcwd()
    if current_path.startswith('/Users'):
        server = '/Volumes'  # local w VPN
    elif current_path.startswith('/home/'):
        server = '/envau'  # niolon
    elif current_path.startswith('/envau'):
        server = '/envau'  # niolon
        
    
    path  = f'{server}/work/comco/nandi.n/IntrinsicTimescales/data/FR_data/SEL/{monkey}/Quality1/{probe}/AllLayers/FR_ST_bin_{bin_sz}_sl_{sl_wd}.pkl'   
    with open (path, 'rb') as f:
        data = pickle.load(f)
    
    ttype_key = list(data.keys())
    dir_key = list(data[1].keys())
    
    timebins = len(data[1][dir_key[0]][0][0]) #[trialtype][direction][neuronNum][trial1] 
    
    min_C = 12

    outpath = f'{server}/work/comco/nandi.n/IntrinsicTimescales/results/bothMonkeys/DecodingAnalysis'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        
    for ttype in ttype_key:
        data_ttype = data[ttype]
        
   #first choose the neurons which have min trials for all 4 directions - let's say 8 per direction  and save the list of those neurons
        idx_nc = neurons_list_minC(data,ttype,min_C)
        print (f'Number of neurons with {min_C} trials for all direction for ttype {ttype} is {len(idx_nc)}')
       
        # Initialize a 3D array to store cross-temporal accuracies for all iterations
        all_iterations_acc = np.zeros((100, timebins, timebins))  
    
    # Run cross-validation 100 times
        for iteration in range(10):
            print(f'Iteration {iteration + 1}/100')
    
            # Initialize a matrix to store cross-temporal accuracies
            cross_temporal_acc = np.zeros((timebins, timebins))
    
            # Create a single train-test split for this iteration
            train_data_dict = {direction: [] for direction in data_ttype.keys()}
            test_data_dict = {direction: [] for direction in data_ttype.keys()}
            for direction in data_ttype.keys():
                for neuron in idx_nc:
                    trials = data_ttype[direction][neuron]
                    
                    # Split trials into training and testing sets without overlap
                    train_trials, test_trials = train_test_split(trials, train_size=10, random_state=iteration)
                    train_data_dict[direction].extend(train_trials)
                    test_data_dict[direction].extend(test_trials)
    
            # Iterate over time bins for training
            for train_tm in range(timebins):
                # Prepare training data and labels for the current training time bin
                train_data = []
                train_labels = []
                for direction in train_data_dict.keys():
                    train_data.extend([trial[train_tm] for trial in train_data_dict[direction]])
                    train_labels.extend([direction] * len(train_data_dict[direction]))
    
                # Ensure data is a 2D array (samples, features)
                train_data = np.array(train_data).reshape(-1, 1)
    
                # Train an SVM classifier
                svm = SVC(kernel='linear')
                svm.fit(train_data, train_labels)
    
                # Iterate over time bins for testing
                for test_tm in range(timebins):
                    # Prepare testing data and labels for the current testing time bin
                    test_data = []
                    test_labels = []
                    for direction in test_data_dict.keys():
                        test_data.extend([trial[test_tm] for trial in test_data_dict[direction]])
                        test_labels.extend([direction] * len(test_data_dict[direction]))
    
                    # Ensure data is a 2D array (samples, features)
                    test_data = np.array(test_data).reshape(-1, 1)
    
                    # Predict the direction of test trials
                    predictions = svm.predict(test_data)
    
                # Calculate and store accuracy for this train-test pair
                accuracy = accuracy_score(test_labels, predictions)
                all_iterations_acc[iteration, train_tm, test_tm] = accuracy
    
            # Print or save the cross-temporal accuracy matrix
            print(f'Cross-temporal accuracy matrix for iteration {iteration + 1}:')
            print(cross_temporal_acc)
        
                # Calculate the average cross-temporal accuracy across all iterations
        mean_cross_temporal_acc = np.mean(all_iterations_acc, axis=0)
        
        # Plot the cross-temporal prediction map
        plt.figure(figsize=(10, 8))
        plt.imshow(mean_cross_temporal_acc, cmap='viridis', vmin=0.2, vmax=0.8, aspect='auto')
        
        
        t_markers = list(np.array([150, 1450, 2750, 4050]) + 1200)

        t_markers_labels = ['SEL', 'SC1', 'SC2', 'SC3']
        # vertical lines
        t_lines = list(
            np.array([0, 300, 1300, 1600, 2600, 2900, 3900, 4200]) + 1200)
        
                # Add vertical and horizontal lines
        for t in t_lines:
            plt.axvline(x=t, color='white', linestyle='--', linewidth=1)
            plt.axhline(y=t, color='white', linestyle='--', linewidth=1)
        
        # Add time marker labels at the appropriate points
        for i, t in enumerate(t_markers):
            plt.text(t, -20, t_markers_labels[i], color='white', ha='center', va='bottom', fontsize=10)  # X-axis labels
            plt.text(-20, t, t_markers_labels[i], color='white', ha='right', va='center', fontsize=10, rotation=0)  # Y-axis labels

      
        
        plt.tick_params(axis = 'both',labelsize = 18)
        
        plt.colorbar(label='Accuracy')
        plt.title(f'Cross-Temporal Decoding Accuracy Map for {ttype} and {monkey} {probe}')
        plt.xlabel('Testing Time Bin')
        plt.ylabel('Training Time Bin')
        plt.tight_layout()
        
        plt.show()
        
        
        plt.savefig(os.path.join(outpath,f'{ttype}.png'))
        
   

