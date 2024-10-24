#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:03:36 2023

@author: nilanjana
"""
from scipy.io import loadmat


def open_matlab_behaviour(filename):
    """
    Open a matlab behaviour structure from OrganizeBehaviour and save it in a python dictionary. One
    SESSION corresponds to one behavioural file - independently of the number of probes.
    :param filename: name of the matlab structure containing the behavioural information
    ----
    :return: behaviour: dictionary containing the same fields as the behaviour structure from matlab
    """

    # import the matlab file
    mat_behaviour = loadmat(filename)
    
    # get the sub-structure name
    data_str_name = list(mat_behaviour.keys())[-1]
    
    # get the fields of the structure and use them as keys for the new dictionary
    fields = []
    for key in mat_behaviour[data_str_name].dtype.fields.keys():
        fields.append(key)
    
    data = mat_behaviour[data_str_name][0][0]  # get the data inside the structure
    behaviour = {field: data[i_field][:, 0] for i_field, field in enumerate(fields)}
    
    return behaviour