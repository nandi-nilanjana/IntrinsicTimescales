#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 10:23:00 2024

@author: nilanjana
Created this file intially to combine the information from two separate excel files - 
combining AP/LAT file of each monkey with the PMd/M1 information file  - earlier as Mourad/TomyLaminarInfo.xlsx and Mourad/Tomy_AreaInfo.csv


"""

#modules import

import pandas as pd
import os 
import re
import numpy as np

monkey = ['Mourad','Tomy']

#current server through which the code runs
current_path = os.getcwd()
if current_path.startswith('/Users'):
    server = '/Volumes'  # local w VPN
elif current_path.startswith('/home/'):
    server = '/envau'  # niolon
elif current_path.startswith('/envau'):
    server = '/envau'  # niolon

for i in len(monkey):
    
    #import file 1 
    df_1 = pd.read_excel(f'{server}/work/comco/nandi.n/IntrinsicTimescales/docs/{monkey[i]}LaminarInfo.xlsx')
    df_2 = pd.read_csv(f'{server}/work/comco/nandi.n/IntrinsicTimescales/docs/{monkey[i]}_AreaInfo.csv')
    
    merged_df = pd.merge(df_1, df_2[['Session', 'Probe_Num','AP', 'LAT']], on=['Session', 'Probe_Num'], how='left')  
    
    # Filter rows in df_2 that are not in df_1 based on 'Session' and 'probe_num'
    extra_in_df2 = df_2[~df_2.set_index(['Session', 'Probe_Num']).index.isin(df_1.set_index(['Session', 'Probe_Num']).index)]
    
    outFile = f'{server}/work/comco/nandi.n/IntrinsicTimescales/docs/{monkey[i]}updated.xlsx'
    merged_df.to_excel(outFile,index=False)

