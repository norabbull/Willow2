# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 14:58:51 2021

@author: norab
"""


import pandas as pd
import os
from os.path import isfile, join

def make_filelist(self, input_files):
    """
    Input: 
        input_files: path to folder
        type: str
        
    Function: 
        creates list of full path to every file in folder
        
    """
    try: 
        if isinstance(input_files, str):
            if '.csv' in input_files: 
                files = pd.read_csv(input_files, header = None)
                files = files[0].values.tolist()
            else: 
                files = [join(input_files, f) for f in os.listdir(input_files) 
                             if isfile(join(input_files, f))]
                files = [f.strip() for f in files]
    except Exception as e:
        print("Error with make_filelist:")
        print(e)
            
    return files

def load_simData(file='C:/Users/norab/Master/data/simulation/simNull_pvals.csv'):
    return pd.read_csv(file)