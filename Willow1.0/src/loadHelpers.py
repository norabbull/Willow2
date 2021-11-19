# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 14:58:51 2021

@author: norab
"""


import pandas as pd

def load_simData(file='C:/Users/norab/Master/data/simulation/simNull_pvals.csv'):
    return pd.read_csv(file)