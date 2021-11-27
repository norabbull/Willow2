# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 06:45:58 2021

@author: norab
"""

from treeHelpers import *
import pandas as pd
import os
from os.path import isfile, join
import re
import shutil
import numpy as np
import copy


#%%
# Load data 

# Folder all files are stored in : 
folder = 'C:/Users/norab/Master/thesis_data/meta_data/'

nz_phydists = load_totdist(file=folder + 'test_nonzero_phydists.csv')  # works
uniqseqs = load_uniqseqs(file=folder + '9381_uniqseqs.txt') #works
GDRs = load_GDRs()

data_all = totdist.merge(SDRs, on = 'gene')
data_all = data_all.merge(SDVs, on = ['gene', 'level'])
data_all = data_all.merge(uniqseqs, on = ['gene'])
data_all.drop_duplicates(inplace=True)

#%%

