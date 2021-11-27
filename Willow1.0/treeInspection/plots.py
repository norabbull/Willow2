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
folder = 'C:/Users/norab/Master/thesis_data/'

nz_phydists = load_totdist(file=folder + 'meta_data/test_nonzero_phydists.csv')  # works
uniqseqs = load_uniqseqs(file=folder + 'meta_data/9381_uniqseqs.txt') #works
GDRs = load_GDRs(file = folder + 'test_result_data/GDR_10genes_all.csv') # works

# Maybe leave out
GDRsuper = load_GDRs(group_category = 'super')
GDRsub = load_GDRs(group_category = 'sub')

data_all = nz_phydists.merge(uniqseqs, on = 'gene')
data_all = data_all.merge(GDRs, on = 'gene')
data_all.drop_duplicates(inplace=True)

#%%

