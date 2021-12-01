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
from inspection.inspectionHelpers import *   # TO DO: change
from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine import *  # TO DO: change
import numpy as np
import matplotlib.pyplot as plt
#from functools import reduce 

#%%
# Load data 

# Folder all files are stored in : 
folder = 'C:/Users/norab/Master/thesis_data/'

nz_phydists = load_nz_phydists(file=folder + 'meta_data/test_nonzero_phydists.csv')  # works
uniqseqs = load_uniqseqs(folder + 'meta_data/9381_uniqseqs.txt') #works
GDRs = load_GDRs(file = folder + 'result_data/GDR/GDR_all.csv') # works
GDRnull = load_GDRnull()   # works

# Maybe leave out
GDRsuper = load_GDRs(group_category = 'super')
GDRsub = load_GDRs(group_category = 'sub')

data_all = nz_phydists.merge(uniqseqs, on = 'gene')
data_all = data_all.merge(GDRs, on = 'gene')
data_all.drop_duplicates(inplace=True)

#%%

# =============================================================================
# Plot interconnecting GDR for super and sub
# =============================================================================


shift = 0.1

def alt_sign(x):
    "Alternate +1/-1 if x is even/odd"
    return (-1) ** x

m1 = aes(x=stage('level', after_scale='x+shift*alt_sign(x)'))              # shift outward
m2 = aes(x=stage('level', after_scale='x-shift*alt_sign(x)'), group='gene')  # shift inward


#
# =============================================================================
# Violine + points + lines between sub and super
# =============================================================================

(ggplot(data_all, aes('level', 'GDR', fill = 'level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.8, show_legend = False)
 + geom_point(m2, color='none', alpha=0.6, size=1.5, show_legend=False)
 + geom_line(m2, color='gray', size=0.65, alpha=0.6)
 + theme_classic()  # change color
 + theme(figure_size=(8, 6))
 + labs(title='SDR values for 8782 genes (all)')
)


# =============================================================================
# Plot interconnecting GDR for super and sub
# =============================================================================



# Violin + points + boxplot
(ggplot(data_all, aes('level', 'GDR', fill = 'level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
 + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
 + scale_fill_manual(values=['dodgerblue', 'darkorange'])
 + theme_classic()   # Change color
 + theme(figure_size=(8, 6))
 + labs(title='SDR values for 8782 genes (all)')
)


# boxplot

(ggplot(data_all, aes('level', 'GDR', fill = 'level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
 + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
 + scale_fill_manual(values=['indigo', 'darkorange'])
 + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='GDR values for all genes')
)
# =============================================================================
# GDR desnity plot, super + sub 
# =============================================================================


(ggplot(data_all, aes(x='GDR', fill='level')) 
     + geom_density(adjust = 1/2, alpha=0.5) 
     + scale_fill_manual(values=['indigo', 'darkorange'])
     + labs(title='GDRs'))