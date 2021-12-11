# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 06:45:58 2021

@author: norab
"""

from treeHelpers import *
import pandas as pd
import os
from os.path import isfile, join

import numpy as np

from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine import *  # TO DO: change
import numpy as np
import matplotlib.pyplot as plt
#from functools import reduce 
import math

#%%
# Load data 

# Folder all files are stored in : 
folder = 'C:/Users/norab/Master/thesis_data/'

nz_phydists = load_nz_phydists(file=folder + 'meta_data/nonzero_phydists.csv')  # works
uniqseqs = load_uniqseqs(folder + 'meta_data/9381_uniqseqs.txt') #works
GDRs = load_GDRs(file = folder + 'result_data/GDR/GDR_all.csv') # works
GDRnull = load_GDRnull()   # works
GDRsuper = load_GDRs(file = folder + 'result_data/GDR/GDR_all.csv', group_category = 'super')
GDRsub = load_GDRs(file = folder + 'result_data/GDR/GDR_all.csv', group_category = 'sub')

data_all = nz_phydists.merge(uniqseqs, on = 'gene')
data_all = data_all.merge(GDRs, on = 'gene')
data_all.drop_duplicates(inplace=True)

#%%
# all_data, super

data_all_super = nz_phydists.merge(uniqseqs, on = 'gene')
data_all_super = data_all_super.merge(GDRsuper, on = 'gene')
data_all_super.drop_duplicates(inplace=True)

# all_data, sub

data_all_sub = nz_phydists.merge(uniqseqs, on = 'gene')
data_all_sub = data_all.merge(GDRsub, on = 'gene')
data_all.drop_duplicates(inplace=True)


#%% Stats

# Get basic statistics of GDR values
GDRs.describe()
GDRsuper.describe()
GDRsub.describe()

below_095 = GDRsuper[GDRsuper['GDR'] < 0.]
# Inspect which genes that did not obtain GDR value
# 1. Get list of all GDR gene labels
# 2. Get list of all input tree labels
# 3. Check non-zero dist 

# 
zero_phydists = nz_phydists[nz_phydists['totdist'] == 0]
nan_phydists = nz_phydists[nz_phydists['totdist'].isnull()]

one_phydists = GDRs[GDRs['GDR'] == 1]
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
# Plot desnity + lines between sub and super
# =============================================================================

(ggplot(GDRs, aes('level', 'GDR', fill = 'level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.8, show_legend = False)
 + geom_point(m2, color='none', alpha=0.6, size=1.5, show_legend=False)
 + geom_line(m2, color='gray', size=0.65, alpha=0.6)
 + scale_fill_manual(values=['indigo', 'darkorange'])
 + theme_classic()  # change color
 + theme(figure_size=(8, 6))
 + labs(title='GDR')
)


# =============================================================================
# Plot density + boxplot
# =============================================================================

(ggplot(GDRs, aes('level', 'GDR', fill = 'level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
 + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
 + scale_fill_manual(values=['indigo', 'darkorange'])
 + theme_classic()  
 + theme(figure_size=(8, 6))
 + labs(title='GDR for all genes')
)


# =============================================================================
# GDR desnity plot, super + sub 
# =============================================================================


(ggplot(GDRs, aes(x='GDR', fill='level')) 
     + geom_density(adjust = 1/2, alpha=0.5) 
     + scale_fill_manual(values=['indigo', 'darkorange'])
     + labs(title='GDR for all genes'))


#%% Combined plot of all null distributions


# =============================================================================
# GDR for significant genes
# =============================================================================

# Inspect GDRnull distribution values 

import glob

l = [pd.read_csv(filename) for filename in glob.glob("C:/Users/norab/Master/thesis_data/result_data/GDRnull/super/*.csv")]
all_GDRnull = pd.concat(l, axis=0)
all_GDRnull.describe()

half_GDRnull = all_GDRnull.iloc[0:30000]
half_GDRnull.columns = ['GDRnull']


(
ggplot(half_GDRnull, aes(x='GDRnull', y=after_stat('density')))
+ geom_histogram(
    colour='darkgreen', # change the outline
    size=2,        # change the thickness of the outline
    alpha=0.7      # change the transparency
    )
+ labs(title='30.000 random GDR values (super)')
)


# SUB

l = [pd.read_csv(filename) for filename in glob.glob("C:/Users/norab/Master/thesis_data/result_data/GDRnull/sub/*.csv")]
all_GDRnull = pd.concat(l, axis=0)
all_GDRnull.describe()

half_GDRnull = all_GDRnull.iloc[0:30000]
half_GDRnull.columns = ['GDRnull']


(
ggplot(half_GDRnull, aes(x='GDRnull', y=after_stat('density')))
+ geom_histogram(
    colour='orange', # change the outline
    size=2,        # change the thickness of the outline
    alpha=0.7      # change the transparency
    )
+ labs(title='30.000 random GDR values (sub)')
)

#%% Other plots


# Unique sequences vs percentage of non-zero phylogenetic distances in tree

(ggplot(data_all, aes('uniqseq', 'totdist', fill = 'GDR'))
 + geom_point()
 + theme_classic()
 + labs(title='non-zero phydist percent vs uniqseq')
)

# =============================================================================
# La stå
# =============================================================================
# Full
(ggplot(data_all, aes('uniqseq', 'totdist', fill = 'SDRsub'))
 + geom_point(alpha = 0.7, stroke= .05, size = 3)
 + theme_classic()
 + labs(title='SUB: non-zero sample distances VS nr. unique sequences')
)

(ggplot(result, aes('uniqseq', 'totdist', fill = 'SDR'))
 + geom_point(alpha = 0.7, stroke= .05, size = 3)
 + theme_classic()
 + labs(title='SUPER: non-zero sample distances VS nr. unique sequences')
)

# Xlimited

(ggplot(result, aes('uniqseq', 'totdist', fill = 'SDRsub'))
 + geom_point(alpha = 0.7, stroke= .05, size = 4)
 + theme_classic()
 + labs(title='SUB: non-zero sample distances VS nr. unique sequences')
 + xlim(0,30)
)

(ggplot(result, aes('uniqseq', 'totdist', fill = 'SDR'))
 + geom_point(alpha = 0.7, stroke= .05, size = 4)
 + theme_classic()
 + labs(title='SUPER: non-zero sample distances VS nr. unique sequences')
 + xlim(0,30)
)


(ggplot(result, aes('uniqseq', 'totdist', fill = 'SDR'))
 + geom_point(alpha = 0.7, stroke= .01, size = 5 )
 + theme_classic()
 + labs(title='SUPER: totdist vs uniqseq')
 + xlim(0,30)
)



