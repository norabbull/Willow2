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

#%% Signifiacnt genes

data_all_sig_SUB = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRsignificant/SDR1005_significantGenesSUB.csv')
data_all_sig_SUPER = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRsignificant/SDR1005_significantGenesSUPER.csv')

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

below_09 = GDRsuper[GDRsuper['GDR'] < 0.9]
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
 + geom_line(m2, color='saddlebrown', size=0.4, alpha=0.3)
 + scale_fill_manual(values=['indigo', 'darkorange'])
 + theme_classic()  # change color
 + theme(figure_size=(8, 6))
 + labs(title='GDR distributions')
 + xlab("Group level")
 + ylab("GDR")
 + theme_classic()
 + theme(
     axis_text_x=element_text(rotation=30, hjust=1),
         plot_title=element_text(color = "black",
                             size = 35,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
          axis_title=element_text(color='black',
                             size=25,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
          axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=25,
     legend_key_height=15,
     legend_key_size=25,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     )
)

#%% MAL

(ggplot(all_simData, aes('permutations', 'group_size', fill = 'adj p-val'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
 + theme(figure_size=(13, 13))
 + labs(title='GDR simulation: adjusted p - values')
 + scale_x_continuous(name="Number of permuted G2-samples from C2 to C1")
 + scale_y_continuous(name="Group size  (Number of samples in G)")
 + scale_color_cmap(cmap_name="inferno")
 + theme_classic()
  + theme(
     plot_title=element_text(color = "black",
                             size = 25,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
     axis_title=element_text(color='indigo',
                             size=15,
                             family='serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
     
     axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=30,
     legend_key_height=15,
     legend_key_size=30,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     #legend_title_align='center',
     
     )
)




# =============================================================================
# Plot density + boxplot
# =============================================================================


(ggplot(GDRs, aes('level', 'GDR', fill = 'level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.9, show_legend = False)
 + geom_boxplot(width = shift, alpha=0.4, size = 0.9, show_legend = False)
 + scale_fill_manual(values=['indigo', 'darkorange'])
 + theme_classic()  # change color
 + theme(figure_size=(8, 6))
 + labs(title='GDR distributions')
 + xlab("Group level")
 + ylab("GDR")
 + theme(
     axis_text_x=element_text(rotation=30, hjust=1),
         plot_title=element_text(color = "black",
                             size = 35,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
          axis_title=element_text(color='black',
                             size=25,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
          axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=25,
     legend_key_height=15,
     legend_key_size=25,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     )
)



# =============================================================================
# GDR desnity plot, super + sub 
# =============================================================================


(ggplot(GDRs, aes(x='GDR', fill='level'))
 + geom_density(adjust = 1/2, alpha=0.5) 
 + scale_fill_manual(values=['indigo', 'darkorange'])
 + theme_classic()  # change color
 + theme(figure_size=(8, 6))
 + labs(title='GDR distributions')
 + xlab("GDR")
 + ylab("Density")
 + theme(
     axis_text_x=element_text(rotation=30, hjust=1),
         plot_title=element_text(color = "black",
                             size = 35,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
          axis_title=element_text(color='black',
                             size=25,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
          axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=30,
     legend_key_height=20,
     legend_key_size=40,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold', size=12),
     legend_position='right'

     )
)



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
+ scale_x_continuous(name="GDR")
+ scale_y_continuous(name="density")
  + theme(
     plot_title=element_text(color = "black",
                             size = 32,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
     axis_title=element_text(color='black',
                             size=22,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
     
     axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=30,
     legend_key_height=15,
     legend_key_size=30,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     #legend_title_align='center',
)
  )


# MAL
(ggplot(all_simData, aes('permutations', 'group_size', fill = 'p-val'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
 + theme(figure_size=(7, 7))
 + labs(title='GDR simulation: p - values')
 + scale_x_continuous(name="Number of G2-samples moved from C2 to C1")
 + scale_y_continuous(name="Group size  (number of samples in G)")
  + theme(
     plot_title=element_text(color = "black",
                             size = 32,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
     axis_title=element_text(color='black',
                             size=22,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
     
     axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=30,
     legend_key_height=15,
     legend_key_size=30,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     #legend_title_align='center',
     
     )
)

# SUB

# l = [pd.read_csv(filename) for filename in glob.glob("C:/Users/norab/Master/thesis_data/result_data/GDRnull/sub/*.csv")]
# all_GDRnull = pd.concat(l, axis=0)
# all_GDRnull.describe()

# half_GDRnull = all_GDRnull.iloc[0:30000]
# half_GDRnull.columns = ['GDRnull']


# (
# ggplot(half_GDRnull, aes(x='GDRnull', y=after_stat('density')))
# + geom_histogram(
#     colour='orange', # change the outline
#     size=2,        # change the thickness of the outline
#     alpha=0.7      # change the transparency
#     )
# + labs(title='30.000 random GDR values (sub)')
# )

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

#%% MAL


(ggplot(all_simData, aes('permutations', 'group_size', fill = 'adj p-val'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
 + theme(figure_size=(13, 13))
 + labs(title='GDR simulation: adjusted p - values')
 + scale_x_continuous(name="Number of permuted G2-samples from C2 to C1")
 + scale_y_continuous(name="Group size  (Number of samples in G)")
 + scale_color_cmap(cmap_name="inferno")
  + theme(
     plot_title=element_text(color = "black",
                             size = 25,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
     axis_title=element_text(color='indigo',
                             size=15,
                             family='serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
     
     axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=30,
     legend_key_height=15,
     legend_key_size=30,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     #legend_title_align='center',
     
     )
)

#%%
nz_phydists = nz_phydists.dropna()

(
ggplot(nz_phydists, aes(x='totdist', y=after_stat('count/sum(count)*100')))
+ geom_histogram(
    colour='darkgreen', # change the outline
    size=2,        # change the thickness of the outline
    alpha=0.7      # change the transparency
    )
+ labs(title='Non - zero values in distance matrices')
 + scale_x_continuous(name="% non-zero distance-values in distance matrices")
 + scale_y_continuous(name="% of all distance matrices")
  + theme(
     plot_title=element_text(color = "black",
                             size = 32,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
     axis_title=element_text(color='black',
                             size=22,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
     
     axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=30,
     legend_key_height=15,
     legend_key_size=30,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     #legend_title_align='center',
)
  )
