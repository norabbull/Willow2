# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 08:28:05 2021

@author: norab

Inspect SDR, SDV, uniqseq and total non-dist

"""

from inspections.inspectionHelpers import *
from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce 

# =============================================================================
# Load stuff 
# =============================================================================
uniqseqs = load_uniq_seqs()
totdist = load_tot_dist()
SDRs = load_SDRs()
SDVs = load_SDVs()
SDRnull = load_SDRnull()


# Single SDRs	gene
SSDR_NFYA = load_singleSDR("ENSG00000110700___RS13")
SSDR_FGR =load_singleSDR("ENSG00000000938___FGR")
random_trees_50 = load_allSingleSDRs()

# Null dist data
nullDistSuper = load_nullDist(level = "super")
nullDistSub = load_nullDist(level = "sub")

# =============================================================================
# Merge into common df
# =============================================================================

# Merge
data1 = totdist.merge(uniqseqs, on = 'gene')
data2 = data1.merge(SDRs, on  = 'gene', how = 'outer')
data_all = data2.merge(SDVs, on = ['gene', 'level'], how = 'outer')

# Wrangle
data_all.drop_duplicates(inplace=True)
data_all.dropna(inplace=True)

# Sort 
data_all = data_all.sort_values(by=['SDR'])


# =============================================================================
# Plotting
# =============================================================================

# =============================================================================
# X vs Y, dotplot, lineplot
# =============================================================================


(ggplot(data_all, aes('SDR', 'uniqseq', fill = 'level'))
 + geom_point()
 + theme_classic()
 + labs(title='SDR vs uniqseq')
)

(ggplot(data_all, aes('SDR', 'totdist', fill = 'level'))
 + geom_point()
 + theme_classic()
 + labs(title='SDR vs uniqseq')
)

(ggplot(data_all, aes('uniqseq', 'totdist', fill = 'level'))
 + geom_point()
 + theme_classic()
 + labs(title='SDR vs uniqseq')
)

(ggplot(data_all, aes('SDR', 'SDV', fill = 'level'))
 + geom_point()
 + theme_classic()
 + labs(title='SDR vs uniqseq')
)

# =============================================================================
# Good inspiration: 
#    https://dputhier.github.io/jgb53d-bd-prog_github/practicals/intro_ggplot/intro_ggplot.html
# =============================================================================


# =============================================================================
# Violin + points + boxplot
# =============================================================================


shift = 0.1

def alt_sign(x):
    "Alternate +1/-1 if x is even/odd"
    return (-1) ** x

m1 = aes(x=stage('level', after_scale='x+shift*alt_sign(x)'))              # shift outward
m2 = aes(x=stage('level', after_scale='x-shift*alt_sign(x)'), group='gene')  # shift inward


(ggplot(data_all, aes('level', 'SDR', fill = 'level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
 + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
 + scale_fill_manual(values=['dodgerblue', 'darkorange'])
 + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR values for 8782 genes (all)')
)


# =============================================================================
# Line plot SDR vs uniqseqs
# =============================================================================

(ggplot(allData,aes('SDR','uniqseq')) 
  + geom_line()
  #+ scale_x_continuous(breaks=us_td_sdr[''])
)

(ggplot(data_all,aes('SDR','totdist')) 
  + geom_line()
  #+ scale_x_continuous(breaks=us_td_sdr[''])
)


(ggplot(allData,aes('SDR','uniqseq',color='level'))
  + geom_line(alpha = 0.8)
  + labs(title='SDRs vs uniqseq (all)')
  #+ theme_light()
)


# =============================================================================
# 
# =============================================================================



# =============================================================================
#  Overlapping distribution plot
# =============================================================================

density = ggplot(data=data_all, 
       mapping=aes(x='SDR', fill='level')) + geom_density(adjust = 1/4, alpha=0.5)

ggplot_build()

# =============================================================================
# Overlapping histogram plot
# =============================================================================
 
ggplot(allData) + 
aes(x="SDR", fill = "level") + stat_bin(bins=100) + geom_bar()


# =============================================================================
# Cumulative distribution function plot
# =============================================================================

# Extraxt levels 
SDRsuper = allData[allData['level'] == "super"]
SDRsub = allData[allData['level'] == "sub"]
SDRpsuedo = allData[allData['level'] == "psuedo"]

  
# No of data points used
N = 500
  
# normal distribution
data = np.random.randn(N)
  
# sort the data in ascending order
x = np.sort(data)
x_super = np.sort(SDRsuper['SDR'])
x_sub = np.sort(SDRsub['SDR'])
x_psuedo = np.sort(SDRpsuedo['SDR'])
  
# get the cdf values of y
y = np.arange(N) / float(N)
y_super = np.arange(len(x_super)) / float(len(x_super))
y_sub = np.arange(len(x_sub)) / float(len(x_sub))
y_psuedo = np.arange(len(x_psuedo)) / float(len(x_psuedo))


# plotting
plt.figure()
plt.xlabel('x-axis')
plt.ylabel('y-axis')
  
plt.title('Cumulative distribution plot')
  
plt.plot(x_super, y_super, marker='_', label = 'super')
plt.plot(x_sub, y_sub, marker='_', label = 'sub', color = 'red')
plt.plot(x_psuedo, y_psuedo, marker = '_', label = 'psuedo')
plt.legend()
plt.xlabel("SDR")
plt.show()

# =============================================================================
# Trashbin
# =============================================================================



# Sort
allData = us_td_sdr_null.sort_values(by=['SDR'])
genes = allData['gene']

(ggplot(allData,aes(x = 'SDR') 
        + geom_histogram()
        + xlim(0.5, 1.5)
))


(
    ggplot(allData, aes(x='SDR', color = 'level', fill = 'level'))
    + geom_histogram(alpha = 0.5, size = 0.7 )
    + xlim(0.2, 1.5)
)


(
    ggplot(SDRs_ps, aes(x='SDR', color = 'level', fill = 'level'))
    + geom_histogram(alpha = 0.7, size = 0.7)
    + xlim(0.2, 1.5)
)


(
    ggplot(SDRs_ps, aes(x='SDR', color = 'level')
    + geom_histogram()
))







