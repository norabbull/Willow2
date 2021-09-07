# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 08:28:05 2021

@author: norab


    - Inspect SDR, SDRnull
    - Apply filters, inspect again
    - Find SDR treshold

"""

from inspection.inspectionHelpers import *   # TO DO: change
from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine import *  # TO DO: change
import numpy as np
import matplotlib.pyplot as plt
#from functools import reduce 


#%% Load stuff

# =============================================================================
# Load stuff 
# =============================================================================

# uniqseqs = load_uniqseqs()
# uniqseq_counts = load_uniqseq_map()
# totdist = load_totdist()
SDRs = load_SDRs()
# SDVs = load_SDVs()
SDRnullAll = load_SDRnull()
SDRnull93 = load_SDRnull(gene_set="93genes")
SDRnull47 = load_SDRnull(gene_set="47genes")
SDRnull47super = SDRnull47[SDRnull47['level'] == 'super']
SDRnull47sub = SDRnull47[SDRnull47['level'] == 'sub']

# Merge all SDRnull data
SDRnullTotal = SDRnullAll.append(SDRnull47)
SDRnullTotal = SDRnullTotal.append(SDRnull93)


#%% Process loaded stuff

# Remove duplicate values SDR93

# SDRnull93 = SDRnull93.drop_duplicates(
#   subset = ['gene', 'SDR']).reset_index(drop = True)

# SDRnull47 = SDRnull47.drop_duplicates(
#   subset = ['gene', 'SDR']).reset_index(drop = True)

# SDRnullAll = SDRnullAll.drop_duplicates(
#   subset = ['gene', 'SDR']).reset_index(drop = True)

# SDRnullTotal = SDRnullTotal.drop_duplicates(
#   subset = ['gene', 'SDR']).reset_index(drop = True)

# Conert to numpy
SDRnullTotalSuper = SDRnullTotal[SDRnullTotal['level'] == "nullSuper"]
SDRnullTotalSub = SDRnullTotal[SDRnullTotal['level'] == "nullSub"]

SDRnullTotalSuper.drop(columns = ['level'], inplace = True)
SDRnullTotalSub.drop(columns = ['level'], inplace = True)

# Save 
SDRnullTotalSuper.to_csv('C:/Users/norab/Master/data/SDRnull/all/SDRnullTotalSuper_22.07.21.csv', index = False)
SDRnullTotalSub.to_csv('C:/Users/norab/Master/data/SDRnull/all/SDRnullTotalSub_22.07.21.csv', index = False)



#%% Merge stuff

# =============================================================================
# Merge stuff
# =============================================================================

# SDRnull total + SDR

SDRnull_SDR = SDRnullTotal.append(SDRs)

# data1 = totdist.merge(SDRs, on = 'gene')
# data2 = data1.merge(SDVs, on  = ['gene', 'level']).drop_duplicates()
# #data3 = data2.merge(SDRnullAll, on = ['gene', 'level'], how = 'outer')
# data_all = data2.merge(SDVs, on = ['gene', 'level'], how = 'outer')

# # Wrangle
# data_all.drop_duplicates(inplace=True)
# data_all.dropna(inplace=True)

# # Sort 
# data_all = data_all.sort_values(by=['SDR'])

# # SDR data
# SDRdata = SDRnull.append(SDRs)


# =============================================================================
# Plotting
# =============================================================================

# # Dotplot, lineplot
# (ggplot(data2, aes('SDV', fill = 'level'))
#  + geom_density()
#  + geom_density(aes('SDR'))
#  + theme_classic()
#  + labs(title='SDV')
# )

# (ggplot(data_all, aes('SDR', 'totdist', fill = 'level'))
#  + geom_point()
#  + theme_classic()
#  + labs(title='SDR vs uniqseq')
# )

# (ggplot(data_all, aes('uniqseq', 'totdist', fill = 'level'))
#  + geom_point()
#  + theme_classic()
#  + labs(title='SDR vs uniqseq')
# )

# (ggplot(data_all, aes('SDR', 'SDV', fill = 'level'))
#  + geom_point()
#  + theme_classic()
#  + labs(title='SDR vs uniqseq')
# )

# (ggplot(data_all, aes('SDV', 'uniqseq', fill = 'level'))
#  + geom_point()
#  + theme_classic()
#  + labs(title='SDV vs uniqseq')
#)

# =============================================================================
# Good inspiration: 
#    https://dputhier.github.io/jgb53d-bd-prog_github/practicals/intro_ggplot/intro_ggplot.html
# =============================================================================



# =============================================================================
# Next level violin plots
# =============================================================================


shift = 0.1

def alt_sign(x):
    "Alternate +1/-1 if x is even/odd"
    return (-1) ** x

m1 = aes(x=stage('Level', after_scale='x+shift*alt_sign(x)'))              # shift outward
m2 = aes(x=stage('Level', after_scale='x-shift*alt_sign(x)'), group='gene')  # shift inward


# Violine + points + lines between sub and super
(ggplot(SDRall, aes('Level', 'value', fill = 'Level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.8, show_legend = False)
 + geom_point(m2, color='none', alpha=0.6, size=1.5, show_legend=False)
 + geom_line(m2, color='gray', size=0.65, alpha=0.6)
 + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR values for 8782 genes (all)')
)

# Violin + points + boxplot
(ggplot(SDRall, aes('Level', 'value', fill = 'Level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
 + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
 + scale_fill_manual(values=['dodgerblue', 'darkorange'])
 + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR values for 8782 genes (all)')
)


# =============================================================================
# Violin + points + boxplot
# =============================================================================

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

#%%  Overlapping distribution plot

# TEMPLATE
(ggplot(data="dataset", 
       mapping=aes(x="column_in_dataset", fill='another_column_color')) 
 + geom_density(alpha=0.5)
 + labs(title="title"))


# Null dist, allGenes
(ggplot(data=SDRnullAll, mapping=aes(x='SDR', fill='level')) 
     + geom_density(adjust = 1/4, alpha=0.5)
     + labs(title="SDRnull, all"))

# Null dist, total (all + 47 + 93)
(ggplot(data=SDRnullTotal, mapping=aes(x='SDR', fill='level')) 
     + geom_density(adjust = 1/4, alpha=0.5)
     + labs(title="SDRnull, total"))

# Null dist, 47 genes
(ggplot(data=SDRnull47, mapping=aes(x='SDR', fill='level')) 
     + geom_density(adjust = 1/4, alpha=0.5)
     + labs(title="SDRnull, 47genes"))

# Null dist, 93 genes
(ggplot(data=SDRnull93, mapping=aes(x='SDR', fill='level')) 
     + geom_density(adjust = 1/4, alpha=0.5)
     + labs(title="SDRnull, 93genes"))


# SDR + SDRnull
(ggplot(data=SDRnull_SDR, 
       mapping=aes(x='SDR', fill='level')) 
     + geom_density(adjust = 1/2, alpha=0.5) 
     + scale_fill_manual(values=['indigo', 'darkorange', 'red', 'green'])
     + labs(title='SDRs'))


#%% Overlapping histogram plot


# Total SDRnull values, no duplicates, SDR rounded to 4 decimals
(ggplot(SDRnullTotal) + aes(x="SDR", fill = "level")
     + stat_bin(bins = 1000, alpha = 0.7)
     + scale_fill_manual(values=['indigo', 'darkorange'])
     + labs(title='SDRnull')
     )

# 93 genes, no duplicates, SDR not rounded
(ggplot(SDRnull93) + aes(x="SDR", fill = "level")
     + stat_bin(bins = 1000, alpha = 0.7)
     + scale_fill_manual(values=['indigo', 'darkorange'])
     + labs(title='SDRnull, 93 genes')
     )

# 47 genes, no duplicates, SDR not rounded
(ggplot(SDRnull47) + aes(x="SDR", fill = "level")
     + stat_bin(bins = 1000, alpha = 0.7)
     + scale_fill_manual(values=['indigo', 'darkorange'])
     + labs(title='SDRnull, 47 genes')
     )



# SDRnull total + SDR, no duplicates
(ggplot(SDRnull_SDR) + aes(x="SDR", fill = "level")
      + stat_bin(bins = 1500, alpha = 0.7)
      + scale_fill_manual(values=['indigo', 'darkorange', 'red', 'green'])
      + labs(title='SDRs')
      )


#%% Cumulative distriubtion function plot

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











#%% Trash
