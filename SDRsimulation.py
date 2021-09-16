# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 13:47:30 2021

@author: norab
"""





import numpy as np
import math
import scipy.special as ss
import itertools
import pandas as pd
import matplotlib.pyplot as plt

#%% 


# # Testing #######
# ss.comb(4, 2, exact = True)
# itertools.permutations(3, 2)

# # Dummies #######
# g1 = 5
# g2 = 10
# dist = 0.1
# perm  = 3         # Number of permutations

# ####################

def simSDR(g1, g2, perm, dist = 0.1):
    
    # Within 
    comb_within_g1 = ss.comb(g1, 2, exact = True)
    comb_within_g2 = ss.comb(g2, 2, exact = True)
    comb_within_g1_perm = g1 * perm 
    tot_dist_within = comb_within_g1_perm * dist
    mean_within = tot_dist_within / (comb_within_g1 + comb_within_g2)

    # Between
    comb_between_g1_nonperm = g1 * (g2 - perm) 
    comb_between_g2_g1 = g1 * g2 
    
    tot_dist_bewteen = comb_between_g1_nonperm * dist
    mean_between = tot_dist_bewteen / comb_between_g2_g1
    
    # SDR
    try: 
        SDR = mean_within / mean_between
    except: 
        SDR = "nan"
        print("Division by zero. ")
        print("G1: ", g1)
        print("G2: ", g2)
        print("P:", p)
        
    return SDR

#%% Run simulation

group1 = list(np.arange(2,100))
group2 = list(np.arange(2,100))

SDRresults = pd.DataFrame(columns = ['Group1', 'Group2', 'Perm', 'SDR'])
SDRresults2 = pd.DataFrame(columns = ['Group1', 'Group2', 'Perm', 'SDR'])
#for g1 in range(2,50):
for g2 in range(2,100):
    p = 1
    while p < (g2/2):
        SDR = simSDR(g2, g2, p)
        SDRresults2.loc[len(SDRresults2)] = [g2, g2, p, SDR]
        p += 1
            
#%% Plot



from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine import *  # TO DO: change


(ggplot(SDRresults2, aes('SDR', 'Perm', fill = 'Group2'))
 + geom_point()
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
# + scale_fill_manual(values=['dodgerblue', 'darkorange'])
 + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR simulation')
)














