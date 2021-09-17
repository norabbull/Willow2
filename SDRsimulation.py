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
from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine import *  # TO DO: change



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

# Function, different gorup sizes possible
def simSDR(g1, g2, perm, dist = 0.1):
    
    # Within 
    comb_within_g1 = ss.comb(g1, 2, exact = True)
    comb_within_g2 = ss.comb(g2, 2, exact = True)
    comb_within_g1_perm = (g1-perm) * perm 
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


#%%

def simSDR2(g, perm, dist = 0.1):
    # g = number of members in group / population
    # Within 
    comb_within_g = ss.comb(g, 2, exact = True)
    comb_within_g_perm = (g-perm) * perm 
    tot_dist_within = comb_within_g_perm * dist
    mean_within = tot_dist_within / (comb_within_g*2)

    # Between
    comb_between_g_nonperm = g * (g - perm) 
    comb_between_g_g = g * g
    
    tot_dist_bewteen = comb_between_g_nonperm * dist
    mean_between = tot_dist_bewteen / comb_between_g_g
    
    # SDR
    try: 
        SDR = mean_within / mean_between
    except: 
        SDR = "nan"
        print("Division by zero. ")
        print("G1: ", g)
        print("P:", p)
        
    return SDR

#%%  Multiple groups

# Dummies #######
# g = 2000
# num_g = 26         # Number of permutations
# perm = 1000

####################

def simSDR3(g, num_g, perm, dist = 0.1):
    # g = number of members in group / population
    # num_g = nuber of groups
    
    # Within 
    comb_within_g = ss.comb(g, 2, exact = True) * num_g
    
    comb_within_g_perm = (g-perm) * perm 
    tot_dist_within = comb_within_g_perm * dist
    mean_within = tot_dist_within / (comb_within_g)

    # Between
    comb_between_g_nonperm = g * (g - perm) * (num_g - 1)
    comb_between_g_g = (g * g) * ss.comb(num_g, 2, exact = True)
    
    tot_dist_bewteen = comb_between_g_nonperm * dist
    mean_between = tot_dist_bewteen / comb_between_g_g
    
    # SDR
    try: 
        SDR = mean_within / mean_between
    except: 
        SDR = "nan"
        print("Division by zero. ")
        print("G1: ", g)
        print("P:", p)
        
    return SDR

#%% Run simulation

SDRresults = pd.DataFrame(columns = ['Group1', 'Group2', 'Perm', 'SDR'])
#SDRresults2 = pd.DataFrame(columns = ['Group1', 'Group2', 'Perm', 'SDR'])
#for g1 in range(2,50):
for g in range(2,100):
    p = 1
    while p < (g/2):
        SDR = simSDR2(g, g, p)
        SDRresults.loc[len(SDRresults)] = [g, g, p, SDR]
        p += 1

#%% 
SDRresults2 = pd.DataFrame(columns = ['Group1', 'Group2', 'Perm', 'SDR'])
#SDRresults2 = pd.DataFrame(columns = ['Group1', 'Group2', 'Perm', 'SDR'])
#for g1 in range(2,50):
for g2 in range(2,100):
    p = 1
    while p < (g2/2):
        SDR = simSDR2(g, p)
        SDRresults2.loc[len(SDRresults2)] = [g, g, p, SDR]
        p += 1
            
#%% Testing sim3

SDRresults3 = pd.DataFrame(columns = ['Group_size', 'Num_groups', 'Perm', 'SDR'])
#for g1 in range(2,50):
for g in range(2,100):
    for ng in range(2,3):
        p = 0
        while p < (g):
            SDR = simSDR3(g, num_g = ng, perm = p)
            SDRresults3.loc[len(SDRresults3)] = [g, ng, p, SDR]
            p += 1
            
#%% Single val calc

g = 10
p = 8
num_g = 6
SDR = simSDR3(g,num_g,perm = p)


#%% Plot


(ggplot(SDRresults3, aes('SDR', 'Group_size', fill = 'Perm'))
 + geom_point()
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
# + scale_fill_manual(values=['dodgerblue', 'darkorange'])
 + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR simulation')
)














