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
from inspection.inspectionHelpers import *

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
for g in range(2,200):
    for ng in range(2,3):
        p = 0
        while p < (g):
            SDR = simSDR3(g, num_g = ng, perm = p, dist = 500)
            SDRresults3.loc[len(SDRresults3)] = [g, ng, p, SDR]
            p += 1
            
#%% Single val calc

g = 10
p = 8
num_g = 6
SDR = simSDR3(g,num_g,perm = p)


#%% Plot


(ggplot(SDRresults3, aes('SDR', 'Group_size', size = 'Perm', color = 'Perm'))
 + geom_point(alpha=1, size=2,stroke = 0.1, color = 'indigo')
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
 #+ scale_fill_manual(values=['dodgerblue', 'darkorange'])
# + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR simulation')
)

(ggplot(SDRresults3, aes('Perm', 'Group_size', fill = 'SDR'))
 + geom_point(alpha=1, size=3,stroke = 0.1, color = 'indigo')
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
 #+ scale_fill_manual(values=['dodgerblue', 'darkorange'])
# + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR simulation')
)


# The simulation demonstartes how group size affects the SDR value to a large degree. 
# So reading an SDR value by itself is not necessarily very meaningful, it must be seen in realtion
# to the group sizes. 
# EG. if You only have 25 samples in total, a few permutations from one calde to another will affect the SDR value alot. 
# But if you have A group sample of 100, one permutation fourt or back wont make the big difference. 
# But this is likely to also turn up in the statistical testing as well.
# I will simulate p values, but hasn't gottne that far. 
# So what this really demonstrates is the well established cocept that sample size is important.


#%% Calculate p - values

"""
    1. Construct distance matrices for 100 hypothetical trees
    with increasing number of permutations. 
    2. Create sample labels that groups them into the 2 separate clades. 
    
    GR1___SIM___IamLeafPURPLE
    GR2___SIM___IamLeafGREEN
    
"""
# Construct names

# Constants
sample_size = 200
num_trees = 100

# Sample labels
samples = ['GR1___SIM___IamLeafPurple' for i in range(int(sample_size/2))]
samples2 = ['GR2___SIM___IamLeafGreen' for i in range(int(sample_size/2))]
samples.extend(samples2)

# Construct data

for tree_num in range(num_trees):
    # loop to make multiple matices
    data = []
    # for r in range(sample_size):
    #     #row = []
    for row in range(tree_num):
        ones = [1 for _ in range(sample_size)]
        data.append(ones)
    for row in range(sample_size - tree_num):
        ones = [1 for _ in range(tree_num)]
        zero = [0 for _ in range(sample_size - len(ones))]
        ones.extend(zero)
        
        data.append(ones)

    for i in range(sample_size):
        data[i][i] = 0

    data = pd.DataFrame(data, columns = [samples], index =  [samples])
    data.to_csv('C:/Users/norab/Master/data/simulation/sim_mat{0}.csv'.format(tree_num), index = True)


# Check mat

sim_mat55 = load_cd_mat('C:/Users/norab/Master/data/job_simulation/cd_simulation/sim_mat55.csv')

#%% Calculate SDR for all sim mats








