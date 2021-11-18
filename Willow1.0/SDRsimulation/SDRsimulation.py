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
#from inspection.inspectionHelpers import *
import random




#%% Visualize  non random

(ggplot(GDRresults5, aes('G2Perm', 'Group_size', fill = 'GDR'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
# + scale_fill_manual(values=['dodgerblue', 'darkorange'])
# + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='GDR simulation')
)





#%% SAVED VERSION TO PRODUCE SIMULATION PLOT: 
    
    
def simGDR(G, G2perm, Random = False, dist = 0.1):
    """
    Simulation procedure: 
        Simulation of GDR values in a tree comprising 2 clades, C1 and C2,
        separated by a distance (dist).
        The leaf nodes of the tree is defined to belong to either of two 
        groups, G1 and G2.
          
        G = group sizes, G = G1 = G2
        C1 = number of samples in clade 1
        C2 = number of samples in clade 2
        G2perm = number of G2 samples permuted from C2 to C1
        
        When number of G2perm = 0, G1 = G2 = C1 = C2.
        G1 and G2 are perfectly separated into the two clades; 
        all G1 samples are in C1, all G2 samples are in C2. 
        
        As G2perm increase, G2 samples are "moved" from C2 to C1.
        
        Number of possible permutations ranges from 0 to G; in which 
        case all samples are clustered together in the same clade, and C2
        no longer exist.
        
        GDR is calculated as a function of group size and number 
        of permutations from G2 to G1. 
        
    Random simulation: 
        Obtain GDR value when partition of G1 and G2 samples into 
        C1 and C2 is random. Sizes of C1 and C2 are determined by the 
        same criteria as above (number of permutations).
        This will simulate the same situation as if leaf nodes are 
        randomly assigned to G1 and G2.
        
        Result: X number of GDR values, each obtained from a random 
        tree clustering
        
    """
    # Clade and group sizes    
    C1 = G + G2perm       # num samples in clade 1
    C2 = G - G2perm       # num samples in clade 2
   
    if Random: 
        G1C2 = random.randint(0,C2)   # G1C2 = num group 1 samples in clade 2
    else:
        G1C2 = 0      

    G2C2 = C2 - G1C2      # G2C2 = num group 2 samples in clade 2
    G1C1 = G - G1C2       # G1C1 = num group 1 samples in clade 1
    G2C1 = G - G2C2       # G2C1 = num group 2 samples in clade 1      
    
    # Total number of pairwise distances between samples of same group
    comb_within_G = ss.comb(G, 2, exact = True) * 2
    
    # Distance between samples from same group with dist != 0:
    dist_within_G1C1_G1C2 = (G1C1 * G1C2) * dist
    dist_within_G2C1_G2C2 = (G2C1 * G2C2) * dist

    tot_dist_within = dist_within_G1C1_G1C2 + dist_within_G2C1_G2C2
    mean_within = tot_dist_within / comb_within_G

    # Total number of pairwise distances between G1 and G2 samples
    comb_between_G1_G2 = G * G
    
    # Distance between samples from different gorups with dist != 0
    dist_between_G1C1_G2C2 = (G1C1 * G2C2) * dist
    dist_between_G1C2_G2C1 = (G1C2 * G2C1) * dist
    
    tot_dist_between = dist_between_G1C1_G2C2 + dist_between_G1C2_G2C1
    mean_between = tot_dist_between / comb_between_G1_G2
    
    # GDR 
    GDR = mean_within / mean_between

    return GDR



#%% Randomize

GDRresults = pd.DataFrame(columns = ['Group_size', 'Permutations', 'GDR'])

for G in range(2,200):
    p = 0
    while p < G:
        GDR = simGDR(G, G2perm = p, Random = False)
        GDRresults.loc[len(GDRresults)] = [G, p, GDR]
        p += 1


#%% Visualize simGDR (not random)

(ggplot(GDRresults, aes('Permutations', 'Group_size', fill = 'GDR'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
 + theme(figure_size=(8, 6))
 + labs(title='GDR simulation')
)


#%% Randomize

GDRresultsRAND = pd.DataFrame(columns = ['Group_size', 'Permutations', 'GDR', 'randomGDR'])

for G in range(2,200):
    p = 0
    while p < G:
        GDR = simGDR(G, G2perm = p, Random = False)
        randomGDRs = []
        for _ in range(100):
            rGDR = simGDR(G, G2perm = p, Random = True)
            randomGDRs.append(rGDR)
            
        GDRresultsRAND.loc[len(GDRresultsRAND)] = [G, p, GDR, randomGDRs]
        p += 1

#%% Visualize simGDR (not random)

(ggplot(GDRresultsRAND, aes('Permutations', 'Group_size', fill = 'randomGDR'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
 + theme(figure_size=(8, 6))
 + labs(title='GDR simulation')
 + theme(text=element_text(size=15))
)

#%% Save values to files

for index, row in GDRresultsRAND.iterrows():
    
    header = 'G:' + str(row['Group_size']) + '_P:' + str(row['Permutations']) + '_GDR:' + str(row['GDR'])
    print(index)

    file = 'C:/Users/norab/Master/thesis_data/simulation/' + str(header) + '.csv'
    with open(file, 'w') as f:
        f.write(header)
        f.write('\n')
        for item in row['randomGDR']:
            f.write("%s\n" % item)



#%% Visualize random


all_simData = load_simData()

(ggplot(all_simData, aes('perm', 'group_size', fill = 'pval'))
 + geom_point(alpha=1, size=2, stroke = 0.1, color = 'indigo')
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
# + scale_fill_manual(values=['dodgerblue', 'darkorange'])
# + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='GDR simulation')
)


(ggplot(all_simData, aes('perm', 'group_size', fill = 'pval_adj'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
# + scale_fill_manual(values=['dodgerblue', 'darkorange'])
# + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='GDR simulation')
)


