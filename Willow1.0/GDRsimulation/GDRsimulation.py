# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 13:47:30 2021

@author: norab
"""



import scipy.special as ss
import pandas as pd
from plotnine import ggplot, aes, theme, geom_point, labs
from plotnine import *
import random

import scipy
scipy.__version__
import plotnine
plotnine.__version__

"""
Function to calculate simulated GDR values 
"""
    
def simGDR(G, p, Random = False, dist = 0.1):
    """
    Simulation procedure: 
        Simulation of GDR values in a tree comprising 2 clades, C1 and C2,
        separated by a distance (dist).
        The leaf nodes of the tree is defined to belong to either of two 
        groups, G1 and G2.
          
        G = group sizes, G = G1 = G2
        C1 = number of samples in clade 1
        C2 = number of samples in clade 2
        p = number of G2 samples permuted from C2 to C1
        
        When number of p = 0, G1 = G2 = C1 = C2.
        G1 and G2 are perfectly separated into the two clades; 
        all G1 samples are in C1, all G2 samples are in C2. 
        
        As p increase, G2 samples are "moved" from C2 to C1.
        
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
    # C1 = G + p       # num samples in clade 1
    C2 = G - p       # num samples in clade 2
   
    if Random: 
        G1C2 = random.randint(0,C2)   # G1C2 = num group 1 samples in clade 2
    else:
        G1C2 = 0      

    G2C2 = C2 - G1C2      # G2C2 = num group 2 samples in clade 2
    G1C1 = G - G1C2       # G1C1 = num group 1 samples in clade 1
    G2C1 = G - G2C2       # G2C1 = num group 2 samples in clade 1      
    
    # Total number of pairwise distances between samples of same group
    comb_within_G = ss.comb(G, 2, exact = True) * 2  # 2 groups of equal sizes
    
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



#%% Old

GDRresults = pd.DataFrame(columns = ['Group_size', 'Permutations', 'GDR'])

for G in range(2,200):
    p = 0
    while p < G:
        GDR = simGDR(G, p = p, Random = False)
        GDRresults.loc[len(GDRresults)] = [G, p, GDR]
        p += 1




#%% Calculate GDRs + random GDRs

"""
GDR simulation: 
    - Calculates simulated GDR values in with simGDR function for
    group size, G, in range 2-200 and all permutations, p, where p < G.

    - Calculates 100 random GDR values for each constellation.

"""

# Table to save results
GDRresults = pd.DataFrame(columns = ['Group_size', 'Permutations', 'GDR', 'randomGDR'])

for G in range(2,200):              # group size in range 2 - 200
    p = 0                           # p = number of permutations
    while p < G:                
        GDR = simGDR(G, p = p, Random = False)      # simGDR calculation
        randomGDRs = []
        for _ in range(100):
            rGDR = simGDR(G, p = p, Random = True)  # randomGDR calculation 
            randomGDRs.append(rGDR)
            
        GDRresults.loc[len(GDRresults)] = [G, p, GDR, randomGDRs]
        p += 1


#%% Visualize simGDR

"""
Plot: 
    Simulated GDR values 
    
"""

(ggplot(GDRresults, aes('Permutations', 'Group_size', fill = 'GDR'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
 + theme(figure_size=(7, 7))
 + labs(title='GDR simulation')
 + scale_x_continuous(name="Number of G2-samples moved from C2 to C1")
 + scale_y_continuous(name="Group size  (number of samples in G)")
 + scale_color_cmap(cmap_name="inferno")
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

#%% MAL


(ggplot(seqdata_HAUS4) 
 + geom_bar(aes(x='seq', fill='super')) 
 + labs(title = "HAUS4, SUPER     GDR = 0.971") 
 #+ scale_x_discrete(name="unique samples")
 #+ scale_y_discrete(name="count")
 + xlab("unique samples")
 + theme(figure_size=(13, 13))
 + theme(axis_text_x=element_text(rotation=30, hjust=1),
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
#%% Save random values to files

"""
Save values to file

"""
for index, row in GDRresults.iterrows():
    
    title = 'G' + str(row['Group_size']) + '_P' + str(row['Permutations']) 
    header = title + '_GDR' + str(row['GDR'])
    print(index)

    file = 'C:/Users/norab/Master/thesis_data/simulation/simData/randSim' + title + '.csv'
    with open(file, 'w') as f:
        f.write(header)
        f.write('\n')
        for item in row['randomGDR']:
            f.write("%s\n" % item)


#%% 



#%% Visualize random

"""
p-value calculation:
    p-values for all random GDRs are calculated in R, see script 
    "GDRsimulationStats".
    p-values are read and visualized with code below.
    
"""

# Read data
pval_file = 'C:/Users/norab/Master/thesis_data/simulation/simNull_alldata_24.11.21.csv'
all_simData = pd.read_csv(pval_file)


#%%
all_simData.rename(columns={'pval':'p-val'}, inplace=True)

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


#%% OLD"
(ggplot(all_simData, aes('permutations', 'group_size', fill = 'p-val'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
 + theme(figure_size=(13, 13))
 + labs(title='GDR simulation: p - values')
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



#---------------------------------------------------
# plot p-values, uncorrected
(ggplot(all_simData, aes('permutations', 'group_size', fill = 'pval'))
 + geom_point(alpha=1, size=3.2, stroke = 0.1, color = 'indigo')
 + theme(figure_size=(8, 6))
 + labs(title='GDR simulation p-values')
)

#%% 

all_simData.rename(columns={'pval_adj_holm':'adj p-val'}, inplace=True)

(ggplot(all_simData, aes('permutations', 'group_size', fill = 'adj p-val'))
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

#---------------------------------------------

# plot p-values, corrected
(ggplot(all_simData, aes('permutations', 'group_size', fill = 'pval_adj_holm'))
 + geom_point(alpha=1, size=3.2, stroke = 0.1, color = 'indigo')
 + theme(figure_size=(8, 6))
 + labs(title='GDR simulation, adjusted p-values')
)




