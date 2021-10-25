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
import random

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
# def simSDR(g1, g2, perm, dist = 0.1):
    
#     # Within 
#     comb_within_g1 = ss.comb(g1, 2, exact = True)
#     comb_within_g2 = ss.comb(g2, 2, exact = True)
#     comb_within_g1_perm = (g1-perm) * perm 
#     tot_dist_within = comb_within_g1_perm * dist
#     mean_within = tot_dist_within / (comb_within_g1 + comb_within_g2)

#     # Between
#     comb_between_g1_nonperm = g1 * (g2 - perm) 
#     comb_between_g2_g1 = g1 * g2 
    
#     tot_dist_bewteen = comb_between_g1_nonperm * dist
#     mean_between = tot_dist_bewteen / comb_between_g2_g1
    
#     # SDR
#     try: 
#         SDR = mean_within / mean_between
#     except: 
#         SDR = "nan"
#         print("Division by zero. ")
#         print("G1: ", g1)
#         print("G2: ", g2)
#         print("P:", p)
        
#     return SDR


#%%

# def simSDR2(g, perm, dist = 0.1):
#     # g = number of members in group / population
#     # Within 
#     comb_within_g = ss.comb(g, 2, exact = True)
#     comb_within_g_perm = (g-perm) * perm 
#     tot_dist_within = comb_within_g_perm * dist
#     mean_within = tot_dist_within / (comb_within_g*2)

#     # Between
#     comb_between_g_nonperm = g * (g - perm) 
#     comb_between_g_g = g * g
    
#     tot_dist_bewteen = comb_between_g_nonperm * dist
#     mean_between = tot_dist_bewteen / comb_between_g_g
    
#     # SDR
#     try: 
#         SDR = mean_within / mean_between
#     except: 
#         SDR = "nan"
#         print("Division by zero. ")
#         print("G1: ", g)
#         print("P:", p)
        
#     return SDR

#%%  Multiple groups

# Dummies #######
# g = 2000
# num_g = 26         # Number of permutations
# perm = 1000

####################

# def simSDR3(g, num_g, perm, dist = 0.1):
#     # g = number of members in group / population
#     # num_g = nuber of groups
    
#     # Within 
#     comb_within_g = ss.comb(g, 2, exact = True) * num_g
    
#     comb_within_g_perm = (g-perm) * perm 
#     tot_dist_within = comb_within_g_perm * dist
#     mean_within = tot_dist_within / (comb_within_g)

#     # Between
#     comb_between_g_nonperm = g * (g - perm) * (num_g - 1)
#     comb_between_g_g = (g * g) * ss.comb(num_g, 2, exact = True)
    
#     tot_dist_bewteen = comb_between_g_nonperm * dist
#     mean_between = tot_dist_bewteen / comb_between_g_g
    
#     # SDR
#     try: 
#         SDR = mean_within / mean_between
#     except: 
#         SDR = "nan"
#         print("Division by zero. ")
#         print("G1: ", g)
#         print("P:", p)
        
#     return SDR

#%% simSDR4 - without number of groups (unnecessary, since only two clades.)


# def simSDR4(G, G2perm, dist = 0.1):
#     """
#         Simulation of SDR values in a tree comprising 2 clades, C1 and C2,
#         separated by a distance (dist).
#         The leaf nodes of the tree is defined to belong to either of two 
#         groups, G1 and G2.
          
#         G = group size of either group. G1 = G2, always.
#         C1 = number of samples in clade 1
#         C2 = number of samples in clade 2
        
#         At start of simulation, G1 and G2 are perfectly separated
#         into the two clades; G1 samples are in C1, G2 samples are in C2. 
#         Then, G1 = G2 = C1 = C2.
        
#         G2perm = number of G2 samples permuted from C2 to C1.
        
#         When number of G2 permutations (G2perm) = 0, the groups are perfectly 
#         separated into C1 and C2, but as G2perm increase, G2 samples are 
#         clustered together with G1 samples in C1, rather than separately in C2. 
        
#         Number of possible permutations ranges from 0 to G; in which 
#         case all samples are clustered together in the same clade, and C2
#         no longer exist.
        
#         SDR is calculated as a function of G and perm; group size 
#         and number of permutations from G2 to G1. 
          
#     """
    
#     C1 = G + P      # num samples in tree clade 1
#     C2 = G - P      # num samples in tree clade 2
    
      
#     # Within 
#     comb_within_G = ss.comb(G, 2, exact = True) * 2   # Total nr. of combinations between samples of the same defined group
#     comb_within_G2_G2perm = C2 * G2perm   # Nr. of sample combinations from same group (G1) with dist != 0
#     tot_dist_within = comb_within_G2_G2perm * dist
#     mean_within = tot_dist_within / comb_within_G

#     # Between
#     comb_between_G1_G2 = G * G
#     comb_between_G1_G2nonperm = G * C2   # Nr. of sample combinations from different groups (G1 and G2) with dist != 0

#     tot_dist_bewteen = comb_between_G1_G2nonperm * dist
#     mean_between = tot_dist_bewteen / comb_between_G1_G2
    
#     # SDR
#     try: 
#         SDR = mean_within / mean_between
#     except: 
#         SDR = "nan"
#         print("Division by zero. ")
#         print("G1: ", g)
#         print("P:", p)
        
#     return SDR

#%% Run simulation

# SDRresults = pd.DataFrame(columns = ['Group1', 'Group2', 'Perm', 'SDR'])
# #SDRresults2 = pd.DataFrame(columns = ['Group1', 'Group2', 'Perm', 'SDR'])
# #for g1 in range(2,50):
# for g in range(2,100):
#     p = 1
#     while p < (g/2):
#         SDR = simSDR2(g, g, p)
#         SDRresults.loc[len(SDRresults)] = [g, g, p, SDR]
#         p += 1

# #%% 
# SDRresults2 = pd.DataFrame(columns = ['Group1', 'Group2', 'Perm', 'SDR'])
# #SDRresults2 = pd.DataFrame(columns = ['Group1', 'Group2', 'Perm', 'SDR'])
# #for g1 in range(2,50):
# for g2 in range(2,100):
#     p = 1
#     while p < (g2/2):
#         SDR = simSDR2(g, p)
#         SDRresults2.loc[len(SDRresults2)] = [g, g, p, SDR]
#         p += 1
            
#%% Testing sim3

# SDRresults3 = pd.DataFrame(columns = ['Group_size', 'Num_groups', 'Perm', 'SDR'])
# #for g1 in range(2,50):
# for g in range(2,200):
#     for ng in range(2,3):
#         p = 0
        
#         while p < g:
#             SDR = simSDR3(g, num_g = ng, perm = p, dist = 500)
#             SDRresults3.loc[len(SDRresults3)] = [g, ng, p, SDR]
#             p += 1

#%% 

# SDRresults4 = pd.DataFrame(columns = ['Group_size', 'Num_groups', 'Perm', 'SDR'])
# #for g1 in range(2,50):
# for g in range(2,200):
#         p = 0
#         while p < g:
#             SDR = simSDR4(g, perm = p, dist = 500)
#             SDRresults4.loc[len(SDRresults4)] = [g, ng, p, SDR]
#             p += 1



#%% 


def simSDR_random(G, G2perm, Random = False, dist = 0.1):
    """
        Simulation of SDR values in a tree comprising 2 clades, C1 and C2,
        separated by a distance (dist).
        The leaf nodes of the tree is defined to belong to either of two 
        groups, G1 and G2.
          
        G = group size of either group. G1 = G2, always.
        C1 = number of samples in clade 1
        C2 = number of samples in clade 2
        
        At start of simulation, G1 and G2 are perfectly separated
        into the two clades; G1 samples are in C1, G2 samples are in C2. 
        Then, G1 = G2 = C1 = C2.
        
        G2perm = number of G2 samples permuted from C2 to C1.
        
        When number of G2 permutations (G2perm) = 0, the groups are perfectly 
        separated into C1 and C2, but as G2perm increase, G2 samples are 
        clustered together with G1 samples in C1, rather than separately in C2. 
        
        Number of possible permutations ranges from 0 to G; in which 
        case all samples are clustered together in the same clade, and C2
        no longer exist.
        
        SDR is calculated as a function of G and perm; group size 
        and number of permutations from G2 to G1. 
        
        Random: 
            Obtain SDR value when partition of G1 and G2 samples into 
            C1 and C2 is random. Sizes of C1 and C2 are determined by the 
            same criteria as above (number of permutations).
            (This will simulate the same situation as if leaf nodes are 
            randomly assigned to G1 and G2.
            
           ((( Result: #X number of SDR values, each obtained from a random 
            tree clustering )))
            
    """
    
    C1 = G + G2perm      # num samples in tree clade 1
    C2 = G - G2perm      # num samples in tree clade 2
    
    if Random: 
        G1C2 = random.randint(0,C2)     # G1C2 = Group 1 samples in clade 2
    else:
        G1C2 = 0
        
    G2C2 = C2 - G1C2                # G2C2 = Group 2 samples in clade 2
    G1C1 = G - G1C2                 # G1C1 = Group 1 samples in clade 1
    G2C1 = G - G2C2                 # G2C1 = Group 2 samples in clade 1
        
    # Total number of pairwise distances between samples of same group
    comb_within_G = ss.comb(G, 2, exact = True) * 2   # Total nr. of combinations between samples of the same defined group
    
    # Distance betweem samples from same group with dist != 0:
    dist_within_G1C1_G1C2 = (G1C1 * G1C2) * dist
    dist_within_G2C1_G2C2 = (G2C1 * G2C2) * dist

    tot_dist_within = dist_within_G1C1_G1C2 + dist_within_G2C1_G2C2
    mean_within = tot_dist_within / comb_within_G


    # Total number of pairwise distances between samples of different groups
    comb_between_G1_G2 = G * G
    
    # Distance between samples form different gorups with dist != 0
    dist_between_G1C1_G2C2 = (G1C1 * G2C2) * dist
    dist_between_G1C2_G2C1 = (G1C2 * G2C1) * dist
    
    tot_dist_between = dist_between_G1C1_G2C2 + dist_between_G1C2_G2C1
    
    mean_between = tot_dist_between / comb_between_G1_G2
    
    # SDR
    try: 
        SDR = mean_within / mean_between
    except: 
        SDR = float("nan")
        print("Division by zero. ")
        print("G1: ", G)
        print("P:", G2perm)
    

    return SDR




#%% Non random
SDRresults5 = pd.DataFrame(columns = ['Group_size', 'G2Perm', 'SDR'])
#for g1 in range(2,50):
for G in range(2,200):
        p = 0
        while p < G:
            SDR = simSDR_random(G, G2perm = p, Random = False)
            SDRresults5.loc[len(SDRresults5)] = [G, p, SDR]
            p += 1


#%% Visualize  non random

(ggplot(SDRresults5, aes('G2Perm', 'Group_size', fill = 'SDR'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
# + scale_fill_manual(values=['dodgerblue', 'darkorange'])
# + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR simulation')
)


#%% Randomize

SDRresults6 = pd.DataFrame(columns = ['Group_size', 'G2Perm', 'SDR', 'randomSDRs'])
#for g1 in range(2,50):
for G in range(2,200):
        p = 0
        while p < G:
            SDR = simSDR_random(G, G2perm = p, Random = False)
            randomSDRs = []
            for _ in range(100):
                rSDR = simSDR_random(G, G2perm = p, Random = True)
                randomSDRs.append(rSDR)
                
            SDRresults6.loc[len(SDRresults6)] = [G, p, SDR, randomSDRs]
            p += 1

#%% Save values to files

for index, row in SDRresults6.iterrows():
    
    header = 'G:' + str(row['Group_size']) + '_P:' + str(row['G2Perm']) + '_SDR:' + str(row['SDR'])
    print(index)

    file = 'C:/Users/norab/Master/data/simulation/simNull/simNull_' + str(index) + '.csv'
    with open(file, 'w') as f:
        f.write(header)
        f.write('\n')
        for item in row['randomSDRs']:
            f.write("%s\n" % item)



#%% Visualize random! 

# Filled Density Plot
ggplot(SDRresults6['randomSDRs'], aes(x = randomSDRs)) + 
       geom_histogram(color = "red", # Curve color
                    fill = "purple",  # Area color
                    alpha = 0.5)   # Area transparency


#%% 


# N = population size
# G = number of samples in groups (G1 = G2)
# n = size of C2 (sample size)
# x = number of G1 elements in C2
# p = number of permutations from G1 to C1

# x = n because the simulation only includes situation of 
# C2 beaing a "clean" G1 clade. 


#____________________________
# M = total number of objects
# n = total number of Group 1 objects
# N = number of samples in random draw. In my case, samples in clade 2. 

# M = Total number of leaf nodes 
# n = total number of G1 objects. G1 = G2.
# N = number of samples in clade 2. ("outgroup")
# x = number of G1 objects in N, drawn without replacement
#     from the total population

# # N = x

# SDRresults4 = pd.DataFrame(columns = ['Group_size', 'Perm', 'SDR', 'p-value', 'x', 'M', 'N'])
# #for g1 in range(2,50):

# for G in range(2,200):
#     p = 0
#     N = G - p
#     n = G
#     M = G * 2
#     x = N
#     while p < G:
#         SDR = simSDR4(G, perm = p)
#         pval = hypergeom.cdf(x, M, n, N)  # (x, M, n, N)
        
#         SDRresults4.loc[len(SDRresults4)] = [G, p, SDR, pval, x, M, N]
#         p += 1


#%% Single val calc

# g = 10
# p = 8
# num_g = 6
# SDR = simSDR3(g,num_g,perm = p)


#%% Plot


# (ggplot(SDRresults3, aes('SDR', 'Group_size', size = 'Perm', color = 'Perm'))
#  + geom_point(alpha=1, size=2,stroke = 0.1, color = 'indigo')
# # + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# # + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
#  #+ scale_fill_manual(values=['dodgerblue', 'darkorange'])
# # + theme_classic()
#  + theme(figure_size=(8, 6))
#  + labs(title='SDR simulation')
# )

# (ggplot(SDRresults3, aes('Perm', 'Group_size', fill = 'SDR'))
#  + geom_point(alpha=1, size=3,stroke = 0.1, color = 'indigo')
# # + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# # + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
#  #+ scale_fill_manual(values=['dodgerblue', 'darkorange'])
# # + theme_classic()
#  + theme(figure_size=(8, 6))
#  + labs(title='SDR simulation')
# )


# (ggplot(SDRresults4, aes('Perm', 'Group_size', fill = 'SDR'))
#  + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
# # + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# # + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
#  #+ scale_fill_manual(values=['dodgerblue', 'darkorange'])
# # + theme_classic()
#  + theme(figure_size=(8, 6))
#  + labs(title='SDR simulation')
# )

(ggplot(SDRresults5, aes('G2Perm', 'Group_size', fill = 'SDR'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
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

# """
#     1. Construct distance matrices for 100 hypothetical trees
#     with increasing number of permutations. 
#     2. Create sample labels that groups them into the 2 separate clades. 
    
#     GR1___SIM___IamLeafPURPLE
#     GR2___SIM___IamLeafGREEN
    
# """


# # Constants
# sample_size = 200
# num_trees = 200

# # Sample labels
# samples = ['GR1___SIM___IamLeafPurple' for i in range(int(sample_size/2))]
# samples2 = ['GR2___SIM___IamLeafGreen' for i in range(int(sample_size/2))]
# samples.extend(samples2)

# # Construct data

# for tree_num in range(num_trees):
#     # loop to make multiple matices
#     data = []
#     # for r in range(sample_size):
#     #     #row = []
#     for row in range(tree_num):
#         ones = [1 for _ in range(sample_size)]
#         data.append(ones)
#     for row in range(sample_size - tree_num):
#         ones = [1 for _ in range(tree_num)]
#         zero = [0 for _ in range(sample_size - len(ones))]
#         ones.extend(zero)
        
#         data.append(ones)

#     for i in range(sample_size):
#         data[i][i] = 0

#     data = pd.DataFrame(data, columns = [samples], index =  [samples])
#     data.to_csv('C:/Users/norab/Master/data/simulation/sim_mat{0}.csv'.format(tree_num), index = True)


# # Check mat

# sim_mat55 = load_cd_mat('C:/Users/norab/Master/data/job_simulation/cd_simulation/sim_mat55.csv')
# sim_mat100 = load_cd_mat('C:/Users/norab/Master/data/job_simulation/cd_simulation/sim_mat100.csv')

# sim_mat199 = load_cd_mat('C:/Users/norab/Master/data/job_simulation/cd_simulation/sim_mat199.csv')

# # #%% Second attempt

# # # Constants
# # sample_size = 200
# # num_trees = 200

# # # Sample labels
# # samples = ['GR1___SIM___IamLeafPurple' for i in range(int(sample_size/2))]
# # samples2 = ['GR2___SIM___IamLeafGreen' for i in range(int(sample_size/2))]
# # samples.extend(samples2)

# # # Construct data

# # data = pd.DataFrame(np.zeros((sample_size, sample_size)) )
# # dist_mat = self.dist_mat.to_numpy()
        
# #         # Create distance sum placeholders
# #         # Dict with pop as key and [dist value, number of dist values added] as value
# #         WithSums = {subtype : [0,0] for subtype in self.subtype_info[0]}
# #         BetSums = {subtype : [0,0] for subtype in self.subtype_info[0]}
        

# #         # Iter upper triangular dist_mat
# #         row_length = len(dist_mat)
# #         row_start = 1
# #         for col in range(row_length):
# #             for row in range(row_start, row_length):
# #                 sample1 = col
# #                 sample2 = row
                
# #                 # Get popName
# #                 Subtype1 = self.sample_info[sample1][1]
# #                 Subtype2 = self.sample_info[sample2][1]
                
# #                 val = dist_mat[sample2][sample1]        # Distance value

# #                 if Subtype1 == Subtype2:        # Same super pop, same sub pop
# #                     WithSums[Subtype1] = [i+j for i, j in zip(WithSums[Subtype1], [val, 1])]  
# #                 else:
# #                     BetSums[Subtype1] = [i+j for i, j in zip(BetSums[Subtype1], [val, 1])]  
# #                     BetSums[Subtype2] = [i+j for i, j in zip(BetSums[Subtype2], [val, 1])]  
                        
# #             row_start += 1
            
# #         #self.pop_dists = {'With': WithSums, 'Bet': BetSums}
# #         self.subtype_dists = {'With': WithSums, 'Bet': BetSums}


#%% Hypergeometric probability function


# from scipy.stats import hypergeom
# import matplotlib.pyplot as plt



# # M = total number of objects
# # n = total number of Group 1 objects
# # N = number of samples in random draw. In my case, samples in clade 2. 

# # M = Total number of leaf nodes 
# # n = total number of G1 objects. G1 = G2.
# # N = number of samples in clade 2. ("outgroup")
# # x = number of G1 objects in N drawn without replacement
# #     from the total population

# """ 

#     Thought process: DIscussion!
#         Calculate the probability that the samples in 
#         one of the clades randomly appear in this clade, 
#         if there was no underlying pattern creating this 
#         group structure. 
#         Really, it is a calculation of the samples being
#         drawn from the population. 
#         Which would be the case if the underlying cause of
#         the observed variation in the sequence data were 
#         completely random. 
#         As if the mutations observed was there for no reason. 
#         ... more on drive. 
# """

# # Probability mass function, pmf
# [M, n, N] = [10, 5, 5]
# rv = hypergeom(M, n, N)
# x = np.arange(0, n+1)
# pmf_g1 = rv.pmf(x)
# cdf_g1 = rv.cdf(x)

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(x, pmf_g1, 'bo')
# ax.vlines(x, 0, pmf_g1, lw=2)
# ax.set_xlabel('# of g1 in c1')
# ax.set_ylabel('hypergeom PMF')
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(x, cdf_g1, 'bo')
# ax.vlines(x, 0, cdf_g1, lw=2)
# ax.set_xlabel('# of g1 in c1')
# ax.set_ylabel('hypergeom PMF')
# plt.show()

# sum(pmf_g[]

# # Cumulative distribution function

# # M = G = total number of objects
# # n = total number of Group 1 objects
# # N = number of samples in random draw. In my case, samples in clade 2. 

# # G = Total number of leaf nodes 
# # n = total number of G1 objects. G1 = G2.
# # N = number of samples in clade 2. ("outgroup")
# # x = number of G1 objects in N drawn without replacement
# #     from the total population
# # x represents the number of 

# x = 2       #(should be 50 % cahnce of drawing each)
# G = 100
# n = 50
# N = x

# #[x, G, n, N] = [100, 50, 100, 30]

# prb = hypergeom.cdf(x, M, n, N)
# #R = hypergeom.rvs(M, n, N, size = 10)


# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(x, prb, 'bo')
# ax.vlines(x, 0, prb, lw=2)
# ax.set_xlabel('# of g1 in c1')
# ax.set_ylabel('hypergeom CDF')
# plt.show()


# # Simulate p values for all simulated values. Put in matrix. 




#%% SAVED VERSION TO PRODUCE SIMULATION PLOT: 
    
    
def simSDR_random(G, G2perm, Random = False, dist = 0.1):
    """
        Simulation of SDR values in a tree comprising 2 clades, C1 and C2,
        separated by a distance (dist).
        The leaf nodes of the tree is defined to belong to either of two 
        groups, G1 and G2.
          
        G = group size of either group. G1 = G2, always.
        C1 = number of samples in clade 1
        C2 = number of samples in clade 2
        
        At start of simulation, G1 and G2 are perfectly separated
        into the two clades; G1 samples are in C1, G2 samples are in C2. 
        Then, G1 = G2 = C1 = C2.
        
        G2perm = number of G2 samples permuted from C2 to C1.
        
        When number of G2 permutations (G2perm) = 0, the groups are perfectly 
        separated into C1 and C2, but as G2perm increase, G2 samples are 
        clustered together with G1 samples in C1, rather than separately in C2. 
        
        Number of possible permutations ranges from 0 to G; in which 
        case all samples are clustered together in the same clade, and C2
        no longer exist.
        
        SDR is calculated as a function of G and perm; group size 
        and number of permutations from G2 to G1. 
        
        Random: 
            Obtain SDR value when partition of G1 and G2 samples into 
            C1 and C2 is random. Sizes of C1 and C2 are determined by the 
            same criteria as above (number of permutations).
            (This will simulate the same situation as if leaf nodes are 
            randomly assigned to G1 and G2.
            
           ((( Result: #X number of SDR values, each obtained from a random 
            tree clustering )))
            
    """
    
    C1 = G + G2perm      # num samples in tree clade 1
    C2 = G - G2perm      # num samples in tree clade 2
    
    if Random: 
        G1C2 = random.randint(0,C2)     # G1C2 = Group 1 samples in clade 2
    else:
        G1C2 = 0
        
    G2C2 = C2 - G1C2                # G2C2 = Group 2 samples in clade 2
    G1C1 = G - G1C2                 # G1C1 = Group 1 samples in clade 1
    G2C1 = G - G2C2                 # G2C1 = Group 2 samples in clade 1
        
    # Total number of pairwise distances between samples of same group
    comb_within_G = ss.comb(G, 2, exact = True) * 2   # Total nr. of combinations between samples of the same defined group
    
    # Distance betweem samples from same group with dist != 0:
    dist_within_G1C1_G1C2 = (G1C1 * G1C2) * dist
    dist_within_G2C1_G2C2 = (G2C1 * G2C2) * dist

    tot_dist_within = dist_within_G1C1_G1C2 + dist_within_G2C1_G2C2
    mean_within = tot_dist_within / comb_within_G


    # Total number of pairwise distances between samples of different groups
    comb_between_G1_G2 = G * G
    
    # Distance between samples form different gorups with dist != 0
    dist_between_G1C1_G2C2 = (G1C1 * G2C2) * dist
    dist_between_G1C2_G2C1 = (G1C2 * G2C1) * dist
    
    tot_dist_between = dist_between_G1C1_G2C2 + dist_between_G1C2_G2C1
    
    mean_between = tot_dist_between / comb_between_G1_G2
    
    # SDR
    try: 
        SDR = mean_within / mean_between
    except: 
        SDR = float("nan")
        print("Division by zero. ")
        print("G1: ", G)
        print("P:", G2perm)
    

    return SDR



#%% Randomize


SDRresults5 = pd.DataFrame(columns = ['Group_size', 'G2Perm', 'SDR'])
#for g1 in range(2,50):
for G in range(2,200):
        p = 0
        while p < G:
            SDR = simSDR_random(G, G2perm = p, Random = False)
            SDRresults5.loc[len(SDRresults5)] = [G, p, SDR]
            p += 1


#%% Visualize 

(ggplot(SDRresults5, aes('G2Perm', 'Group_size', fill = 'SDR'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
# + scale_fill_manual(values=['dodgerblue', 'darkorange'])
# + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR simulation')
)


#%% Visualize p val data

all_simData = load_simData()

(ggplot(all_simData, aes('perm', 'group_size', fill = 'pval'))
 + geom_point(alpha=1, size=2, stroke = 0.1, color = 'indigo')
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
# + scale_fill_manual(values=['dodgerblue', 'darkorange'])
# + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR simulation')
)


(ggplot(all_simData, aes('perm', 'group_size', fill = 'pval_adj'))
 + geom_point(alpha=1, size=3, stroke = 0.1, color = 'indigo')
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
# + scale_fill_manual(values=['dodgerblue', 'darkorange'])
# + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR simulation')
)


