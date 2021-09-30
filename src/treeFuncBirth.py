# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 11:48:34 2021

@author: norab
"""

"""

    Idea: Make calcPopDists work for one subtype only. 
    Time difference of looping once and twice. 
    
"""

import pandas as pd
# NEW

def calcOnePopDists(self):
        """
        Calculate inter and intra population distances
        for both classification types (Super and sub).
        These are merged into one function for efficiency reasons.
    
        """
        dist_mat = self.dist_mat.to_numpy()
        
        # Create distance sum placeholders
        # Dict with pop as key and [dist value, number of dist values added] as value
        WithSums = {pop : [0,0] for pop in self.pop_info[0]}
        BetSums = {pop : [0,0] for pop in self.pop_info[0]}

        # Iter upper triangular dist_mat
        row_length = len(dist_mat)
        row_start = 1
        for col in range(row_length):
            for row in range(row_start, row_length):
                sample1 = col
                sample2 = row
                
                # Get popName
                Pop1 = self.sample_info[sample1][1]
                Pop2 = self.sample_info[sample2][1]
                
                val = dist_mat[sample2][sample1]        # Distance value

                if Pop1 == Pop2:        # Same super pop, same sub pop
                    WithSums[Pop1] = [i+j for i, j in zip(WithSums[Pop1], [val, 1])]  
                else:
                    BetSums[Pop1] = [i+j for i, j in zip(BetSums[Pop1], [val, 1])]  
                    BetSums[Pop2] = [i+j for i, j in zip(BetSums[Pop2], [val, 1])]  
                        
            row_start += 1
            
        self.pop_dists = {'With': WithSums, 'Bet': BetSums}



# OLD
def calcPopDists(self):
        """
        Calculate inter and intra population distances
        for both classification types (Super and sub).
        These are merged into one function for efficiency reasons.
    
        """
        dist_mat = self.dist_mat.to_numpy()
        
        # Create distance sum placeholders
        # Dict with pop as key and [dist value, number of dist values added] as value
        supWithSums = {pop : [0,0] for pop in self.pop_info[0]}
        supBetSums = {pop : [0,0] for pop in self.pop_info[0]}
        subWithSums = {pop : [0,0] for pop in self.pop_info[1]}   
        subBetSums = {pop : [0,0] for pop in self.pop_info[1]}

        # Iter upper triangular dist_mat
        row_length = len(dist_mat)
        row_start = 1
        for col in range(row_length):
            for row in range(row_start, row_length):
                sample1 = col
                sample2 = row
                
                # Get popName
                supPop1 = self.sample_info[sample1][1]
                supPop2 = self.sample_info[sample2][1]
                subPop1 = self.sample_info[sample1][2]
                subPop2 = self.sample_info[sample2][2]
                
                val = dist_mat[sample2][sample1]        # Distance value
    # =============================================================================
    #             Within
    # =============================================================================
                if subPop1 == subPop2:        # Same super pop, same sub pop
                    if supPop1 == supPop2:
                        # Super
                        supWithSums[supPop1] = [i+j for i, j in zip(supWithSums[supPop1], [val, 1])]  
                        # Sub
                        subWithSums[subPop1] = [i+j for i, j in zip(subWithSums[subPop1], [val, 1])] 
                    else:
                        # Super
                        supBetSums[supPop1] = [i+j for i, j in zip(supBetSums[supPop1], [val, 1])]  
                        supBetSums[supPop2] = [i+j for i, j in zip(supBetSums[supPop2], [val, 1])]  
                        # Sub
                        subWithSums[subPop1] = [i+j for i, j in zip(subWithSums[subPop1], [val, 1])] 
                        
    # =============================================================================
    #             Within and Between
    # =============================================================================
                elif supPop1 == supPop2:      # Same super pop, different sub-pop
                    # Super within
                    supWithSums[supPop1] = [i+j for i, j in zip(supWithSums[supPop1], [val, 1])]
                    
                    # Sub between
                    subBetSums[subPop1] = [i+j for i, j in zip(subBetSums[subPop1], [val, 1])]
                    subBetSums[subPop2] = [i+j for i, j in zip(subBetSums[subPop2], [val, 1])] 
    # =============================================================================
    #             Between
    # =============================================================================
                else:                           # All different: add all to between - pops            
                    # Super between
                    supBetSums[supPop1] = [i+j for i, j in zip(supBetSums[supPop1], [val, 1])]
                    supBetSums[supPop2] = [i+j for i, j in zip(supBetSums[supPop2], [val, 1])]
                    
                    # Sub between
                    subBetSums[subPop1] = [i+j for i, j in zip(subBetSums[subPop1], [val, 1])]
                    subBetSums[subPop2] = [i+j for i, j in zip(subBetSums[subPop2], [val, 1])]
            row_start += 1
            
        self.pop_dists = {'supWith': supWithSums, 'supBet': supBetSums, 
                'subWith': subWithSums, 'subBet': subBetSums}
        
        
        
#%% SetPopInfo: 
   
# New

def addSubtypeInfo(subtype_info_file):
    """
    file: string filepath to subtype type info-file.
          pop types = super and sub
    function: organize information into list of two sets 
            of contained populations in super- and sub pops respectivly
    """
    subtype_info = pd.read_csv(subtype_info_file, delimiter='\t')
    subtype_levels = len(subtype_info['ClassificationType'].unique())

    if subtype_levels == 2:     # Old version - still runable.
        subtype_levels = 2 # self
        super_pops = subtype_info.loc[subtype_info['ClassificationType'] == 'SUPER']
        super_pops = set(super_pops['ClassificationName'])
        sub_pops = subtype_info.loc[subtype_info['ClassificationType'] == 'SUB']
        sub_pops = set(sub_pops['ClassificationName'])
        self.pop_info = [super_pops, sub_pops]   # to be replaced/deleted
        self.super_pops = super_pops
        self.sub_pops = sub_pops
        
    elif subtype_levels == 1:    # New version
        subtype_levels = 1 # self
        subtypes = subtype_info.loc[subtype_info['ClassificationType'] == 'SUPER']
        subtype_names = set(subtypes['ClassificationName'])
        
        # EDIT IN NEW VERSION:
        self.subtype_info = [subtype_names] # self
        self.subtype_names = subtype_names
        
    else: 
        raise ValueError("subtype levels error. Number of levels found: ", subtype_levels)
        
        
    #%% 
# Test

subtype_info_file = 'C:/Users/norab/Master/data/job_simulation/job_input/simulation_group_classes.tsv'
test_popinfo = addPopInfo(test_file)

#%%
# Old
def setPopInfo(self, pop_info_file):
    """
    file: string filepath to pop type info-file.
          pop types = super and sub
    function: organize information into list of two sets 
            of contained populations in super- and sub pops respectivly
    """
    pop_info = pd.read_csv(pop_info_file, delimiter='\t')

    super_pops = pop_info.loc[pop_info['ClassificationType'] == 'SUPER']
    super_pops = set(super_pops['ClassificationName'])
    sub_pops = pop_info.loc[pop_info['ClassificationType'] == 'SUB']
    sub_pops = set(sub_pops['ClassificationName'])
    
    self.pop_info = [super_pops, sub_pops]   # to be replaced/deleted
    self.super_pops = super_pops
    self.sub_pops = sub_pops
    

    
        
        
        
        
        
        
        
