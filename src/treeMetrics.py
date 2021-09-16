# -*- coding: utf-8 -*-
"""
Created on Thu May  6 11:33:03 2021

@author: norab

"""

import pandas as pd
import numpy as np
from treeInformation import treeInfo

class treeMetrics(treeInfo):
    
    def __init__(self):
        
        self.pop_dists = None
        self.mean_pop_dists = None
        self.mean_type_dists = None
        self.SDRsuper = None
        self.SDRsub = None
        self.SDVsuper = None
        self.SDVsub = None
        self.singleSuperSDRs = None
        self.singleSubSDRs = None
        
        treeInfo.__init__(self)
    
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

        
    def calcMeanTypeDists(self):
        """
        Calculates mean distance between all pairwise samples within a
        classification type.
        """
        if not self.pop_dists: self.calcPopDists()
        
        dist_summary = {'supBet':0, 'supWith':0, 'subBet':0, 'subWith':0}
        count_summary = {'supBet':0, 'supWith':0, 'subBet':0, 'subWith':0}
        
        for key, val in self.pop_dists.items():
            for pop, dist_count in val.items():
                dist_summary[key] += dist_count[0]
                count_summary[key] += dist_count[1]
        
        self.mean_type_dists = {key: round(dist_summary[key] / count_summary[key], 8)  # Was 8 
                                for key, val in dist_summary.items()}
             
        
    def calcMeanPopDists(self):
        """
        Calculates mean distance between pairwise samples within a population,
        for each classification type. That is, for both super and sub populations,
        the mean inter- and intra distances are calculated and stored in 
        mean_pop_dists. 
        
        """
        
        if not self.pop_dists: self.calcPopDists()
        
        mean_pop_dists = {}
        
        for key, val in self.pop_dists.items():
            type_means = {}
            
            for pop, dist_count in val.items():
                if dist_count[1]:  
                    type_means[pop] = round(dist_count[0] / dist_count[1], 6) # was 6
                else:   # No distance calculated for this pop
                    type_means[pop] = 0
             
            mean_pop_dists[key] = type_means
             
        self.mean_pop_dists = mean_pop_dists
        
        
    def calcSDR(self):
        """
        Input: 
            group_dists: Numpy arrayList of dataframes containing info of total group distances and
            number of comparisons.
            calc_SDRgroupwise: If SDR for all single groups should be calculated. This is 
                required to calculate SDVs. 
        Function: 
            Calculates group mean for either full group (calc_SDRgroupwise = False)
                or every single group (.. = True.)
        Testing: OK.   
        """

        if not self.mean_type_dists: 
            self.calcMeanTypeDists()
        elif (self.mean_type_dists and self.random_pops):
            self.calcMeanTypeDists()
        
            
        if self.mean_type_dists['supBet']:
            self.SDRsuper = round(self.mean_type_dists['supWith'] / 
                                  self.mean_type_dists['supBet'], 6) # was 6
        else: 
            self.SDRsuper = float('NaN')
        
        if self.mean_type_dists['subBet']:
            self.SDRsub = round(self.mean_type_dists['subWith'] / 
                                self.mean_type_dists['subBet'], 6) # was 6
        else:
            self.SDRsub = float('NaN')
        
    
    def calcSingleSDRs(self):
        """
        Ratio of mean intra population distance to mean inter population distance
        for each population. 
        
        Example: Africa
            mean(distances between all African samples) /
            mean(distances between all African and non-African samples )
            
        Used to calculate SDV and as indicational measure.
        
        Test: OK.
        """
        
        if not self.mean_pop_dists : self.calcMeanPopDists()
        
        self.singleSuperSDRs = {}
        self.singleSubSDRs = {}
        #frame = pd.DataFrame(self.mean_pop_dists)
        
       
        for pop in self.super_pops:
            ratio = float('NaN')
            if self.mean_pop_dists['supBet'][pop]:
                ratio = (self.mean_pop_dists['supWith'][pop] / 
                        self.mean_pop_dists['supBet'][pop])

            self.singleSuperSDRs[pop] = round(ratio, 3)
        
        for pop in self.sub_pops: 
            ratio = float('NaN')
            if self.mean_pop_dists['subBet'][pop]:
                ratio = (self.mean_pop_dists['subWith'][pop] / 
                         self.mean_pop_dists['subBet'][pop])

            self.singleSubSDRs[pop] = round(ratio, 3)

    
    def calcSDV(self):
        """
        Measure of variability of single SDRs.
        Gives measure of difference in cluster structures. 
        Testing: to do
        """
        if not self.singleSuperSDRs: self.calcSingleSDRs()
        
        # Test if this works
        self.SDVsuper = round(np.array(list(self.singleSuperSDRs.values())).var(ddof=1), 4)
        self.SDVsub = round(np.array(list(self.singleSubSDRs.values())).var(ddof=1), 4)
    
    def calcNonZeroTotal(self, 
                         percent = True):
        """

        Testing: to do
        Input:
            dist_mat : DataFrame with cophenetic distances for a tree
    
        Returns:
            totDist: float. total amount of pairdistances being non-zero as 
            percent or count. 
            
        Note: This does the same as what u did in R. Have all values alread. 
        """
        nonZero_row = pd.DataFrame((self.dist_mat != 0).astype(int).sum(axis=1))
        nonZeroCount = int(nonZero_row.sum(axis=0))
        
        num_entries = (self.dist_mat.shape[0] * self.dist_mat.shape[1]) - self.dist_mat.shape[0]
        nonZeroPercent = nonZeroCount / num_entries
        if percent:    
            return nonZeroPercent
        else: 
            return nonZeroCount
        
    
    def calcNonZerosForSamples(self, percent = True):
        """
        Input:
            dist_mat : DataFrame with cophenetic distances for a tree
    
        Returns:
            totDist: dataframe with total amount of pairdistances being non-zero as 
            percent and count per row. 
        """
         
        nonZero_row = pd.DataFrame((self.dist_mat != 0).astype(int).sum(axis=1), columns = ['nonZero_count'])
        nonZero_row['percent'] = round(nonZero_row['nonZero_count'] / (nonZero_row.shape[0] - 1), 4)
        
        return nonZero_row
    
    # test = calcNonZerosForSamples(MPP5_cd, percent=False)
    # test2 = calcNonZeroTotal(MPP5_cd,)
    
    def calcNonZerosForPops(self, popType='all', percent = True):
        """
        Input: 
            dist_mat : DataFrame with cophenetic distances for a tree
            sample_info: a dict with numeric index values as keys and sample
            infor as value, stored on format [gene_name, superpop, subpop]
            pop_info: Defines populations in popType. 
        Function:
            COunt number of 
        
        Returns:
            Amount of non-zero distances for each subpop
            
        TO DO: 
            create pop_info and pops with functions within this function.
            Have to fix imports and project structure first. 
        """        
        # find total sample cells in matrix
        sampleCount = pd.DataFrame(self.countPopSamples(self.sample_info, popType=popType))
        sampleCount['totalPopComp'] = (sampleCount['sampleCount'] * self.dist_mat.shape[0]) - sampleCount['sampleCount']    
        
        # Make dict counter for pops --- TEST!!!
        if popType=='super':
            sup_sums = {pop : 0 for pop in self.pop_info[0]}   # Placeholder for sum values
            sup_sample_info = dict([val[0], val[1]] for val in self.sample_info.values())
            
        elif popType=='sub':
            sub_sums = {pop : 0 for pop in self.pop_info[1]}   # Placeholder for sum values
            sub_sample_info = dict([val[0], val[2]] for val in self.sample_info.values())
            
        elif popType=='all':
            sup_sums = {pop : 0 for pop in self.pop_info[0]}   # Placeholder for sum values
            sup_sample_info = dict([val[0], val[1]] for val in self.sample_info.values())
            
            sub_sums = {pop : 0 for pop in self.pop_info[1]}   # Placeholder for sum values
            sub_sample_info = dict([val[0], val[2]] for val in self.sample_info.values())
            
        #     sums1 = {pop : 0 for pop in pop_info[0]}   # Placeholder for sum values
        #     sums2 = {pop : 0 for pop in pop_info[1]}
        #     sums = dict(sums1, **sums2)
            
        #     sample_info1 = dict([val[0], val[1]] for val in sample_info.values())
        #     sample_info2 = dict([val[0], val[2]] for val in sample_info.values())
            
        #     sample_info = {**sample_info1, **sample_info2}
            
        #     return sample_info1
        else: 
            raise ValueError("Invalid calssification type. ")
            
        # map nonZero and sample_info together into matrix
        
        nonZeros = self.calcNonZerosForSamples(self.dist_mat, percent = True)  
        if popType == 'all':
            nonZeros['supPop'] = nonZeros.index.map(sup_sample_info)
            nonZeros['subPop'] = nonZeros.index.map(sub_sample_info)
            for key, val in nonZeros.iterrows():
                sup_sums[val[2]] += val[0]
                sub_sums[val[3]] += val[0]
        
            sums = {**sup_sums,**sub_sums}
            
        elif popType == 'super':
            raise ValueError("This option is not yet available")
        else: 
            raise ValueError("This option is not yet available")
            
        # Divide non-zero values for pop on total pop compariosns    
        sampleCount['popSums'] = sampleCount.index.map(sums)
        sampleCount['percentNonZeroForPop'] = sampleCount.popSums / sampleCount.totalPopComp
        sampleCount['gene'] = self.getGeneName()
        
        return sampleCount
        
    def getPopDists(self):
        if not self.pop_dists: self.calcPopDists()
        return self.pop_dists
    
    def getMeanTypeDists(self):
        if not self.mean_type_dists: self.calcMeanTypeDists()
        return self.mean_type_dists
        
    def getMeanPopDists(self): 
        if not self.mean_pop_dists: self.calcMeanPopDists()
        return self.mean_pop_dists
    
    def getSuperPops(self): return self.super_pops
    def getSubPops(self): return self.sub_pops
    def getSingleSuperSDRs(self): return self.singleSuperSDRs
    def getSingleSubSDRs(self): return self.singleSubSDRs
    def getSDRsuper(self): return self.SDRsuper
    def getSDRsub(self): return self.SDRsub
    def getSDVsuper(self): return self.SDVsuper
    def getSDVsub(self): return self.SDVsub

