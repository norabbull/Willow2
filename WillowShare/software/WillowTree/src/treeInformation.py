# -*- coding: utf-8 -*-
"""
Created on Fri May 21 08:32:27 2021

@author: norab
"""

import pandas as pd
import re
from random import shuffle

class treeInfo:
    

    def __init__(self):
        
        """
        dist_mat = matrix containing cophenetic distances
        sample_info = indices of sample names 
        pop_info = list of sets with population names.Index 0 = super, index 1 = sub 

        """
        self.dist_mat = None
        self.pop_info = None
        self.super_pops = None
        self.sub_pops = None
        self.sample_info = None
        self.gene_name = None
        self.pop_dists = None
        self.mean_pop_dists = None
        self.mean_type_dists = None    
        self.random_pops = False
        self.uniqseq_map = None
        self.uniqseq_count= None
                
    def setup(self, dist_mat_file, pop_info_file):
        self.setDistMat(dist_mat_file)
        self.setPopInfo(pop_info_file)
        self.setGeneName(dist_mat_file)
        self.setSampleInfo()
    
    def setDistMat(self, dist_mat_file):
        self.dist_mat = pd.read_csv(dist_mat_file, index_col = 0, dtype={'a': str})
    
    def setGeneName(self, dist_mat_file):
        """
        Input: 
            file: string filepath to distance matrix-file. 
        Function: 
            Filters out name of gene the tree represents and assign to
            class variable "gene_name".
            Both Ensembl and gene name identifiers included on the form: 
                'ENSG00000000938___FGR'
        """
        subName = re.sub('^.*ENS', 'ENS', dist_mat_file)
        self.gene_name = re.sub('___CopD.csv$','', subName)
    
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
     
    def setSampleInfo(self):
        """
        Input:
        dist_mat: dataframe of cophenetic distances between all pairs in tree
        
        Function:
        structure sample_info: dict with created int as index and sample info stored
        in list as value: [gene name, super pop name, sub pop name]
        
        Reason for numerical indices: Can use to iterate triangular matrix in 
        calcTypeDists
        """
        self.sample_info = {}
        
        for ind, sample in enumerate(self.dist_mat.columns):
            subName = re.sub('^...___', '', sample)
            subName = re.sub('___.*$','', subName)
            superName = re.sub('_.*$','', sample)
            
            self.sample_info[ind] = [sample, superName, subName]
    
    def shuffleSampleInfo(self):
        """
        Make randomly defined population groups. 
        To be used for creating null-distribution.
            
        """
          
        sup = [val[1] for val in self.sample_info.values()]
        sub = [val[2] for val in self.sample_info.values()]
        shuffle(sup)
        shuffle(sub)

        for key, val in self.sample_info.items():
            val[1] = sup.pop()
            val[2] = sub.pop()
            self.sample_info[key] = val
        
        self.random_pops = True
        self.mean_pop_dists = None
        self.mean_type_dists = None
        self.pop_dists = None

    def setUniqseqMap(self, uniqseq_map_file):
        """
        Set info on number of sequences there are of each unique sequence
        that makes up the tree. 
        """
        # Read file
        file = pd.read_csv(uniqseq_map_file, header = None)
        
        # Counts
        self.uniqseq_count = len(file[0])
        self.uniqseq_map = sorted(file.count(axis=1), reverse=True)
        
    def getSampleInfo(self): return self.sample_info
    def getPopInfo(self): return self.pop_info
    def getGeneName(self): return self.gene_name
    def getDistMat(self):return self.dist_mat
    def getUniqseqMap(self): return self.uniqseq_map
    def getUniqseqCount(self): return self.uniqseq_count


if __name__ == '__main__':
  pass



