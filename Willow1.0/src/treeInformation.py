# -*- coding: utf-8 -*-
"""
Created on Fri May 21 08:32:27 2021

@author: norab
"""

import pandas as pd
import re
from random import shuffle

"""
Notes:
    
    - Remove super_pops and sub_pops self variables. Replace
    by looping through group_info variable.
    - Before program was planned and developed, samples were given on the form 
    SUP___SUB___SOME-INFO, and the first program made to run SDR calculation
    was hardcoded to take these group levels into account. 
    Since all input data of the thesis is structured this way, the parts
    of the program considering these must remain in this version. 
    However, the plan for future development is to restructure these parts
    of the code for a more flexible group level option.
    - In group info file, give a sample that says the category name instead
    of the group of the instance. 
"""

class treeInfo:
    

    def __init__(self):
        
        """
        Function: 
            initiate class variables
        """
        self.name = None
        self.dist_mat = None
        self.categories = None 
        self.group_info = None
        self.sample_info = None
        self.random_groups = False
        
    def setup(self, dist_mat_file, group_info_file, categories):
        """
        Input: 
            dist_mat_file = matrix containing pairwise sample distances (string)
            group_info_file = contain info on group levels and groups in 
                each group level (string)
            categories = info on group categories
                format: 'category1___category2___category3___...____' (string)        
        """
        
        self.setDistMat(dist_mat_file)
        self.setName(dist_mat_file)
        self.setCategories(categories)
        self.setGroupInfo(group_info_file)
        self.setSampleInfo()
        
    def setDistMat(self, dist_mat_file):
        """
        Input:
            dist_mat_file: filepath to distance matrix
                type: string
                fileformat: CSV
        
        Function: 
            read distance matrix, assign to class variable
        """
        self.dist_mat = pd.read_csv(dist_mat_file, index_col = 0, dtype={'a': str})
    
    def setName(self, dist_mat_file):
        """
        Input: 
            dist_mat_file: filepath to distance matrix
                type: string
                fileformat: CSV
        
        Function: 
            read name of distance matrix, assign to class variable
            
        Note: 
            hard-coded for customized format used in master project.            
            both ENSEBL and gene name identifiers included on the form
            'ENSEBL___GENENAME'
            example: 'ENSG00000000938___FGR'
            
        """
        subName = re.sub('^.*ENS', 'ENS', dist_mat_file)
        self.name = re.sub('___CopD.csv$','', subName)
    
    def setCategories(self, categories):
        """
        Input: 
            group category informaiton on same format as sample names are 
            formatted in distance matrix 
            type: string
            format: 'category1___category2___category3___...____'
            example: 'SUPER___SUB'
            
        Function: 
            assign categories class variable
            
        Develop:
            make separator optional
        """
        self.categories = categories.split('___')
        
    
    def setSampleInfo(self):
        """ 
        Function:
            assign sample information to class variable
            samples are stored in nested dictionary
            keys to nested dictionary are the same integer indexes as 
            samples are referred to in rows and columns of distance matrix
                     
            example: 
                0: {'sample': 'HG02545', 'SUPER': 'AFR', 'SUB': 'ACB'}
                1: {'sample': 'HG00334', 'SUPER': 'EUR', 'SUB': 'FIN'}      
        
        """
        self.sample_info = {}
        
        for ind, sample_in in enumerate(self.dist_mat.columns):
            sample_in = sample_in.split('___')
            # Add sample info
            sample_trans = {'sample':sample_in[-1]}
            
            # Add all categories
            for cat_ind, cat in enumerate(self.categories):
                sample_trans[cat] = sample_in[cat_ind]
    
            # Store sample info
            self.sample_info[ind] = sample_trans

    def shuffleSampleInfo(self):
        """
        Function: 
            shuffle sample information in each category.
            Rows and column indices of distance matrix will no longer
            correspond correctly to their sample informaiton. 

        """

        for cat in self.categories:
            samples = [sample[cat] for sample in self.sample_info.values()]
            shuffle(samples)
            
            for ind in self.sample_info:
                new = self.sample_info[ind]
                new[cat] = samples.pop()
                self.sample_info[ind] = new

        self.random_groups = True

        
    def setUniqseqMap(self, uniqseq_map_file):
        """
        Input: 
            uniqseq_map_file: filepath containing unique sequences. 
            There is one unique sequence per row in the file, and all samples
            having that sequence is listed in the same row. 
            
        Function: 
            set info on number of samples having each of the unique 
            sequences
        
        """
        # Read file
        file = pd.read_csv(uniqseq_map_file, header = None)
        
        # Counts
        self.uniqseq_count = len(file[0])   # What is this for
        self.uniqseq_map = sorted(file.count(axis=1), reverse=True)


    def setGroupInfo(self, group_info_file):
        """
        Input: 
            group_info_file: string filepath to group info-file
        function: 
            organize information into a list, containig a set of defined 
            groups for each group level. 
            E.g: 
                two group levels - "super" and "sub"
                list = [{AFR, EUR, EAS}, {FIN, YRI, GBR}]
        """
        group_info = pd.read_csv(group_info_file, delimiter='\t')
        
        self.group_info = dict()
        for cat in self.categories:
            groups = group_info[group_info['GroupCategory'] == cat]
            groups = set(groups['GroupName'])
            self.group_info[cat] = groups


    def getGroupInfo(self): return self.group_info
    def getSampleInfo(self): return self.sample_info
    def getName(self): return self.name
    def getDistMat(self):return self.dist_mat
    def getUniqseqMap(self): return self.uniqseq_map
    def getUniqseqCount(self): return self.uniqseq_count
    def getCategories(self): return self.categories
    

if __name__ == '__main__':
  
    # Test with simple case
    dist_mat_file = 'C:/Users/norab/Master/Willow1.0/jobs/testOneGene/job_input/geneDists/ENSG00000000938___FGR___CopD.csv'
    group_info_file = 'C:/Users/norab/Master/Willow1.0/jobs/testOneGene/job_input/phydist_population_classes.tsv'
    categories = 'SUPER___SUB'
    test_tree = treeInfo()
    test_tree.setup(dist_mat_file, group_info_file, categories)
    
    sample_info = test_tree.getSampleInfo()
    group_info = test_tree.getGroupInfo()
    categories = test_tree.getCategories()


