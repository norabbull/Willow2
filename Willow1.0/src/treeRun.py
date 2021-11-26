# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 05:48:23 2021

@author: norab
"""

import csv
from datetime import datetime
import pandas as pd
# import yaml
# import sys
from src.treeMetrics import treeMetrics


class treeRun:
    """
        class to organize tree instance calculations for files
    """
    
    def __init__(self):
        
        """
        Function: 
            initiate class variables
        """
        
        self.config()
        self.file_list = self.make_filelist(self.input_dist_folder)
        
        # Run! 
        self.run_func()

    
    def config(self, config):
        
        """
        Input: 
            config: configuration filepath  (yaml)
            type: str
                       
        Function: 
            configure variables from file
        
        """
        self.config = config
        
        # Make datetime stamp
        for key, val in self.config.items():
            if 'datetime' in str(val): 
                self.config[key] = val.replace('datetime', datetime.now().strftime("%d.%m.%Y_%H.%M"))
        
        # Folder/file paths
        #self.input_folder = config.get('input_folder').strip()
        self.input_dist_folder = config.get('input_cd_folder').strip()
        self.output_folder = config.get('output_folder').strip()
        self.output_values_file = self.output_folder + config.get('output_values').strip()
        self.output_unprocessed = self.output_folder + config.get('output_unprocessed').strip()
               
        # Input information
        self.func = config.get('func').strip()
        self.group_info = config.get('input_group_info').strip()
        self.filter_select = config.get('filter_select').strip()
        self.filter_skip = config.get('filter_skip').strip()
        self.categories = config.get('categories').strip()
        self.num_random_values = config.get('num_random_values').strip()
            
                    
    def run_calcGDR(self):
        """
        Funcition: 
            - calculate GDR value for all given input distance matrix files
            - save values to file
            - if error with file, filepath is saved to separate file and
            calculation proceeds with next file
        """
        
        try:        
            if '.csv' in self.filter_select:
                filter_select = pd.read_csv(self.filter_select, header=None)
                filter_select.columns = ['name']
                filter_select= list(filter_select['name'])
                
                keep_files = []
                for file in self.file_list:
                    if any(gene in file for gene in filter_select):
                        keep_files.append(file)
                        self.file_list = keep_files
            
            if '.csv' in self.filter_skip:
                filter_skip= pd.read_csv(self.filter_skip, header=None)
                filter_skip.columns = ['name']
                filter_skip= list(filter_skip['name'])
                
                remove_files = []
                for file in self.file_list:
                    if any(gene in file for gene in filter_skip):
                        remove_files.append(file)
        
                self.file_list = [f for f in self.file_list if f not in remove_files]
            
            ind = 1
            #ind_len = len(self.file_list)
            
            if self.random: 
                num_iter = int(self.num_random_values)
            else:
                num_iter = 1

            for i in range(num_iter):
                for dist_file in self.file_list:
                    try:         
                        ind +=1
    
                        tree = treeMetrics()
                        tree.setup(dist_file.strip(), self.group_info, categories)
                        
                        if self.random: 
                            tree.shuffleSampleInfo()
    
                        tree.calcGDR()
                        for cat in tree.getCategories():
                            
                            GDR = tree.getGDR()[]
                            GDR_save = [tree.getGeneName(), tree.getGDR()]
                        
                            with open(self.output_GDR, 'a', newline='') as f:   # write to file    
                                writer = csv.writer(f)
                                writer.writerow(GDR_save)

        
                    except Exception: 
                        
                        # Ssve filepath to unprocessed file
                        file = str(dist_file)
                        with open(self.output_unprocessed, 'a') as f: 
                            writer = csv.writer(f)
                            writer.writerow(file)
                        
                        # continue with next
                        pass
                
        except Exception as e: 
            print("error with runCalc:\n")
            print(e)


       
    def run_calcTest(self):
        """
        Test function to check functionality of program.
        """
        return 


    def run_func(self):
        if self.func == "calcGDR":
            self.run_calcGDR()
        elif self.func == "calcGDRnull":
            self.random = True
            self.run_calcGDR()
        elif self.func == "calcTest":
            return self.run_calcTest()
        else:
            print("No valid function called ", self.func)
   
        #logger.info('Session finished successfully: {0}'.format(datetime.now().strftime("%d.%m.%Y_%H.%M")))

#%%

import timeit

if __name__ == '__main__':
    # Test with simple case
   
    dist_mat_file = 'C:/Users/norab/Master/data/real_tree_data/dist_mat_subset/ENSG00000001626___CFTR___CopD.csv'
    #dist_mat_file = 'C:/Users/norab/Master/Willow1.0/jobs/testOneGene/job_input/geneDists/ENSG00000000938___FGR___CopD.csv'
    group_info_file = 'C:/Users/norab/Master/Willow1.0/jobs/testOneGene/job_input/phydist_population_classes.tsv'
    
    test_tree = treeMetrics()
    test_tree.setup(dist_mat_file, group_info_file, categories)
    
    sample_info = test_tree.getSampleInfo()
    group_info = test_tree.getGroupInfo()
    categories = test_tree.getCategories()
    
    # Calcl group dists
    starttime = timeit.default_timer()
    timeGroupDists = test_tree.calcGroupDists()
    endtime = timeit.default_timer()
    cat_dists = test_tree.getGroupDists()
    
    test_tree.calcMeanGroupDists()
    test_tree.calcGDR()    
    
    