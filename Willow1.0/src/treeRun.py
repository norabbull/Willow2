# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 05:48:23 2021

@author: norab
"""

import csv
from datetime import datetime
import pandas as pd
import yaml
import sys
import os
from os.path import isfile, join
from src.treeMetrics import treeMetrics
#import logging
#import treeLogger

#logger = treeLogger.get_logger(__name__)

class RunStuff:
    
    def __init__(self):
        
        """

        """
        
        self.config()
        self.file_list = self.make_filelist(self.input_dist_folder)
        
        # Run! 
        self.run_func()
        #self.logger = treeLogger.get_logger(__name__)
        #self.logger.info('RunStuff initiated')
        #self.logger.info("It made it here! Hurray!")

    
    def config(self, config):
        
        """
        Input: 
            yaml format configuration file 
            
        ---
     func: 'calcTest'
     input_folder: 'C:/Users/norab/Master/Willow1.0/jobs/testOneGene/job_input/'
     input_group_info: 'phydist_population_classes.tsv'
     input_dist_folder: 'C:/Users/norab/Master/Willow1.0/jobs/testOneGene/job_input/geneDists/'
     output_folder: 'C:/Users/norab/Master/Willow1.0/jobs/testOneGene/job_output/'
     output_log: 'testOneGene_log_datetime.log'
     output_values: 'testOneGene_datetime.csv'
     output_unprocessed: 'testOneGene_unprocessed_files_datetime.csv'
     num_rand_trees: 1
     skip_genes: 'None'
     select_genes: 'None'
     categories: 'SUPER___SUB'
        
        """
        self.config = config
        
        # Make datetime stamp
        for key, val in self.config.items():
            if 'datetime' in str(val): 
                self.config[key] = val.replace('datetime', datetime.now().strftime("%d.%m.%Y_%H.%M"))
        
        # Folder/file paths
        self.input_folder = config.get('input_folder').strip()
        self.input_dist_folder = config.get('input_cd_folder').strip()
        self.output_folder = config.get('output_folder').strip()
        self.output_values_file = self.output_folder + config.get('output_values').strip()
        self.output_unprocessed = self.output_folder + config.get('output_unprocessed').strip()
               
        # Input information
        self.func = config.get('func').strip()
        self.group_info = self.input_folder + config.get('input_group_info').strip()
        self.filter_select = self.input_folder + config.get('filter_select').strip()
        self.filter_skip = self.output_folder + config.get('filter_skip').strip()
        self.categories = config.get('categories').strip()
        self.num_random_values = config.get('num_random_values').strip()
            
    @classmethod
    def make_filelist(self, input_files):
        try: 
            if isinstance(input_files, str):
                if '.csv' in input_files: 
                    files = pd.read_csv(input_files, header = None)
                    files = files[0].values.tolist()
                else: 
                    files = [join(input_files, f) for f in os.listdir(input_files) 
                                 if isfile(join(input_files, f))]
                    files = [f.strip() for f in files]
        except Exception as e:
            #logger.error("Failed to make filelist.")
            print("Error with make_filelist")
            print(e)
                
        return files
                    
    def run_calcSDR(self):
        """
        Funcition: 
            calculate SDR value for all given input distance matrix files
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
            #logger.info(self.file_list)
            
            ind = 1
            ind_len = len(self.file_list)
            
            if self.random: 
                num_iter = int(self.num_random_values)
            else:
                num_iter = 1

            for i in range(num_iter):
                for dist_file in self.file_list:
                    try:         
                        
                        #logger.info("File processed: {0}".format(dist_file))
                        #logger.info("File number: {0} / {1}".format(ind, ind_len))
                        #logger.info("Iter round: {0}".format(i))
                        ind +=1
    
                        tree = treeMetrics()
                        tree.setup(dist_file.strip(), self.group_info, categories)
                        
                        if random: 
                            tree.shuffleSampleInfo()
    
                        tree.calcSDR()
                        SDR_save = [tree.getGeneName(), SDR]
                    
                        with open(output_SDR, 'a', newline='') as f:   # write to file    
                            writer = csv.writer(f)
                            writer.writerow(SDR_save)

        
                    except Exception: 
                       
                        #logger.exception("File disrupted:", dist_file)
                        file = str(dist_file)
                        
                        with open(output_unprocessed, 'a') as f: 
                            writer = csv.writer(f)
                            writer.writerow(file)

                        pass
                
        except Exception as e: 
            #logger.exception(e)
            print("error with runCalc")
            print(e)


       
    def run_calcTest(self):
        """
        Test function to check functionality of program.
        """
        return 


    def run_func(self):
        if self.func == "calcSDR":
            self.run_calcSDR()
        elif self.func == "calcSDRnull":
            self.random = True
            self.run_calcSDR()
        elif self.func == "calcTest":
            return self.run_calcTest()
        else:
            print("No function called ", self.func)
            #logger.error("Not a valid function option. Change in main_config file." +
                              #"options are: calcSDR, calcSDV, calcSingleSDRs, calcNullSDR and calcTest")
   
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
    cat_dists = test_tree.getCategoryDists()
    
    test_tree.calcMeanCategoryDists()
    test_tree.calcSDR()    
    
    