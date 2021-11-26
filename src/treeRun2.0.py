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
    
    def __init__(self, config):
        
        """
            confige_file: dicitonary of 

        """
        
        #self.logger = treeLogger.get_logger(__name__)
        #self.logger.info('RunStuff initiated')
        self.config = config
        
        # Make datetime stamp
        for key, val in self.config.items():
            if 'datetime' in str(val): 
                self.config[key] = val.replace('datetime', datetime.now().strftime("%d.%m.%Y_%H.%M"))
        
        self.input_folder = self.config.get('input_folder').strip()
        self.output_folder = self.config.get('output_folder').strip()
        self.skip_genes = self.output_folder + self.config.get('skip_genes').strip()
        self.func = self.config.get('func').strip()
        #self.logger.info("It made it here! Hurray!")
        
        # Run! 
        self.run_func()
        

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
        
    def run_calcSDR_old(self, random=False):
        """
        random: If groups should be random or not. Used to calculate null distribution.
        skip_genes: list of gene names to be skiped from calculation. 

        """
        
        try:
            
            input_cd_folder = self.config.get('input_cd_folder').strip()
            group_info = self.input_folder + self.config.get('input_group_info').strip()
            
            output_super = self.output_folder + self.config.get('output_super').strip()
            output_sub = self.output_folder + self.config.get('output_sub').strip()
            output_unprocessed = self.output_folder + self.config.get('output_unprocessed').strip()
            
            file_list = self.make_filelist(input_cd_folder)
            skip_genes = self.input_folder + self.config.get('skip_genes').strip()
            select_genes = self.input_folder + self.config.get('select_genes').strip()
            
            
            if '.csv' in select_genes:
                select_genes= pd.read_csv(select_genes, header=None)
                select_genes.columns = ['gene']
                select_genes= list(select_genes['gene'])
                
                keep_files = []
                for file in file_list:
                    if any(gene in file for gene in select_genes):
                        keep_files.append(file)
                        file_list = keep_files
            
            if '.csv' in skip_genes:
                skip_genes= pd.read_csv(skip_genes, header=None)
                skip_genes.columns = ['gene']
                skip_genes= list(skip_genes['gene'])
                
                remove_files = []
                for file in file_list:
                    if any(gene in file for gene in skip_genes):
                        remove_files.append(file)
        
                file_list = [f for f in file_list if f not in remove_files]
            #logger.info(file_list)
            
            ind = 1
            ind_len = len(file_list)
            
            if random: 
                num_iter = int(self.config.get('num_rand_trees'))
            else:
                num_iter = 1

            for i in range(num_iter):
                for cd_file in file_list:
                    try:         
                        
                        #logger.info("File processed: {0}".format(cd_file))
                        #logger.info("File number: {0} / {1}".format(ind, ind_len))
                        #logger.info("Iter round: {0}".format(i))
                        ind +=1
    
                        tree = treeMetrics()
                        tree.setup(cd_file.strip(), group_info)
                        
                        if random: 
                            tree.shuffleSampleInfo()
    
                        tree.calcSDR()
                        
                        supSDR = tree.getSDRsuper()
                        subSDR = tree.getSDRsub()
                        
                        supSDR = [tree.getGeneName(), supSDR]
                        subSDR = [tree.getGeneName(), subSDR]
                    
                        with open(output_super, 'a', newline='') as f:   # write to file    
                            writer = csv.writer(f)
                            writer.writerow(supSDR)
                            
                        with open(output_sub, 'a', newline='') as f:   # write to file    
                            writer = csv.writer(f)
                            writer.writerow(subSDR)
        
                    except Exception: 
                       
                        #logger.exception("File disrupted:", cd_file)
                        file = str(cd_file)
                        
                        with open(output_unprocessed, 'a') as f: 
                            f.write(file)
                            #open text file
                        pass
                
        except Exception as e: 
            #logger.exception(e)
            print("error with runCalc")
            print(e)
            
            
    # NEW IN 3.0:
    def run_calcSDR(self, random=False):
        """
        random: If groups should be random or not. Used to calculate null distribution.
        skip_genes: list of gene names to be skiped from calculation. 

        """
        
        try:
            
            input_cd_folder = self.config.get('input_cd_folder').strip()
            subtype_info = self.input_folder + self.config.get('input_subtype_info').strip()
            output_SDR = self.output_folder + self.config.get('output_SDR').strip()
            output_unprocessed = self.output_folder + self.config.get('output_unprocessed').strip()
            
            file_list = self.make_filelist(input_cd_folder)
            skip_genes = self.input_folder + self.config.get('skip_genes').strip()
            select_genes = self.input_folder + self.config.get('select_genes').strip()
            
            
            if '.csv' in select_genes:
                select_genes= pd.read_csv(select_genes, header=None)
                select_genes.columns = ['gene']
                select_genes= list(select_genes['gene'])
                
                keep_files = []
                for file in file_list:
                    if any(gene in file for gene in select_genes):
                        keep_files.append(file)
                        file_list = keep_files
            
            if '.csv' in skip_genes:
                skip_genes= pd.read_csv(skip_genes, header=None)
                skip_genes.columns = ['gene']
                skip_genes= list(skip_genes['gene'])
                
                remove_files = []
                for file in file_list:
                    if any(gene in file for gene in skip_genes):
                        remove_files.append(file)
        
                file_list = [f for f in file_list if f not in remove_files]
            #logger.info(file_list)
            
            ind = 1
            ind_len = len(file_list)
            
            if random: 
                num_iter = int(self.config.get('num_rand_trees'))
            else:
                num_iter = 1

            for i in range(num_iter):
                for cd_file in file_list:
                    try:         
                        
                        #logger.info("File processed: {0}".format(cd_file))
                        #logger.info("File number: {0} / {1}".format(ind, ind_len))
                        #logger.info("Iter round: {0}".format(i))
                        ind +=1
    
                        tree = treeMetrics()
                        tree.setup(cd_file.strip(), subtype_info)
                        
                        if random: 
                            tree.shuffleSampleInfo()
    
                        tree.calcSDR()
                        
                        # Old: 
                        # supSDR = tree.getSDRsuper()
                        # subSDR = tree.getSDRsub()
                        
                        # supSDR = [tree.getGeneName(), supSDR]
                        # subSDR = [tree.getGeneName(), subSDR]
                        
                        # New in 3.0
                        SDR = tree.getSDR()
                        SDR_save = [tree.getGeneName(), SDR]
                    
                        with open(output_SDR, 'a', newline='') as f:   # write to file    
                            writer = csv.writer(f)
                            writer.writerow(SDR_save)

        
                    except Exception: 
                       
                        #logger.exception("File disrupted:", cd_file)
                        file = str(cd_file)
                        
                        with open(output_unprocessed, 'a') as f: 
                            f.write(file)
                            #open text file
                        pass
                
        except Exception as e: 
            #logger.exception(e)
            print("error with runCalc")
            print(e)
       
    def run_calcSingleSDRs(self, random = False):
        
        """
        Save one file for each gene. 
        File example: 
            
            Level   population    SingleSDR
            Super   AFR           0.4
            Super   EUR           0.5
            Sub     FIN           0.8
            ...     ...           ...
            
        """

        input_files = self.config.get('input_files')
        group_info = self.config.get('input_group_info')
        SSDR_output_dir = self.config.get('output_SSDR')
        output_unprocessed = self.config.get('output_unprocessed')
               
        file_list = self.make_filelist(input_files)

        ind = 1
        ind_len = len(file_list)
            
        logger.info("Files to process:\n", file_list)

        try:
            while file_list:
                dist_file = file_list.pop().strip()
                
                logger.info("File processed: {0}".format(dist_file))   # TO DO: convert to log  
                logger.info("Number: {0} / {1}".format(ind, ind_len))
                
                tree = treeMetrics()
                tree.setup(dist_file, group_info)
                if random: 
                    tree.shuffleSampleInfo()
                    
                tree.calcSingleSDRs()
                
                supSDRs = tree.getSingleSuperSDRs()
                subSDRs = tree.getSingleSubSDRs()

                df_sup = pd.DataFrame([['super', group, val] for group, val in supSDRs.items()],
                                   columns=['level', 'group', 'singleSDR'])
                df_sub = pd.DataFrame([['sub', group, val] for group, val in subSDRs.items()],
                                   columns=['level', 'group', 'singleSDR'])
                
                # Write to file
                df = pd.concat([df_sup, df_sub], ignore_index = True, axis = 0)
                SSDR_output_file = (SSDR_output_dir + 'SSDR_' + tree.getGeneName()) +'.csv'                
                df.to_csv(SSDR_output_file, index = False, header = True,na_rep = "NULL")  # Test
                
                ind +=1
                

        except Exception: 
           
            logger.exception(f"File disrupted: {dist_file} ")    
            with open(output_unprocessed, 'w+') as f: 
                write = csv.writer(f) 
                write.writerow(dist_file)
                
            pass

       
    def run_calcSDV(self, random = False):
        """
        Input: 
            input_files: path to files containing cophenetic distance matrices
            group_info: path to group information file
            output_path: path to file SDR values are to be written to. 
            unprocessed: path to file for saving list of unproccessed file upon disruption
            
        Function: 
            Pipeline of function for SDR calculations. 
            Writes SDR values to file with specified gene name.
            
            If function is disrupted, a list of the remaining files to 
            calculate SDR for from the filelist is written to a file. 
        
        Output: 
            Writes SDR to file.         
        """
        input_files = self.config.get('input_files').strip()
        group_info = self.config.get('input_group_info').strip()
        SDVsuper_output = self.config.get('output_SDRsuper').strip()
        SDVsub_output = self.config.get('output_SDRsuper').strip()
        save_unprocessed = self.config.get('output_unprocessed').strip()
            
        file_list = self.make_filelist(input_files)

        ind = 1
        ind_len = len(file_list)
            
        print("Files to procescs:\n", file_list)

        try:
            while file_list:
                dist_file = file_list.pop().strip()
                
                print("File processed: ", dist_file)   # TO DO: convert to log  
                print(f"Number: {ind} / {ind_len}")
                
                tree = treeMetrics()
                tree.setup(dist_file, group_info)
                
                
                if random: 
                    tree.shuffleSampleInfo()
                    
                tree.calcSDV()
                
                supSDV = tree.getSDVsuper()
                subSDV = tree.getSDVsub()
                
                print("Super SDV: ", supSDV)
                print("Sub SDV: ", subSDV)
                
                ind += 1
                
                supSDV = [tree.getGeneName(), supSDV]
                subSDV = [tree.getGeneName(), subSDV]
        
                with open(SDVsuper_output, 'a', newline='') as f:   # write to file    
                    writer = csv.writer(f)
                    writer.writerow(supSDV)
                    
                with open(SDVsub_output, 'a', newline='') as f:   # write to file    
                    writer = csv.writer(f)
                    writer.writerow(subSDV)
                
        except: 
            # Write remaining filelist to file
            print("Disrupted.")
            print("Last file processed:", dist_file)
            
            with open(save_unprocessed, 'a') as f: 
                write = csv.writer(f) 
                write.writerow(file_list) 

       
    def run_calcTest(self):
        """
        Test function to check functionality of program.
        """
        return self.run_calcSDR(random=True)
        

    # def run_CalcNonZeroForPop(self):
    #     """
    #     Input: input file for cd files
    #     Output: file path to store values. 
    #         Three columns: gene name, pop and nonZero pop value. 
    #     Return: None
    
    #     """
    #     pass
    #     # # List of files
    #     # file_list = self.make_filelist(input_files)
    
    #     # ind = 0
    #     # print("Files to process:\n", file_list)
        
    #     # try:
    #     #     while file_list:
    #     #         dist_file = file_list.pop().strip()
    #     #         print("File processed: ", dist__file)
    #     #         print("Number: ", ind)
                
    #     #         # make tree
    #     #         tree = treeInfo(dist_file,
    #     #                         pop_info)
        
    #     #         result = tree.calcNonZerosForPops(popType='all', percent = True)
                
    #     #         with open(output_file, 'a') as f:
    #     #             result.to_csv(output_file, mode = 'a', index_label = 'pop', header=f.tell()==0)
    #     #         ind += 1
    #     #         except: 
    #         # # Write remaining filelist to file
    #         # print("Disrupted.")
    #         # print("Last file processed:", dist_file)
            
    #         # with open(save_unprocessed, 'w') as f: 
    #         #     write = csv.writer(f) 
    #         #     write.writerow(file_list) 
    

    def run_func(self):
        if self.func == "calcSDR":
            self.run_calcSDR()
        elif self.func == "calcSDV":
            self.run_calcSDV()
        elif self.func == "calcSingleSDRs":
            self.run_calcSingleSDRs()
        elif self.func == "calcSDRnull":
            self.run_calcSDR(random = True)
        elif self.func == "calcTest":
            return self.run_calcTest()
        else:
            print("Error with func")
            #logger.error("Not a valid function option. Change in main_config file." +
                              #"options are: calcSDR, calcSDV, calcSingleSDRs, calcNullSDR and calcTest")
   
        #logger.info('Session finished successfully: {0}'.format(datetime.now().strftime("%d.%m.%Y_%H.%M")))



if __name__ == '__main__':
    
   pass