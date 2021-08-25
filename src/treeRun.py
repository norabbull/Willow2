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
from treeMetrics import treeMetrics
import logging
import logging.config

    
class RunStuff:
    
    def __init__(self, config_file):
        
        self.config_file = config_file
        
        for key, val in self.config_file.items():
            if 'datetime' in str(val): 
                self.config_file[key] = val.replace('datetime', datetime.now().strftime("%d.%m.%Y_%H.%M"))
        
        # Set function to run
        self.func = self.config_file.get('func').strip()
        self.input_folder = self.config_file.get('input_folder').strip()
        self.output_folder = self.config_file.get('output_folder').strip()
        self.skip_genes = self.output_folder + self.config_file.get('skip_genes').strip()
        
        logger.info("RUN INITIATED SUCCESSFULLY")

    def make_filelist(self, input_files):

        if isinstance(input_files, str):
            if '.csv' in input_files: 
                files = pd.read_csv(input_files, header = None)
                files = files[0].values.tolist()
            else: 
                files = [join(input_files, f) for f in os.listdir(input_files) 
                             if isfile(join(input_files, f))]
                files = [f.strip() for f in files]
                
        return files
        
    def run_calcSDR(self, skip_genes = None, random=False):
        """
        random: If groups should be random or not. Used to calculate null distribution.
        skip_genes: list of gene names to be skiped from calculation. 

        """
        
        try:
            
            input_cd_folder = self.config_file.get('input_cd_folder').strip()
            group_info = self.input_folder + self.config_file.get('input_group_info').strip()
            
            output_super = self.output_folder + self.config_file.get('output_super').strip()
            output_sub = self.output_folder + self.config_file.get('output_sub').strip()
            output_unprocessed = self.output_folder + self.config_file.get('output_unprocessed').strip()
                     
            file_list = self.make_filelist(input_cd_folder)

            ind = 1
            ind_len = len(file_list)
                
            print("Files to process:\n", file_list)
      
            for cd_file in file_list:
                try:         
                  #  dist_file = dist_file.strip()
                                  
                    print("File processed: ", cd_file)
                    print(f"Number: {ind} / {ind_len}")
                    ind += 1
                    
                    tree = treeMetrics()
                    
                    if skip_genes:
                        
                        tree.setGeneName(cd_file)
                        gene_name = tree.getGeneName()
                        skip_gene = any(gene_name == gene for gene in skip_genes)    
                    
                    tree.setup(cd_file, group_info)
                    
                    
                    if random: 
                        tree.shuffleSampleInfo()
                    
                    tree.calcSDR()
                    
                    supSDR = tree.getSDRsuper()
                    subSDR = tree.getSDRsub()
                    
                    print("Super SDR: ", supSDR)
                    print("Sub SDR: ", subSDR)
                    print("1st sample: ", tree.getSampleInfo()[0])
                    
                    supSDR = [tree.getGeneName(), supSDR]
                    subSDR = [tree.getGeneName(), subSDR]
            
                    with open(output_super, 'a', newline='') as f:   # write to file    
                        writer = csv.writer(f)
                        writer.writerow(supSDR)
                        
                    with open(output_sub, 'a', newline='') as f:   # write to file    
                        writer = csv.writer(f)
                        writer.writerow(subSDR)
    
                except Exception: 
                   
                    #self.logger.exception(f"File disrupted: {dist_file} ")    
                    with open(output_unprocessed, 'w+') as f: 
                        write = csv.writer(f) 
                        write.writerow(cd_file)
                        
                    pass
                
        except Exception as e: 
            self.logger.exception(e)

       
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

        input_files = self.config_file.get('input_files')
        group_info = self.config_file.get('input_group_info')
        SSDR_output_dir = self.config_file.get('output_SSDR')
        output_unprocessed = self.config_file.get('output_unprocessed')
               
        file_list = self.make_filelist(input_files)

        ind = 1
        ind_len = len(file_list)
            
        logger.info("Files to process:\n", file_list)

        try:
            while file_list:
                dist_file = file_list.pop().strip()
                
                logger.info("File processed: ", dist_file)   # TO DO: convert to log  
                logger.info(f"Number: {ind} / {ind_len}")
                
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
        input_files = self.config_file.get('input_files').strip()
        group_info = self.config_file.get('input_group_info').strip()
        SDVsuper_output = self.config_file.get('output_SDRsuper').strip()
        SDVsub_output = self.config_file.get('output_SDRsuper').strip()
        save_unprocessed = self.config_file.get('output_unprocessed').strip()
            
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
        try:
            logging.basicConfig(filename = self.config_file.get('output_folder') + self.config_file.get('output_log'), 
                                level=logging.DEBUG,
                                format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                                filemode = 'w+')
    
            input_files = self.config_file.get('input_files').strip()
            group_info = self.config_file.get('input_group_info').strip()
            output_folder = self.config_file.get('output_folder').strip()
            measure_output = self.config_file.get('output_SDRsuper').strip()
            save_unprocessed = self.config_file.get('output_unprocessed').strip()
                
            file_list = self.make_filelist(input_files)
            
            print("File list: ")
            print(file_list)
    
            dist_file = file_list.pop().strip()
            print("File processed: ", dist_file)   # TO DO: convert to log  
    
            tree = treeMetrics()
            tree.setup(dist_file, group_info)
            print("Gene name: ", tree.getGeneName())
            print("Sample info: ", tree.getSampleInfo())
            
            rand_num = [15.05]
            
            with open(output_folder + measure_output, 'a', newline='') as f:   # write to file    
                writer = csv.writer(f)
                writer.writerow(rand_num)
            
            self.logger.info("Successful run.")
            
        except Exception as e:
            self.logger.exception(e)
            
        

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
    

    def main(self):

        # logging.info('New session...', extra={'function':self.func})
        
        if self.skip_genes: 
            sg = pd.read_csv(self.skip_genes)
            self.skip_genes = list(sg['gene'])
        if self.func == "calcSDR":
            self.run_calcSDR(self.skip_genes)
        elif self.func == "calcSDV":
            self.run_calcSDV()
        elif self.func == "calcSingleSDRs":
            self.run_calcSingleSDRs()
        elif self.func == "calcSDRnull":
            num_rand_trees = int(self.config_file.get('num_rand_trees'))                
            for _ in range(num_rand_trees):
                self.run_calcSDR(skip_genes = self.skip_genes, random = True)
        elif self.func == "calcTest":
            self.run_calcTest()
        else:
            self.logger.error("Not a valid function option. Change in main_config file." +
                              "options are: calcSDR, calcSDV, calcSingleSDRs, calcNullSDR and calcTest")

        
        # logging.info('Session finished: ', datetime.now().strftime("%d.%m.%Y_%H.%M"))

if __name__ == '__main__':
    
    # Test on lenovo computer
    #configFilepath = 'E:/Master/jobs/debug_calcSDRnullDist_28.06.yml'
    
    # Config file argument
    configFilepath = sys.argv[1]
    
    # Load config arguments
    configFilepath = configFilepath.strip()
    
    with open(configFilepath, 'r') as c:
        config_file = yaml.safe_load(c)
        
    # Initiate logger
            # Quick fix logger
    file_handler = logging.FileHandler(filename=config_file.get('output_folder') + config_file.get('output_log'))
    stdout_handler = logging.StreamHandler(sys.stdout)
    handlers = [file_handler, stdout_handler]
    
    logging.basicConfig(
        level=logging.DEBUG, 
        format='[%(asctime)s] {%(name)s:%(lineno)d} %(levelname)s - %(message)s',
        handlers=handlers
    )

    logger = logging.getLogger(__name__)
    
    run = RunStuff(config_file)
    run.main()
    
    