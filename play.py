# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 13:42:43 2021

@author: norab
"""

import pandas as pd
import numpy as np
import os
from os.path import join, isfile
from src.treeInformation import treeInfo
from src.treeMetrics import treeMetrics
from datetime import datetime
from src.treeRun import RunStuff

test_gene_small = 'C:\\Users\\norab\\MasterDisaster\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
pop_info = '    \n    C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv'

test_files = {'pop_info': pop_info,
           'dist_mat':test_gene_small}

tree = treeMetrics()
tree.setup(test_files['dist_mat'], test_files['pop_info'])
dist_mat = tree.getDistMat()
sample_info = tree.getSampleInfo()
# Understand iteration of upper matrix
dist_mat = dist_mat.to_numpy()

k = test_files.get('pop_info').strip()

if '.csv' in 'E:/Master/jobs/job_calcSDRnull/job_output/unprocessed_genes_30.06.2021_13.36.csv':
    print('yes')
    
    
    

import sys
  
# stdout assigned to a variable
var = sys.stdout
arr = ['geeks', 'for', 'geeks']
  
# printing everything in the same line
for i in arr:
    var.write(i)
  
# printing everything in a new line
for j in arr:
    var.write('\n'+j)
    
    
# Test if this skip_genes thing work: 
    
configFilepath = 'E:/Master/jobs_template/job_input/main_config_template_calcTest.yml'

file_list = RunStuff.make_filelist('E:/Master/cophenetic_dists_47genes/')
read_genes = pd.read_csv('E:/Master/external_runs/data_software_SDRnull_allGenes_x1000/data/job_input/skip_genes.csv')
skip_genes = list(read_genes['gene'])[0:5]
file_list = file_list[0:5]
skip_genes.append('ENSG00000108849___PAHO')
for cd_file in file_list:

    tree = treeMetrics()
    
    if skip_genes:
        
        tree.setGeneName(cd_file)
        gene_name = tree.getGeneName()
        print(gene_name)
        skip_gene = any(gene_name == gene for gene in skip_genes)    
        print(skip_gene)
        
        if skip_gene:
            continue
    
    print("LOOOOL")
            
        # if skip_gene: 
        #     continue

                    