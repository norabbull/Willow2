# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 06:45:58 2021

@author: norab
"""

from inspection.inspectionHelpers import *
from src.treeInformation import treeInfo
import pandas as pd
import os
from os.path import isfile, join
import re
import shutil
import numpy as np


#%% Helpers

def make_filelist(input_files):

    if isinstance(input_files, str):     
        files = [join(input_files, f) for f in os.listdir(input_files) 
                     if isfile(join(input_files, f))]
    return files

def geneName(dist_mat_file):
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
    gene_name = re.sub('___CopD.csv$','', subName)
    
    return gene_name

#%% Load data

test_genes = ['C:/Users/norab/Master/data/real_tree_data/dist_mat_subset/ENSG00000000419___DPM1___CopD.csv',
             'C:/Users/norab/Master/data/real_tree_data/dist_mat_subset/ENSG00000000938___FGR___CopD.csv',
             'C:/Users/norab/Master/data/real_tree_data/dist_mat_subset/ENSG00000000971___CFAH___CopD.csv']

pop_info = 'C:/Users/norab/Master/data/real_tree_data/phydist_population_classes.tsv'
    

#%% make test dataframes

# =============================================================================
#  Make distance matrix subset, 10 x 10 samples
# =============================================================================

tree_FGR = treeInfo()
tree_FGR.setup(test_genes[1], pop_info)
mat_FGR = tree.getDistMat()

submat_FGR = mat.iloc[0:10,0:10]
submat_FGR.to_csv('C:/Users/norab/Master/data/real_tree_data/dist_mat_test/FGR_10x10.csv')

# =============================================================================
# Make distance matrix subset, 5 x 5 samples
# =============================================================================

# TO DO: 
# Add MIO gene to subset on local machine
tree_MIO = treeInfo()
tree_MIO.setup(test_genes[0], pop_info)   # Change this to correct
mat_MIO = tree_MIO.getDistMat()

# Create matrix with values good for testing. Not reeal values. 
submat_MIO = mat_MIO.iloc[132:137,132:137]
submat_MIO['EUR___GBR___HG00236']['EUR___GBR___HG00159'] = 0.005
submat_MIO['EUR___GBR___HG00159']['EUR___GBR___HG00236'] = 0.005
submat_MIO['AMR___MXL___NA19758']['AMR___MXL___NA19795'] = 0.003
submat_MIO['AMR___MXL___NA19795']['AMR___MXL___NA19758'] = 0.003

submat_MIO.to_csv('C:/Users/norab/Master/data/real_tree_data/dist_mat_test/MIO_5x5.csv')

# Check if equal: 
submat4 = pd.read_csv('C:/Users/norab/Master/data/real_tree_data/dist_mat_test/MIO_5x5.csv', index_col = 0)
(submat3 == submat4).eq(True).all().all()



#%% Prosess and extract genelists 

# =============================================================================
#  TEMPLATE: Make list of unprocces files from value-file (eg. SDRs)
# =============================================================================

processed_genes = pd.read_csv('path_to_file', index_col = 0, header = None)
processed_genes = pd.Series(processed_genes.index)

# all genes
all_genes = make_filelist('E:\Master\cophenetic_dists')
unprocessed_genes = all_genes.copy()

for file in all_genes: 
    for gene in processed_genes: 
        if gene in file: 
            try:
                unprocessed_genes.remove(file)
            except:
                print(f"{file} not removed. ")
            

unprocessed_genes = pd.Series(unprocessed_genes)
unprocessed_genes.to_csv('save_path', index = False)

# =============================================================================
# TEMPLATE: Multiple files to put into same list
# =============================================================================

all_files = make_filelist('folder_path')
all_files = [set(f for f in all_files if 'some_common_filename' in f)]

# =============================================================================
# 
# =============================================================================
# Interruption 09.06.14.40
processed_genes = pd.read_csv('E:\Master\current_run\SDRsuper_runDate08.06.2021_15.57.csv', index_col = 0)
processed_genes = pd.Series(processed_genes.index) 

# Add ful path
processed_genes.to_csv('E:\\Master\\current_run\\processed_genes_09.06_14.40.csv', index = False)

# all genes
all_genes = make_filelist('E:\Master\cophenetic_dists')
unprocessed_genes = all_genes.copy()

for file in all_genes: 
    name = geneName(file)
    for g in processed_genes: 
        if name in g: 
            unprocessed_genes.remove(file)

unprocessed_genes = pd.Series(unprocessed_genes)
unprocessed_genes.to_csv('E:\\Master\\current_run\\unprocessed_genes_09.06_14.40.csv', index = False)

# 10.06, 13.16 ----------------------------------------------------------------
processed_genes = list(pd.read_csv('E:\Master\current_run\SDRsub_runDate08.06.2021_15.57.csv', index_col = 0, header = None).index)
processed_genes2 = list(pd.read_csv('E:\Master\current_run\SDRsub_runDate09.06.2021_16.57.csv', index_col = 0, header = None).index)
processed_genes3 = list(pd.read_csv('E:\Master\current_run\SDRsub_runDate10.06.2021_03.54.csv', index_col = 0, header = None).index)
processed_genes4 = list(pd.read_csv('E:\Master\current_run\SDRsub_runDate10.06.2021_13.46.csv', index_col = 0, header = None).index)

all_processed_genes = list(set(processed_genes + processed_genes2 + processed_genes3 + processed_genes4))
#all_processed_genes = processed_genes + processed_genes2 + processed_genes3 + processed_genes4
# Add ful path

pd.Series(all_processed_genes).to_csv('E:\\Master\\current_run\\processed_genes_10.06.2021_13.35.csv', index = False)

# all genes
all_genes = list(set(make_filelist('E:\Master\cophenetic_dists')))
unprocessed_genes = all_genes.copy()

for file in all_genes: 
    #name = geneName(file)
    for gene in all_processed_genes: 
        if gene in file: 
            try:
                unprocessed_genes.remove(file)
            except:
                print(f"{file} not removed. ")
            

ug = pd.Series(unprocessed_genes)
ug.to_csv('E:\\Master\\current_run\\unprocessed_genes_11.06.21_11.40.csv', index = False)

# =============================================================================
# Make list of unprocessed files from SDR file, null dist
# =============================================================================
# Interruption 09.06.14.40
processed_genes = pd.read_csv('E:\\Master\\jobs\\job_calcSDRnull\\job_output\\nullDistSDRsuper_30.06.2021_13.36.csv', index_col = 0, header = None)
processed_genes = pd.Series(processed_genes.index) 

# all genes
all_genes = make_filelist('E:\Master\cophenetic_dists')
unprocessed_genes = all_genes.copy()

for file in all_genes: 
    #name = geneName(file)
    for gene in processed_genes: 
        if gene in file: 
            try:
                unprocessed_genes.remove(file)
            except:
                print(f"{file} not removed. ")
            

unprocessed_genes = pd.Series(unprocessed_genes)
unprocessed_genes.to_csv('E:\\Master\\jobs\\job_calcSDRnull\\job_output\\unprocessed_genes_30.06.2021_13.36.csv', index = False)


### SAME ######################################################################

# Interruption, 07.07.32
processed_genes = pd.read_csv('E:\\Master\\jobs\\job_calcSDRnull\\job_output\\nullDistSDRsuper_06.07.2021_12.10.csv', index_col = 0, header = None)
processed_genes = pd.Series(processed_genes.index) 
# Add ful path

#processed_genes.to_csv('E:\\Master\\current_run\\processed_genes_09.06_14.40.csv', index = False)

# all genes
all_genes = make_filelist('E:\Master\cophenetic_dists')
unprocessed_genes = all_genes.copy()

for file in all_genes: 
    #name = geneName(file)
    for gene in processed_genes: 
        if gene in file: 
            try:
                unprocessed_genes.remove(file)
            except:
                print(f"{file} not removed. ")
            

unprocessed_genes = pd.Series(unprocessed_genes)
unprocessed_genes.to_csv('E:\\Master\\jobs\\job_calcSDRnull\\job_output\\unprocessed_genes_30.06.2021_13.36.csv', index = False)

# =============================================================================
# Make list of unprocessed files from unprocessed files list file
# =============================================================================

redhood_input_files_continue2 = 'C:\\Users\\norab\\Master\\data\\runstop_save\\unprocessed_files_SDRcalc_09.06.2021_16.57.csv' 
liste = pd.read_csv(redhood_input_files_continue2 , header = None, sep = ",",).transpose()
liste.to_csv('C:\\Users\\norab\\Master\\data\\runstop_save\\unprocessed_files_SDRcalc_09.06.2021_16.57_update.csv' , 
             index = False,
             header = False)


#%% Process datafiles to prepare for inspections


# =============================================================================
# Merge SDR calculation files to one
# =============================================================================

# Extract SDRsub and SDRsuper file paths


all_files = make_filelist('E:\\Master\\data\\SDRnull')
SDRsuper_files = [f for f in all_files if 'SDRsuper' in f]
SDRsub_files = [f for f in all_files if 'SDRsub' in f]

# Load
SDRsuper = pd.DataFrame(columns = ['gene', 'val'])
SDRsub = pd.DataFrame(columns = ['gene', 'val'])
for f in SDRsuper_files: 
    SDRsuper = SDRsuper.append(pd.read_csv(f, names = ['gene', 'val']))
for f in SDRsub_files: 
    SDRsub = SDRsub.append(pd.read_csv(f, names = ['gene', 'val']))

SDRsuper.drop_duplicates(inplace=True)
SDRsub.drop_duplicates(inplace=True)

# Save
SDRsuper.to_csv('E:\\Master\\current_run\\SDRsuper_all.csv', index = False, header = False)
SDRsub.to_csv('E:\\Master\\current_run\\SDRsub_all.csv', index = False, header = False)

# =============================================================================
# Merge SDR null dist calc files to one
# =============================================================================

SDRsuper_files = [f for f in all_files if 'SDRsuper' in f]
SDRsub_files = [f for f in all_files if 'SDRsub' in f]

# Load
SDRsuper = pd.DataFrame(columns = ['gene', 'val'])
SDRsub = pd.DataFrame(columns = ['gene', 'val'])
for f in SDRsuper_files: 
    SDRsuper = SDRsuper.append(pd.read_csv(f, names = ['gene', 'val']))
for f in SDRsub_files: 
    SDRsub = SDRsub.append(pd.read_csv(f, names = ['gene', 'val']))

SDRsuper.dropna(inplace=True)
SDRsub.dropna(inplace=True)

# Save
SDRsuper.to_csv('E:\\Master\\data\\SDRnull\\nullDistSDRsuper_all_07.07.21.csv', index = False, header = False)
SDRsub.to_csv('E:\\Master\\data\\SDRnull\\nullDistSDRsub_all_07.07.21.csv', index = False, header = False)

# =============================================================================
# Filter genes with SDR <= 0.7, 0.6 and 0.5 - NOT RUN YET CUZ SDV IS RUNNING
# =============================================================================

# Super:

SDRsuper_all = pd.read_csv('C:\\Users\\norab\\Master\\data\\SDR\\SDRsuper_all.csv',
                           index = False, header = False, names = ['gene', 'val'])
SDRsub_all = pd.read_csv('C:\\Users\\norab\\Master\\data\\SDR\\SDRsub_all.csv',
                         index = False, header = False, names = ['gene', 'val'])

SDRsuper_07 = SDRsuper_all[SDRsuper_all['val'] <= 0.7]
SDRsub_07 = SDRsub_all[SDRsub_all['val'] <= 0.7]

SDRsuper_06 = SDRsuper_all[SDRsuper_all['val'] <= 0.6]
SDRsub_06 = SDRsub_all[SDRsub_all['val'] <= 0.6]

SDRsuper_05 = SDRsuper_all[SDRsuper_all['val'] <= 0.5]
SDRsub_05 = SDRsub_all[SDRsub_all['val'] <= 0.5]

# Save 
SDRsuper_07.to_csv('C:\\Users\\norab\\Master\\data\\SDR\\SDRsuper_07.csv', 
                   index = False, header = False)
SDRsub_07.to_csv('C:\\Users\\norab\\Master\\data\\SDR\\SDRsub_07.csv', 
                   index = False, header = False)
SDRsuper_06.to_csv('C:\\Users\\norab\\Master\\data\\SDR\\SDRsuper_06.csv', 
                   index = False, header = False)
SDRsub_06.to_csv('C:\\Users\\norab\\Master\\data\\SDR\\SDRsub_06.csv', 
                   index = False, header = False)
SDRsuper_05.to_csv('C:\\Users\\norab\\Master\\data\\SDR\\SDRsuper_05.csv', 
                   index = False, header = False)
SDRsub_05.to_csv('C:\\Users\\norab\\Master\\data\\SDR\\SDRsub_05.csv', 
                   index = False, header = False)

#%%  Make data files ready for inspecton 

"""
From .csv file on the form: 
    
    ENSG00000288516___KAPCA,0.0934
    ENSG00000288513___OR4F4,0.0
    ENSG00000288508___ISLR,0.1869
    ENSG00000288499___FA20C,0.0053
    ENSG00000288490___GALT9,0.3585
    ...

To .csv file on the form: 
    
    level,gene,val
    sub,ENSG00000288516___KAPCA,0.0934
    sub,ENSG00000288513___OR4F4,0.0
    sub,ENSG00000288508___ISLR,0.1869
    sub,ENSG00000288499___FA20C,0.0053
    sub,ENSG00000288490___GALT9,0.3585
    ...
    super,ENSG00000101199___ARFG1,0.0751
    super,ENSG00000101191___DIDO1,0.001
    super,ENSG00000101189___MRGBP,0.0462
    ...
    
"""

# =============================================================================
# SDR
# =============================================================================

all_files = make_filelist('C:/Users/norab/Master/data/SDR')
SDRsuper_files = [f for f in all_files if 'SDRsuper' in f]
SDRsub_files = [f for f in all_files if 'SDRsub' in f]

# Load
SDRsuper = pd.DataFrame(columns = ['gene', 'SDR'])
SDRsub = pd.DataFrame(columns = ['gene', 'SDR'])
for f in SDRsuper_files: 
    SDRsuper = SDRsuper.append(pd.read_csv(f, names = ['gene', 'SDR']))
for f in SDRsub_files: 
    SDRsub = SDRsub.append(pd.read_csv(f, names = ['gene', 'SDR']))

SDRsuper.dropna(inplace=True)
SDRsub.dropna(inplace=True)

# Add level column
SDRsuper.insert(0, 'level', 'super')
SDRsub.insert(0, 'level', 'sub')
#SDRpsuedo.insert(0, 'level', 'psuedo')

SDRall = pd.concat([SDRsub, SDRsuper])

# Save
SDRsuper.to_csv('C:/Users/norab/Master/data/SDR/SDRsuper.csv', index = False, header = True)
SDRsub.to_csv('C:/Users/norab/Master/data/SDR/SDRsub.csv', index = False, header = True)
SDRall.to_csv('C:/Users/norab/Master/data/SDR/SDR_all.csv', index = False, header = True)


# =============================================================================
# SDV
# =============================================================================

all_files = make_filelist('C:/Users/norab/Master/data/SDV')
SDVsuper_files = [f for f in all_files if 'SDVsuper' in f]
SDVsub_files = [f for f in all_files if 'SDVsub' in f]

# Load
SDVsuper = pd.DataFrame(columns = ['gene', 'SDV'])
SDVsub = pd.DataFrame(columns = ['gene', 'SDV'])
for f in SDVsuper_files: 
    SDVsuper = SDVsuper.append(pd.read_csv(f, names = ['gene', 'SDV']))
for f in SDVsub_files: 
    SDVsub = SDVsub.append(pd.read_csv(f, names = ['gene', 'SDV']))

SDVsuper.dropna(inplace=True)
SDVsub.dropna(inplace=True)

# Add level column
SDVsuper.insert(0, 'level', 'super')
SDVsub.insert(0, 'level', 'sub')
#SDVpsuedo.insert(0, 'level', 'psuedo')

SDVall = pd.concat([SDVsub, SDVsuper])

# Save
SDVsuper.to_csv('C:/Users/norab/Master/data/SDV/SDVsuper.csv', index = False, header = True)
SDVsub.to_csv('C:/Users/norab/Master/data/SDV/SDVsub.csv', index = False, header = True)
SDVall.to_csv('C:/Users/norab/Master/data/SDV/SDV_all.csv', index = False, header = True)

# =============================================================================
# Null dist
# =============================================================================

all_files = make_filelist('C:/Users/norab/Master/data/SDRnull/runs/run_47genes')
SDRsuper_files = [f for f in all_files if 'SDRsuper' in f]
SDRsub_files = [f for f in all_files if 'SDRsub' in f]

# Load
SDRsuper = pd.DataFrame(columns = ['gene', 'SDR'])
SDRsub = pd.DataFrame(columns = ['gene', 'SDR'])
for f in SDRsuper_files: 
    SDRsuper = SDRsuper.append(pd.read_csv(f, names = ['gene', 'SDR']))
for f in SDRsub_files: 
    SDRsub = SDRsub.append(pd.read_csv(f, names = ['gene', 'SDR']))

SDRsuper.dropna(inplace=True)
SDRsub.dropna(inplace=True)

# Add level column
SDRsuper.insert(0, 'level', 'nullSuper')
SDRsub.insert(0, 'level', 'nullSub')
#SDRpsuedo.insert(0, 'level', 'psuedo')

SDRall = pd.concat([SDRsub, SDRsuper])

# Save
SDRsuper.to_csv('C:/Users/norab/Master/data/SDRnull/all/SDRnullSuper_47genes.csv', index = False, header = True)
SDRsub.to_csv('C:/Users/norab/Master/data/SDRnull/all/SDRnullSub_47genes.csv', index = False, header = True)
SDRall.to_csv('C:/Users/norab/Master/data/SDRnull/all/SDRnullAll_47genes.csv', index = False, header = True)

# =============================================================================
# Extract value for interspaced null Distribution 
# =============================================================================

SDRnullSuper_files = ['C:/Users/norab/Master/data/SDRnull/all/nullDistSDRsuper_47genes_12.07.2021_12.38.csv']
SDRnullSub_files = ['C:/Users/norab/Master/data/SDRnull/all/nullDistSDRsub_47genes_12.07.2021_12.38.csv']

# Load
SDRnullSuper = pd.DataFrame(columns = ['gene', 'SDR'])
SDRnullSub = pd.DataFrame(columns = ['gene', 'SDR'])
for f in SDRnullSuper_files: 
    SDRnullSuper = SDRnullSuper.append(pd.read_csv(f, names = ['gene', 'SDR']))
for f in SDRnullSub_files: 
    SDRnullSub = SDRnullSub.append(pd.read_csv(f, names = ['gene', 'SDR']))

SDRnullSuper.dropna(inplace=True)
SDRnullSub.dropna(inplace=True)

# Add level column
SDRnullSuper.insert(0, 'level', 'nullSuper')
SDRnullSub.insert(0, 'level', 'nullSub')
#SDVpsuedo.insert(0, 'level', 'psuedo')

SDRnullAll = pd.concat([SDRnullSuper, SDRnullSub])

# Save

SDRnullAll.to_csv('C:/Users/norab/Master/data/SDRnull/all/SDRnull47_all_12.07.2021.csv', index = False, header = True)


# =============================================================================
# Separate genes for repeated SDR calculation
# =============================================================================

LIN37_super = SDRnullSuper[SDRnullSuper['gene'] == 'ENSG00000267796___LIN37']
KRA22_super = SDRnullSuper[SDRnullSuper['gene'] == 'ENSG00000214518___KRA22']
DAXX_super = SDRnullSuper[SDRnullSuper['gene'] == 'ENSG00000227046___DAXX']

LIN37_sub = SDRnullSub[SDRnullSub['gene'] == 'ENSG00000267796___LIN37']
KRA22_sub = SDRnullSub[SDRnullSub['gene'] == 'ENSG00000214518___KRA22']
DAXX_sub = SDRnullSub[SDRnullSub['gene'] == 'ENSG00000227046___DAXX']

LIN37_SDRnull_all = pd.concat([LIN37_super, LIN37_sub])
KRA22_SDRnull_all = pd.concat([KRA22_super, KRA22_sub])
DAXX_SDRnull_all = pd.concat([DAXX_super, DAXX_sub])

LIN37_SDRnull_all['level'].replace('nullSuper', 'nullSuperLIN37', inplace = True)
KRA22_SDRnull_all['level'].replace('nullSuper', 'nullSuperKRA22', inplace = True)
DAXX_SDRnull_all['level'].replace('nullSuper', 'nullSuperDAXX', inplace = True)
LIN37_SDRnull_all['level'].replace('nullSub', 'nullSubLIN37', inplace = True)
KRA22_SDRnull_all['level'].replace('nullSub', 'nullSubKRA22', inplace = True)
DAXX_SDRnull_all['level'].replace('nullSub', 'nullSubDAXX', inplace = True)

LIN37_SDRnull_all.to_csv('C:/Users/norab/Master/data/SDRnull/all/LIN37_nullSDR_all.csv', index = False, header = True)
KRA22_SDRnull_all.to_csv('C:/Users/norab/Master/data/SDRnull/all/KRA22_nullSDR_all.csv', index = False, header = True)
DAXX_SDRnull_all.to_csv('C:/Users/norab/Master/data/SDRnull/all/DAXX_nullSDR_all.csv', index = False, header = True)

# Make gene as level
LIN37_SDRnull_all['level'].replace('nullSuper', 'nullSuperLIN37', inplace = True)
LIN37_SDRnull_all['level'].replace('nullSuper', 'nullSuperLIN37', inplace = True)
LIN37_SDRnull_all['level'].replace('nullSuper', 'nullSuperLIN37', inplace = True)


# =============================================================================
# Load data
# =============================================================================

# (Removed read-in...)

nullSDRsuper_all.sort_values(by = 'SDR', inplace = True, ignore_index = True)
nullSDRsub_all.sort_values(by = 'SDR', inplace = True, ignore_index = True)

# =============================================================================
# Get 20 values from each distribution, interspaced through whole range
# =============================================================================

idx = np.linspace(0, len(nullSDRsuper_all)-1, 20)
superGenes20 = nullSDRsuper_all.iloc[idx]
subGenes20 = nullSDRsub_all.iloc[idx]

# =============================================================================
# Get 20 values from each distribution, with SDR < 0.9
# =============================================================================

# Make 0.9 cutoff
nullSDRsuper09 = nullSDRsuper_all[nullSDRsuper_all['SDR'] <= 0.8]
nullSDRsub09 = nullSDRsub_all[nullSDRsub_all['SDR'] <= 0.8]

# Select 10 from each
idx = np.linspace(1, len(nullSDRsuper09)-1, 5)
superGenes5 = nullSDRsuper_all.iloc[idx]
subGenes5 = nullSDRsub_all.iloc[idx]

allSuperGenes = superGenes20.append(superGenes5)
allSubGenes = subGenes20.append(subGenes5)

# Extract gene list: 

allSuperGenes = allSuperGenes['gene']    
allSubGenes = allSubGenes['gene']    

allGenes = allSuperGenes.append(allSubGenes)

allGenes.drop_duplicates(inplace = True)

# Make folder containing cophenetic dists files with those genes:
    
all_files = make_filelist('E:\\Master\\cophenetic_dists')

allGenes_fullpath = pd.DataFrame([f for f in all_files for g in allGenes if g in f])

# save list of genes: 

allGenes_fullpath.to_csv('E:\\Master\\data\\SDRnull\\other\\nullSDR_47genes.csv', index = False, header = False)



# =============================================================================
# Get 70 genes with lowest SDR
# =============================================================================

SDRs_super = load_SDRs(level='super').sort_values(by = 'SDR')
SDRs_sub = load_SDRs(level='sub'.sort_values(by= 'SDR')
SDRs_super70 = SDRs_super.iloc[0:70]
SDRs_sub70 = SDRs_sub.iloc[0:70]
genes = SDRs_super70['gene'].append(SDRs_sub70['gene'])
genes.drop_duplicates(inplace=True)
genes = genes.tolist()

all_files = make_filelist('E:\\Master\\cophenetic_dists')
genes_fullpath = pd.DataFrame([f for f in all_files for g in genes if g in f])
genes_fullpath.to_csv('E:\\Master\\data\\SDRnull\\other\\nullSDR_93lowestSDRGenes.csv', index = False, header = False)

# Movie files from a folder to another

for file in all_files:
    for gene in genes:
        if gene in file: 
            shutil.copy(file, 'E:\\Master\\external_runs\\nora_data_software2\\data\\job_input\\cophenetic_dists_93genes')
            
            
#%% Null distribution calculations

import math
prob = 0.05
var_p = prob * (1 - prob)
num_perm_samples = 100      # 70 = empirical p val = 0.1
var_p_n = var_p / num_perm_samples
sd_p_n = math.sqrt(var_p_n)

# 0.025 quantile
max(0, (prob-1.96*sd_p_n))

# 0.0975 quantile of true empirical p-val
prob + 1.96 * sd_p_n



# Find genes with already enough values: 

SDRs = load_SDRs()
# SDVs = load_SDVs()
SDRnullAll = load_SDRnull()
SDRnull93 = load_SDRnull(gene_set="93genes")
SDRnull47 = load_SDRnull(gene_set="47genes")
# SDRnull47super = SDRnull47[SDRnull47['level'] == 'super']
# SDRnull47sub = SDRnull47[SDRnull47['level'] == 'sub']

# Merge all SDRnull data
SDRnullTotal = SDRnullAll.append(SDRnull47)
SDRnullTotal = SDRnullTotal.append(SDRnull93)

SDR_counts = SDRnullTotal['gene'].value_counts().to_frame()
SDR_counts['gene'] = SDR_counts.Index
enough_vals = SDR_counts[SDR_counts['gene'] > 8 ]
enough_vals_genes = enough_vals.index.to_list()


file = 'E:/Master/external_runs/nora_data_software/data/job_input/cophenetic_dists_47genes/ENSG00000110628___S22AI___CopD.csv'

for gene in enough_vals_genes:
    print(gene in file)

for 
skip_gene = any(gene in dist_file for gene in genes)

input_files =  'E:\\Master\\cophenetic_dists'
file_list = make_filelist(input_files)
file_list_stripped = [f.strip() for f in file_list]


gene_names = [geneName(gene) for gene in file_list_stripped]
skip_genes = [file_list_stripped[i] for i in np.arange(0,200,10)]

gene_name == gene





