# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 07:05:03 2021

@author: norab

DESCRIPTION: 
    Formatting of result data for easier data inspection and plotting.
    Data output can be read directly by functions in "treeHelpers.py":

"""

from treeHelpers import *
from treeMetrics import treeMetrics
import pandas as pd




#%% Create concatenated file with all GDR values for super- and sub populations

# load result GDR values
GDRsuper = pd.read_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_SUPER_27.11.2021.csv', names = ['gene', 'GDR'])
GDRsub = pd.read_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_SUB_27.11.2021.csv', names = ['gene', 'GDR'])

GDRsuper.dropna(inplace=True)
GDRsub.dropna(inplace=True)

# Add level column
GDRsuper.insert(0, 'level', 'super')
GDRsub.insert(0, 'level', 'sub')

# Concatenate
GDRall = pd.concat([GDRsub, GDRsuper])

# Save
#GDRall.to_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_10genes_all.csv', index = False, header = True)
#GDRsuper.to_csv('C:/Users/norab/Master/data/GDR/GDRsuper.csv', index = False, header = True)
#GDRsub.to_csv('C:/Users/norab/Master/data/GDR/GDRsub.csv', index = False, header = True)

#%% Filter genes to calculate GDR null distributions for 

"""
    
----------------------------------------------------------------
Filtering process: 

- Goal: End up with approximatly 1000 genes with the lowest GDRs across 
        both super and sub populations. 

1. Filtered out the genes with less than 20 samples that exhibited a distance
to any other sample in the tree (since such small population 
outgroups were not of interest in the master project)

2. Selected the 1600 genes with lowest GDR, for both super and sub added. 

3. As this list contained GDR values for both super and sub, some genes were 
duplicated in the list. After dropping duplicates, 1055 genes remained. 


"""

# Load data
GDRs = load_GDRs()   
nz_phydists = load_nz_phydists()  
data_all = nz_phydists.merge(GDRs, on = 'gene')

# Filter genes with less than 20 samples
selection1 = data_all[data_all['totdist'] > 0.016]

# Sort by GDR
selection1.sort_values(by=['GDR'], inplace=True, ignore_index=True)

# Select 1600 genes
selection2 = selection1[selection1.index < 1600]

# Drop genes appearing twice in list (included for both super and sub population)
gene_selection = list(set(list(selection2['gene'])))
gene_selection = pd.DataFrame(gene_selection, columns = ['gene'])
gene_selection = gene_selection.drop_duplicates()

# Save list of selected genes
# gene_selection.to_csv('C:/Users/norab/Master/thesis_data/result_data/GDRnull/GDRnull_1055genes_selection.csv', index = False)


#%% Concatenate super and sub GDRnull values

"""
    Create one file containing all GDRnull values for both super- and sub
    population level
    
"""

GDRnullSuper = pd.read_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_random_SUPER.csv', names = ['gene', 'GDRnull'])
GDRnullSub = pd.read_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_random_SUB.csv', names=['gene', 'GDRnull'])

# Load
GDRnullSuper = pd.DataFrame(columns = ['gene', 'GDRnull'])
GDRnullSub = pd.DataFrame(columns = ['gene', 'GDRnull'])

SDRnullSuper.dropna(inplace=True)
SDRnullSub.dropna(inplace=True)

# Add level column
GDRnullSuper.insert(0, 'level', 'nullSuper')
GDRnullSub.insert(0, 'level', 'nullSub')

# concatenate
GDRnullAll = pd.concat([GDRnullSuper, GDRnullSub])

# Save
GDRnullAll.to_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDRnull_all.csv', index = False, header = True)



#%% Create concatenated file of super- and sub GDRnull distribution data

"""
    Separate GDRnull values for each gene into separate files (for 
    super- and sub population categories respectivly).
    Used for p-val calculation (in R).
"""

# Files produced from GDRnullCalculation, containing all calculated GDRnull values
GDRnullSuper = pd.read_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_random_SUPER_27.11.2021.csv', names = ['gene', 'GDRnull'])
GDRnullSub = pd.read_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_random_SUB_27.11.2021.csv', names = ['gene', 'GDRnull'])

GDRsuper.dropna(inplace=True)
GDRsub.dropna(inplace=True)

GDRnullSuper.index = GDRnullSuper['gene']
GDRnullSub.index = GDRnullSub['gene']


# Check and create list of uniq genes with random values
print(len(set(GDRnullSuper.index))) # Correct. there are 1055 genes in the df.

# Fill gene selection dict with values.
# Extract genes, save GDR null values as values to each gene key
superValues = dict()
for gene in unique:
    values = list(GDRnullSuper[GDRnullSuper.index == gene]['GDRnull'])
    superValues[gene] = values
    
subValues = dict()
for gene in unique:
    values = list(GDRnullSub[GDRnullSub.index == gene]['GDRnull'])
    subValues[gene] = values

# Save GDR null values for each gene to spearate file    
for key, val in superValues.items():
    df = pd.DataFrame(val)
    filename = 'C:/Users/norab/Master/thesis_data/test_result_data/super/GDRnullSuper_' + key + '.csv'
    df.to_csv(filename, index = False, header = ['GDRnull'])

for key, val in subValues.items():
    df = pd.DataFrame(val)
    filename = 'C:/Users/norab/Master/data/GDRnull/refined_values/sub/GDRnullSub_' + key + '.csv'
    df.to_csv(filename, index = False, header = ['GDRnull'])

# Results: two folders (super and sub) with 1055 files in each, containing 
# GDRnull values for each gene respectivly. 


