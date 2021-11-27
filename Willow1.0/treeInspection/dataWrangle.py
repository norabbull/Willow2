# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 07:05:03 2021

@author: norab
"""

from treeHelpers import *
from treeMetrics import treeMetrics


"""
    Script contains formattings and calculations to produce files
    that can be read directly by functions in "treeHelpers.py":
        
        - non-zero phydists
        - unique_seq
        - GDR values for super- and sub populations
        - GDR null-distribution values
"""

#%% Create non-zero distance values file 

"""
    Non-zero phydist values    
    Create file containing total amount of non-zero values in each matrix

"""

# folder with all phydist files 
phydist_files = 'C:/Users/norab/Master/thesis_data/test_data/phydists_10samples'

# calculate percentage of non-zero values in phydist matrices that is non-zero
test = calc_nonZero_phydists(phydist_files)    # function defined above

# save to file 
saveto = 'C:/Users/norab/Master/thesis_data/meta_data/test_nonzero_phydists.csv'
test.to_csv(saveto, index = False, header = True)

#%% Create concatenated file with all GDR values for super- and sub populations

# load calculated GDR values

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
GDRall.to_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_10genes_all.csv', index = False, header = True)
#GDRsuper.to_csv('C:/Users/norab/Master/data/GDR/GDRsuper.csv', index = False, header = True)
#GDRsub.to_csv('C:/Users/norab/Master/data/GDR/GDRsub.csv', index = False, header = True)

#%% Create concatenated file of super- and sub GDRnull distribution data

"""
    Separate random GDR values fro each gene into separate files for 
    super- and sub population categories respectivly.
    Ie. ranodm values for each gene is stored in separate file. 
"""

# Files produced from GDRnullCalculation, containing all calculated GDRnull values
GDRnullSuper = pd.read_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_random_SUPER_27.11.2021.csv', names = ['gene', 'GDRnull'])
GDRnullSub = pd.read_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_random_SUB_27.11.2021.csv', names = ['gene', 'GDRnull'])

GDRsuper.dropna(inplace=True)
GDRsub.dropna(inplace=True)

GDRnullSuper.index = GDRnullSuper['gene']
GDRnullSub.index = GDRnullSub['gene']


# Check and create list of uniq genes with random values
unique = list(set(GDRnullSuper.index)) # Correct. there are 1055 genes in the df.

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

#%% Data Wrangle functions

def calc_nonZero_phydists(folder):
    """
    Input:    folder with phylgenetic distance files
    Function: calculate total amount of non-zero values for each file in folder
    Return:   pandas dataframe with values for all genes
        
    
    """
    
    # make list of files
    phydist_files = make_filelist(folder)

    # initiate df 
    all_nz_phydists = pd.DataFrame(columns = ['gene', 'nz_phydists'])
    
    # process files in filelist
    for file in phydist_files: 
        phydists = load_phydists(file)
        nz_phydists = treeMetrics.calcNonZeroPhydists(phydists)
        geneID = geneName(file)
        all_nz_phydists = all_nz_phydists.append({'gene': geneID, 'nz_phydists':nz_phydists}, ignore_index = True)
        
    
    return all_nz_phydists
    


if __name__ == '__main__':
    pass

