# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 07:05:03 2021

@author: norab
"""

from treeHelpers import *
from treeMetrics import treeMetrics
import pandas as pd



"""
    Script contains formattings and calculations to produce files
    that can be read directly by functions in "treeHelpers.py":
        
        - non-zero phydists
        - unique_seq
        - GDR values for super- and sub populations
        - GDR null-distribution values
"""

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

#%% Filter genes to calculate GDR null distributions for 
"""
    REWRITE: 
        
    Select top 1600 genes with regards to SDR. 
    Selection is made for sub-and super level simultaneously, since
    there will be a large overlap and the goal is to make a rough selection
    based on computitional capacity. 
    Largest SDR after selection is still ~ 0.8, which is high. 
    
    
    ----------------------------------------------------------------
    Filtering process: 
        - Many more than the found ones are likely significant.
        - WHat I did: 
                What I wanted: To end up with approximatly 1000 genes to test with
                    the lowest SDRs for both super and sub populations. 
                1. Filtered out the genes ith less than 20 samples having a distance
                to any other sample in the tree. This is 1/3 of the samllest sbu-population group. 
                2. Selected the 1600 genes with lowest SDR, for both super and sub added. 
                3. As this list conatined SDR values for both super and sub, some genes were 
                duplicated in the list. After dropping duplicates, 1055 genes remained, containing the genes with 
                lowest SDRs in an added list of super and sub populations. 
                And since I were to run SDRnull distribution for these genes anyways, I could do both super and sub for
                both genes, regardless of which of the SDR values it was selected for. 
                If I were to select the gene list of 1000 lowest super SDRs AND 1000 lowest sub SDR, 
                I would end up with far more genes, being computitionally too expensive with repsect to time I hade
                to calculate this. 
                So the genes in the SUPER list that end up non-significant are likely selected to the lit
                because of a low SUB - pop SDR. 
        - Interesting: Which ones are not significant?

"""

# Load GDRs
GDRs = load_GDRs()   # works
nz_phydists = load_nz_phydists()  # works
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


#%% Concatenate super and sub null values


GDRnullSuper = pd.read_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_random_SUPER_27.11.2021.csv', names = ['gene', 'GDRnull'])
GDRnullSub = pd.read_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDR_random_SUB_27.11.2021.csv', names=['gene', 'GDRnull'])

# Load
GDRnullSuper = pd.DataFrame(columns = ['gene', 'GDRnull'])
GDRnullSub = pd.DataFrame(columns = ['gene', 'GDRnull'])

SDRnullSuper.dropna(inplace=True)
SDRnullSub.dropna(inplace=True)

# Add level column
GDRnullSuper.insert(0, 'level', 'nullSuper')
GDRnullSub.insert(0, 'level', 'nullSub')
#SDVpsuedo.insert(0, 'level', 'psuedo')

# concatenate
GDRnullAll = pd.concat([GDRnullSuper, GDRnullSub])

# Save
GDRnullAll.to_csv('C:/Users/norab/Master/thesis_data/test_result_data/GDRnull_all_27.11.21.csv', index = False, header = True)



#%% Create concatenated file of super- and sub GDRnull distribution data

"""
    USE FOR p val calc
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


#%% CAN DELETE THIS

# Fjerne feil i sub
p_val_super = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRnull/GDRnull_pval_super.csv',header=0, names = ['gene','pval'])

# Denne er helt på jordet
p_val_sub = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRnull/GDRnull_pval_sub.csv',header=0, names = ['gene','pval'])
super_list = list(p_val_super['gene'])

p_val_sub2 = p_val_sub[~p_val_sub['gene'].isin(super_list)]

# trikse
p_val_sub = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRsignificant/gene_pval_padj_GDRsub.csv',header=0, names = ['gene','pval','pval_adj','GDRsub'])
p_val_sub_store= p_val_sub[['gene', 'pval']]
p_val_sub_store.to_csv('C:/Users/norab/Master/thesis_data/result_data/GDRnull/GDRnull_pval_sub.csv',header=['gene', 'pval'], index=False)


# trikse
p_val_super = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRsignificant/gene_pval_padj_GDRsuper.csv',header=0, names = ['gene','pval','pval_adj','GDRsub'])
p_val_super_store= p_val_sub[['gene', 'pval']]
p_val_super_store.to_csv('C:/Users/norab/Master/thesis_data/result_data/GDRnull/GDRnull_pval_super2.csv',header=['gene', 'pval'], index=False)


if __name__ == '__main__':
    pass

