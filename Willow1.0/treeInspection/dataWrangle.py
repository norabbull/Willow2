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
#SDRsuper.to_csv('C:/Users/norab/Master/data/SDR/SDRsuper.csv', index = False, header = True)
#SDRsub.to_csv('C:/Users/norab/Master/data/SDR/SDRsub.csv', index = False, header = True)




#%% Data Wranglew functions

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

