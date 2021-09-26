# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 11:55:22 2021

@author: norab
"""


from inspection.inspectionHelpers import *
from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce 
from src.treeInformation import treeInfo
from src.treeRun import RunStuff

#%% load all data table

uniqseqs = load_uniqseqs()
totdist = load_totdist()
SDRs = load_SDRs()
SDVs = load_SDVs()
data_all = totdist.merge(SDRs, on = 'gene')
data_all = data_all.merge(SDVs, on = ['gene', 'level'])
data_all = data_all.merge(uniqseqs, on = ['gene'])
data_all.drop_duplicates(inplace=True)

#%% 
"""
    RBL1and RBL2 are key proteins in the DREAM complex. 
    High degree of conservacy
"""

RBL1 = get_gene_entry("RBL1", data_all)
RBL2 = get_gene_entry("RBL2", data_all)

seqdata_RBL1 = find_uniqseq_groups('ENSG00000080839___RBL1')
seqdata_RBL2 = find_uniqseq_groups('ENSG00000103479___RBL2')

ggplot(seqdata_RBL1) + geom_bar(aes(x='seq', fill='super')) + labs(title = "RBL1, SUPER") + theme(axis_text_x=element_text(rotation=90, hjust=1))
ggplot(seqdata_RBL2) + geom_bar(aes(x='seq', fill='super')) + labs(title = "RBL2, SUPER") + theme(axis_text_x=element_text(rotation=90, hjust=1))

gene_selection.iloc[:0]['ENSG00000080839___RBL1']
'ENSG00000080839___RBL1' in gene_selection



"""
    ENSG00000164916___FOXK1: 
        - 837 uniqseqm SDR = 7972 (sub) 
        - Ubiquitously expressed in brain, testis and 25 other tissues
        - Central gene
"""


"""
    ENSG00000274736___CCL23
    This gene is selected to gene list because of sub SDR. Have a high super SDR.
    Interesting to see distribution of populations. 
    
"""
CCL23 = get_gene_entry("CCL23", data_all)
seqdata_CCL23 = find_uniqseq_groups('ENSG00000274736___CCL23')
ggplot(seqdata_CCL23) + geom_bar(aes(x='seq', fill='super')) + labs(title = "CCL23, SUPER") + theme(axis_text_x=element_text(rotation=90, hjust=1))
ggplot(seqdata_CCL23) + geom_bar(aes(x='seq', fill='sub')) + labs(title = "CCL23, SUb") + theme(axis_text_x=element_text(rotation=90, hjust=1))

# Filter for Africa (want a closer look at if a sub - population ospecifically is only in the outgroup clade)

CCL23_AFR = seqdata_CCL23.loc[seqdata_CCL23['super'] == 'AFR']
ggplot(CCL23_AFR) + geom_bar(aes(x='seq', fill='sub')) + labs(title = "CCL23, AFR, SUB") + theme(axis_text_x=element_text(rotation=90, hjust=1))

# Result: All subpops included. 
CCL23_outgroup = seqdata_CCL23.loc[seqdata_CCL23['seq'] == 'ENSG00000274736___CCL23_HUMAN__2']
# 146 samples in outgroup, of which 141 from Africa.  
# So for super, this is like 146 / 2548 = 5.5 % of all samples in an outgroup
# And it is 21 % of the African suptyoe that is "permutated" from the rest. 
# But in the sub-population case, where the populaitons are smaller, 

# How many of each sub-population
CCL23_outgroup['sub'].value_counts()

# For LWK, 47 / 103 = 46 % is in the outgorup,
# and for GWD, 24 / 113 = 21 % in the outgroup
# and for ESN, 21 / 100 = 21 % in the outgroup. 

# So in this case many of the 
CCL23_distmat = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000274736___CCL23___CopD.csv')

# If distance is high, will it affect ratio? 


# Simulate with regards to increasing distance betweem samples. 
# 10 in outgorup, 100 in in group, 
# Distance 




RS6 = get_gene_entry("RS6", data_all)





