# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 15:41:40 2021

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

def find_uniqseq_groups(gene):
    
    file_path = "E:/Master/Data/other/uniqseq_maps/maps/" + gene + "_HUMAN__uniq_samplemap.tsv"
    file = pd.read_csv(file_path, header = None, sep = '\s', index_col = 0, engine='python')
    gene_name = re.sub('_HUMAN.*$','', file.index[0])
    
    seq_data = pd.DataFrame(columns = ['gene', 'seq', 'super', 'sub'])
    
    for row in file.iterrows():
        seq = row[0]
        samples = list(row[1])[0].split(',')
        
        for s in samples: 
            sup = s[0:3]
            sub = s[6:9]
            seq_data.loc[len(seq_data)] = [gene_name, seq, sup, sub]
    
    return seq_data
    
#%% seqdata genes

seqdata_OPLA = find_uniqseq_groups(gene = 'ENSG00000178814___OPLA')
seqdata_SELB = find_uniqseq_groups('ENSG00000284869___SELB')
seqdata_KRA22 = find_uniqseq_groups('ENSG00000214518___KRA22')
seqdata_PTK6 = find_uniqseq_groups('ENSG00000101213___PTK6')
seqdata_TAGL2 = find_uniqseq_groups('ENSG00000158710___TAGL2')
seqdata_CADM3 = find_uniqseq_groups('ENSG00000162706___CADM3')
seqdata_RN135 = find_uniqseq_groups('ENSG00000181481___RN135')
seqdata_NDUS5 = find_uniqseq_groups('ENSG00000168653___NDUS5')


#%% ggplots


ggplot(seqdata_OPLA) + geom_bar(aes(x='seq', fill='sub')) 
ggplot(seqdata_SELB) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_KRA22) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_KRA22) + geom_bar(aes(x='seq', fill='sub')) 
ggplot(seqdata_PTK6) + geom_bar(aes(x='seq', fill='sub')) 
ggplot(seqdata_TAGL2) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_CADM3) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_CADM3) + geom_bar(aes(x='seq', fill='sub'))
ggplot(seqdata_RN135) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_RN135) + geom_bar(aes(x='seq', fill='sub')) 
ggplot(seqdata_NDUS5) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_NDUS5) + geom_bar(aes(x='seq', fill='sub')) 

#%% get gene entries

get_gene_entry("RN135", data_all)
NDUS5 = get_gene_entry("NDUS5", data_all)


#%% 

"""
    Count populations samples in uniqseqs with low sample count.
    Are there any popluations that show high degree of SNVs?

"""

# Load gene list
significant_genes = pd.read_csv('C:/Users/norab/Master/data/SDR/SDR1005_significantGenesFullname.csv',
                                header = None, names = ['gene'])
test = seqdata_OPLA.count(axis = 'columns')
#

seqdata_counts_sup = pd.DataFrame(columns = ['seq', 'super', 'entries'])
seqdata_counts_sup = seqdata_counts_sup.append(test)

test = seqdata_OPLA.groupby(["seq", "super"]).size().reset_index(name="entries") 

gene_list = list(significant_genes['gene'])
seqdata_counts_sup = pd.DataFrame(columns = ['seq', 'super', 'entries'])
for gene in gene_list: 
    seqdata = find_uniqseq_groups(gene)
    counts = seqdata.groupby(["seq", "super"]).size().reset_index(name="entries")
    seqdata_counts_sup = seqdata_counts_sup.append(counts)


# Save 
seqdata_counts_sup.to_csv('C:/Users/norab/Master/data/SDR/seqdata_counts_sup.csv', index = False)












