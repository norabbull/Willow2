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

#%% SDR plot


(ggplot(data_all,aes('SDR','uniqseq', fill = 'level')) 
  + geom_point()
  + labs(title = 'SDR vs unique sequences')
  + geom_label(label = "SDR", y = "unique sequences")
  #+ scale_x_continuous(breaks=us_td_sdr[''])
)

(ggplot(SDRresults3, aes('SDR', 'Group_size', fill = 'Perm'))
 + geom_point(alpha=1, size=2,stroke = 0.1, color = 'indigo')
# + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
# + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
 #+ scale_fill_manual(values=['dodgerblue', 'darkorange'])
 + theme_classic()
 + theme(figure_size=(8, 6))
 + labs(title='SDR simulation')
 + geom_label(x = "SDR", y = "Group size")
)


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
seqdata_AMPH = find_uniqseq_groups('ENSG00000078053___AMPH')
seqdata_NDKB = find_uniqseq_groups('ENSG00000243678___NDKB')



#%% ggplots


ggplot(seqdata_OPLA) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_SELB) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_KRA22) + geom_bar(aes(x='seq', fill='super')) + labs(title = "KRA22, SUPER") + theme(axis_text_x=element_text(rotation=90, hjust=1))
ggplot(seqdata_KRA22) + geom_bar(aes(x='seq', fill='sub')) + labs(title = "KRA22, SUB") + theme(axis_text_x=element_text(rotation=90, hjust=1))
ggplot(seqdata_PTK6) + geom_bar(aes(x='seq', fill='sub')) 
ggplot(seqdata_TAGL2) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_CADM3) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_CADM3) + geom_bar(aes(x='seq', fill='sub'))
ggplot(seqdata_RN135) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_RN135) + geom_bar(aes(x='seq', fill='sub')) 
ggplot(seqdata_NDUS5) + geom_bar(aes(x='seq', fill='super')) 
ggplot(seqdata_NDUS5) + geom_bar(aes(x='seq', fill='sub')) 
ggplot(seqdata_AMPH) + geom_bar(aes(x='seq', fill='super'))
ggplot(seqdata_NDKB) + geom_bar(aes(x='seq', fill='super')) 


#%% get gene entries

get_gene_entry("RN135", data_all)
NDUS5 = get_gene_entry("NDUS5", data_all)

get_gene_entry("RT35", data_all)

# DPOA2 inspection
DPOA2data = get_gene_entry("DPOA2", data_all)
seqdata_DPOA2 = find_uniqseq_groups('ENSG00000014138___DPOA2')
ggplot(seqdata_DPOA2) + geom_bar(aes(x='seq', fill='super'))


# APOE
APOEdata = get_gene_entry("APOE", data_all)
seqdata_APOE = find_uniqseq_groups('ENSG00000130203___APOE')
ggplot(seqdata_APOE) + geom_bar(aes(x='seq', fill='super'))

# GAPD1
GAPD1data = get_gene_entry("GAPD1", data_all)
seqdata_GAPD1 = find_uniqseq_groups('ENSG00000165219___GAPD1')
ggplot(seqdata_GAPD1) + geom_bar(aes(x='seq', fill='super'))

# GAPD1
GAPD1data = get_gene_entry("GAPD1", data_all)
seqdata_GAPD1 = find_uniqseq_groups('ENSG00000165219___GAPD1')
ggplot(seqdata_GAPD1) + geom_bar(aes(x='seq', fill='super'))


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
seqdata_counts_sub = pd.DataFrame(columns = ['seq', 'sub', 'entries'])
seqdata_counts_sub = seqdata_counts_sub.append(test)

test = seqdata_OPLA.groupby(["seq", "super"]).size().reset_index(name="entries") 

gene_list = list(significant_genes['gene'])
seqdata_counts_sub = pd.DataFrame(columns = ['seq', 'sub', 'entries'])
for gene in gene_list: 
    seqdata = find_uniqseq_groups(gene)
    counts = seqdata.groupby(["seq", "sub"]).size().reset_index(name="entries")
    seqdata_counts_sub = seqdata_counts_sub.append(counts)


# Save 
seqdata_counts_sub.to_csv('C:/Users/norab/Master/data/SDR/seqdata_counts_sub.csv', index = False)


#%%




#%% Inspect results from above

# Read

seqdata_count_sup = pd.read_csv('C:/Users/norab/Master/data/SDR/seqdata_counts_sup.csv')
seqdata_count_sub = pd.read_csv('C:/Users/norab/Master/data/SDR/seqdata_counts_sub.csv')

(ggplot(seqdata_count_sup, mapping = aes(x = "super") 
        + geom_bar()))
ggplot(seqdata_count_sup) + geom_point(aes(x='super', y = 'entries')) 


(ggplot(data="dataset", 
       mapping=aes(x="column_in_dataset", fill='another_column_color')) 
 + geom_density(alpha=0.5)
 + labs(title="title"))





