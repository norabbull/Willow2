# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 15:41:40 2021

@author: norab
"""

from functools import reduce 


from treeHelpers import *
import pandas as pd
import os
from os.path import isfile, join

import numpy as np

from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine import *  # TO DO: change
import numpy as np
import matplotlib.pyplot as plt
#from functools import reduce 
import math

#%% load all data table

uniqseqs = load_uniqseqs()
totdist = load_totdist()
SDRs = load_SDRs()
data_all = totdist.merge(SDRs, on = 'gene')
data_all = data_all.merge(SDVs, on = ['gene', 'level'])
data_all = data_all.merge(uniqseqs, on = ['gene'])
data_all.drop_duplicates(inplace=True)


#%% Inspect significant genes

data_all_sig_SUB = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRsignificant/SDR1005_significantGenesSUB.csv')
data_all_sig_SUPER = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRsignificant/SDR1005_significantGenesSUPER.csv')


#%% SDR plot


(ggplot(data_all,aes('GDR','uniqseq', fill = 'level')) 
  + geom_point()
  + labs(title = 'SDR vs unique sequences')
  + geom_label(label = "SDR", y = "unique sequences")
  #+ scale_x_continuous(breaks=us_td_sdr[''])
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


seqdata_KRA22 = find_uniqseq_groups('ENSG00000214518___KRA22')
seqdata_NDUS5 = find_uniqseq_groups('ENSG00000168653___NDUS5')
seqdata_HAUS4 = find_uniqseq_groups('ENSG00000092036___HAUS4')

k = seqdata_KRA22[(seqdata_KRA22['seq'] == 'ENSG00000214518___KRA22_HUMAN__2') & (seqdata_KRA22['sub'] == 'FIN')]
l = seqdata_KRA22[(seqdata_KRA22['seq'] == 'ENSG00000214518___KRA22_HUMAN__2')]
#%% ggplots throw

 
ggplot(seqdata_KRA22) + geom_bar(aes(x='seq', fill='super')) + labs(title = "KRA22, SUPER") + theme(axis_text_x=element_text(rotation=90, hjust=1))
ggplot(seqdata_KRA22) + geom_bar(aes(x='seq', fill='sub')) + labs(title = "KRA22, SUB") + theme(axis_text_x=element_text(rotation=90, hjust=1))
ggplot(seqdata_NDUS5) + geom_bar(aes(x='seq', fill='super')) # High SDR sub
ggplot(seqdata_HAUS4) + geom_bar(aes(x='seq', fill='super')) # Many african



#%% 

# Get GDR info of genes
data_all_sig_SUPER[data_all_sig_SUPER['gene'].str.contains('HAUS4')]
data_all_sig_SUB[data_all_sig_SUB['gene'].str.contains('HAUS4')]

data_all_sig_SUPER[data_all_sig_SUPER['gene'].str.contains('NDUS')]
data_all_sig_SUB[data_all_sig_SUB['gene'].str.contains('NDUS')]


#%% Plot KRA22 gene
(ggplot(seqdata_KRA22) 
 + geom_bar(aes(x='seq', fill='super')) 
 + labs(title = "KRA22, SUPER     GDR = 0.394") 
 #+ scale_x_discrete(name="unique samples")
 #+ scale_y_discrete(name="count")
 + xlab("unique 3'UTR variants")
 + theme(figure_size=(13, 13))
 + theme(axis_text_x=element_text(rotation=30, hjust=1),
         plot_title=element_text(color = "black",
                             size = 35,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
          axis_title=element_text(color='black',
                             size=25,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
          axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=25,
     legend_key_height=15,
     legend_key_size=25,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
         )
)

(ggplot(seqdata_KRA22) 
 + geom_bar(aes(x='seq', fill='sub')) 
 + labs(title = "KRA22, SUB     GDR = 0.371") 
 #+ scale_x_discrete(name="unique samples")
 #+ scale_y_discrete(name="count")
 + xlab("unique 3'UTR variants")
 + theme(figure_size=(13, 13))
 + theme(axis_text_x=element_text(rotation=30, hjust=1),
         plot_title=element_text(color = "black",
                             size = 35,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
          axis_title=element_text(color='black',
                             size=25,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
          axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=25,
     legend_key_height=15,
     legend_key_size=25,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
         )
)


#%% Plot NDUS5 gene


(ggplot(seqdata_NDUS5) 
 + geom_bar(aes(x='seq', fill='super')) 
 + labs(title = "NDUS5, SUPER     GDR = 0.804") 
 #+ scale_x_discrete(name="unique samples")
 #+ scale_y_discrete(name="count")
 + xlab("unique 3'UTR variants")
 + theme(figure_size=(13, 13))
 + theme(axis_text_x=element_text(rotation=30, hjust=1),
         plot_title=element_text(color = "black",
                             size = 35,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
          axis_title=element_text(color='black',
                             size=25,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
          axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=25,
     legend_key_height=15,
     legend_key_size=25,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
         )
)


#%% Plot HAUS4 gene


(ggplot(seqdata_HAUS4) 
 + geom_bar(aes(x='seq', fill='super')) 
 + labs(title = "HAUS4, SUPER     GDR = 0.971") 
 #+ scale_x_discrete(name="unique samples")
 #+ scale_y_discrete(name="count")
 + xlab("unique samples")
 + theme(figure_size=(13, 13))
 + theme(axis_text_x=element_text(rotation=30, hjust=1),
         plot_title=element_text(color = "black",
                             size = 35,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
          axis_title=element_text(color='black',
                             size=25,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
          axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=25,
     legend_key_height=15,
     legend_key_size=25,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
         )
)



#%% IF TIME: 

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


