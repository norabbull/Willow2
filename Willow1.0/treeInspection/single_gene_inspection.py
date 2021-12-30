# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 15:41:40 2021

@author: norab



"""

 
import re
import pandas as pd
from plotnine import *  

#%% Load significant genes

data_all_sig_SUB = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRsignificant/SDR1005_significantGenesSUB.csv')
data_all_sig_SUPER = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRsignificant/SDR1005_significantGenesSUPER.csv')


# Get GDRs

data_all_sig_SUPER[data_all_sig_SUPER['gene'].str.contains('KRA22')]
data_all_sig_SUB[data_all_sig_SUB['gene'].str.contains('KRA22')]

data_all_sig_SUPER[data_all_sig_SUPER['gene'].str.contains('HAUS4')]
data_all_sig_SUB[data_all_sig_SUB['gene'].str.contains('HAUS4')]

data_all_sig_SUPER[data_all_sig_SUPER['gene'].str.contains('NDUS')]
data_all_sig_SUB[data_all_sig_SUB['gene'].str.contains('NDUS')]

#%% Data wrange unique sequence information


def find_uniqseq_groups(gene):
    """
    Input: 
        gene: gene to inspect
    
    Funtion: 
        load file containing information of sample distributions across 
        the unique sequences detected for each gene. Organize into 
        pandas dataframe. 
    
    Return: 
        pandas dataframe of loaded info. 
    """
    
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

# Find sample distribution across unique sequences for the genes
# KRA22, NDUS5 and HAUS4, using function created above
seqdata_KRA22 = find_uniqseq_groups('ENSG00000214518___KRA22')
seqdata_NDUS5 = find_uniqseq_groups('ENSG00000168653___NDUS5')
seqdata_HAUS4 = find_uniqseq_groups('ENSG00000092036___HAUS4')

# Inspection of number of Finnish samples of the ENSG00000214518___KRA22_HUMAN__2-
# unique sequence for the KRA22 gene
num_finnish = seqdata_KRA22[(seqdata_KRA22['seq'] == 'ENSG00000214518___KRA22_HUMAN__2') & (seqdata_KRA22['sub'] == 'FIN')]
num_all = seqdata_KRA22[(seqdata_KRA22['seq'] == 'ENSG00000214518___KRA22_HUMAN__2')]



#%% KRA22 

# =============================================================================
# Plot KRA22 gene sample distribution 
# =============================================================================


# figure 6.11

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


# figure 6.12

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


#%% NDUS5

# =============================================================================
# Plot NDUS5 gene sample distribution 
# =============================================================================


# figure 6.9

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


#%% HAUS4 

# =============================================================================
# Plot HAUS4 gene sample distribution 
# =============================================================================

# figure 6.10

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

