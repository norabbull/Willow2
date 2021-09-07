# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 08:28:05 2021

@author: norab

    - Inspect uniqseq count map
    - 
"""

from inspection.inspectionHelpers import *
from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce 
from src.treeInformation import treeInfo
from src.treeRun import RunStuff


#%% Load stuff

uniqseqs = load_uniqseqs()
# uniqseq_counts = load_uniqseq_map()  # Should specify gene set. Otherwise far too many. 
totdist = load_totdist()
SDRs = load_SDRs()
SDVs = load_SDVs()
SDRnullAll = load_SDRnull()
SDRnull93 = load_SDRnull(gene_set="93genes")
SDRnull47 = load_SDRnull(gene_set="47genes")
# SDRnull47super = SDRnull47[SDRnull47['level'] == 'super']
# SDRnull47sub = SDRnull47[SDRnull47['level'] == 'sub']

SDRsuper = SDRs[SDRs['level'] == 'super']
SDRsub = SDRs[SDRs['level'] == 'sub']


# Merge all SDRnull data
SDRnullTotal = SDRnullAll.append(SDRnull47)
SDRnullTotal = SDRnullTotal.append(SDRnull93)
SDRnullTotal.drop_duplicates(inplace=True)
SDRnullTotalSuper = SDRnullTotal[SDRnullTotal['level'] == 'nullSuper']
SDRnullTotalSub = SDRnullTotal[SDRnullTotal['level'] == 'nullSub']

#%% Test

mat = treeInfo.getDistMat("E:/Master/cophenetic_dists/ENSG00000063177___RL18___CopD.csv")


#%% Process loaded stuff

# Merge dataframes
data_all = totdist.merge(SDRs, on = 'gene')
data_all = data_all.merge(SDVs, on = ['gene', 'level'])
data_all = data_all.merge(uniqseqs, on = ['gene'])
data_all.drop_duplicates(inplace=True)

#%% Threshold
selection1 = data_all[data_all['totdist'] > 0.0157]    # Genes with less than 20 samples with values filtered out
selection1.sort_values(by=['SDR'], inplace=True, ignore_index=True)


"""
    Select top 1600 genes with regards to SDR. 
    Selection is made for sub-and super level simultaneously, since
    there will be a large overlap and the goal is to make a rough selection
    based on computitional capacity. 
    Largest SDR after selection is still ~ 0.8, which is high. 

"""
selection2 = selection1[selection1.index < 1600]
gene_selection = list(set(list(selection2['gene'])))
gene_selection = pd.DataFrame(gene_selection)

# Save to computer and to disk
gene_selection.to_csv('E:/Master/data/SDRnull/other/SDRnull_genes_selection.csv', index = False)
gene_selection.to_csv('C:/Users/norab/Master/data/SDRnull/other/SDRnull_gene_selection.csv', index = False)

# Filter out genes with already enough values: 
read_genes = pd.read_csv('E:/Master/external_runs/data_software_SDRnull_allGenes_x1000/data/job_input/skip_genes.csv')
skip_genes = list(read_genes['gene'])
genes_to_run_list = [g for g in genes if not g in skip_genes]
genes_to_run = pd.DataFrame(genes_to_run)
genes_to_run.to_csv('E:/Master/data/SDRnull/other/SDRnull_genes_toRun.csv', index = False)


#%% Treshold 2 (Outdated)

"""
    Create genesets based on SDRnull quantiles.
    Find threshold before any filtering --> No reason to not include genes in this process. 
    1. Find 99% quantile SDRs for super and sub
    3. Create genesets
    4. Explore genesets in relation to uniqseq, totdist etc --> Filter
    5. GO!
        
"""


# Quantiles 
SDRsuper_q99 = SDRnullTotalSuper.quantile(.01)  # SDR = 0.910205
SDRsub_q99 = SDRnullTotalSub.quantile(.01)      # SDR = 0.935637

# Create genesets
genes_q99_super = SDRsuper[SDRsuper['SDR'] <= float(SDRsuper_q99)]
genes_q99_sub = SDRsub[SDRsub['SDR'] <= float(SDRsub_q99)]

#%% Filter genes 

# =============================================================================
# Filter for totdist < 30 sample (1.7739 % of all samples have values)
# =============================================================================
"""
    Filter out genes with less samples than 1/3 of the smallest population
    sample group. 
    
    Super: Ad Mixed Americas  - 348 samples
    Num samples: 348 / 3 = 116
    Percentage: 11600 / 2548 = 4.553
    
    Sub: Americans of African Ancestry in SW USA (ASW) - 61 samples
    61 / 2 = 30
    Percentage = 2000 / 2548 = 1.2205

"""


# totdist_SDR = SDRs.merge(totdist, on = 'gene')
# totdist_SDR_filtered = totdist_SDR[totdist_SDR['totdist'] >= 0.012205]
# totdist_SDR_filtered_2 = totdist_SDR[totdist_SDR['totdist'] <= 0.012205]

# Visualize
# SDR + SDRnull
(ggplot(totdist_SDR_filtered, aes(x='SDR', fill='level')) 
     + geom_density(adjust = 1/2, alpha=0.5) 
     + scale_fill_manual(values=['indigo', 'darkorange'])
     + labs(title='SDRs, filtered'))

(ggplot(totdist_SDR, aes(x='SDR', fill='level')) 
     + geom_density(adjust = 1/2, alpha=0.5) 
     + scale_fill_manual(values=['indigo', 'darkorange'])
     + labs(title='SDRs'))

(ggplot(totdist_SDR_filtered) + aes(x="SDR", fill = "level")
      + stat_bin(bins = 1000, alpha = 0.7)
      + scale_fill_manual(values=['indigo', 'darkorange'])
      + labs(title='SDRs, filtered'))

(ggplot(totdist_SDR) + aes(x="SDR", fill = "level")
      + stat_bin(bins = 1000, alpha = 0.7)
      + scale_fill_manual(values=['indigo', 'darkorange'])
      + labs(title='SDRs'))

(ggplot(totdist_SDR_filtered_2) + aes(x="SDR", fill = "level")
      + geom_density(adjust = 1/2, alpha=0.5) 
      + scale_fill_manual(values=['indigo', 'darkorange'])
      + labs(title='SDRs, < 30 sample values'))


# Filter quantile genesets

totdist_super_q99 = genes_q99_super.merge(totdist, on = 'gene', how = 'inner')
totdist_sub_q99 = genes_q99_sub.merge(totdist, on = 'gene', how = 'inner')
super_q99_totdistFiltered = totdist_super_q99[totdist_super_q99['totdist'] >= 0.04553]
sub_q99_totdistFiltered = totdist_sub_q99[totdist_sub_q99['totdist'] >= 0.012205]



#%% Filter2: More than many uniqs

"""
    Investigate: If samples with many uniqseqs also have low SDR. Correlation?
    
"""

# Correlation between uniqseq and SDR

(ggplot(data_all,aes('uniqseq','SDR',fill='level'))
  + geom_point(alpha = 0.6, stroke = 0.2, size = 3)
  + scale_fill_manual(values=['indigo', 'darkorange'])
  + labs(title='SDRs vs uniqseq')
  + theme_light()
)



#%% 
# Inspo
# data = SDRs.merge(uniqseq_counts_subset, on = 'gene')
# data_filter1 = test[test['uniqseq_count'] >= 30]        # Filter out uniqseqs with less than 20 samples having the unique sequence
# data_filter2 = data_filter1.groupby('gene').filter(lambda x: len(x) > 1)# Filter out single standing uniqseq groups
# data_filter3 = data_filter2.groupby('gene').filter(lambda x: len(x) < 20)# Filter out single standing uniqseq groups
# data_filter2
# data2 = SDRs.merge(data_filter2, on='gene')
# data3 = SDRs.merge(data_filter3, on='gene')
# =============================================================================
# Visualize
# =============================================================================


(ggplot(data, aes('gene', 'uniqseq_count', fill = 'SDR') )
 #+ geom_jitter(stroke = 0.2, size = 3)
 + geom_point(stroke = 0.2, size = 3)
 + theme(axis_text_x=element_text(rotation=90, hjust=1))
 )



(ggplot(data, aes('gene', 'SDR', fill = 'uniqseq_count') )
 + geom_jitter(stroke = 0.2, size = 2)
 #+ geom_point(stroke = 0.2, size = 3)
 + theme(axis_text_x=element_text(rotation=90, hjust=1))
 )


(ggplot(data2, aes('gene', 'uniqseq_count', fill = 'SDR') )
 + geom_jitter(stroke = 0.2, size = 2)
 #+ geom_point(stroke = 0.2, size = 3)
 + theme(axis_text_x=element_text(rotation=90, hjust=1))
 )

8


