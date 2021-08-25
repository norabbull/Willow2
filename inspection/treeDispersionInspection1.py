# -*- coding: utf-8 -*-
"""
Created on Wed May 12 10:05:17 2021
@author: norab
"""

from inspections.inspectionHelpers import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import re
import seaborn as sns
from plotnine import ggplot, aes, geom_line, geom_violin
from plotnine import *


"""
Inspections: 
    - Single gene distance matrices
    - totdist, uniqseq, SDR, SDV epxloration and play

"""

#%% Load single gene matrix

GNAT2_cd= load_cd_mat('E:/Master/cophenetic_dists/ENSG00000134183___GNAT2___CopD.csv')
COPE_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000105669___COPE___CopD.csv')
HSP7C_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000109971___HSP7C___CopD.csv')
S29A4_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000164638___S29A4___CopD.csv')
VAV2_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000160293___VAV2___CopD.csv')
FGR_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv')
CCNA1_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000133101___CCNA1___CopD.csv')
LAS1L_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000001497___LAS1L___CopD.csv')
TRY3_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000010438___TRY3___CopD.csv')
OR6C6_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000188324___OR6C6___CopD.csv')
MPP5_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000072415___MPP5___CopD.csv')
LY96_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000154589___LY96___CopD.csv')
HNRPR_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000282958___HNRPR___CopD.csv')
DAXX_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000227046___DAXX___CopD.csv')
DRB_cd = load_cd_mat('E:/Master/cophenetic_dists/ENSG00000196126___DRB___CopD.csv')
SMALL_FAKE_cd = load_cd_mat('C:/Users/norab/Master/Data/test_data/test_fakeTrees/test_align_2supx2subX4samples___CopD.csv')


#%% Get all info on one gene

LY96_data = all_data1[all_data1['gene'] == 'ENSG00000154589___LY96']
LY96_data = all_data5[all_data5['gene'] == 'ENSG00000154589___LY96']

# Something strang. Inspect: 

# calculate groupwise SDR for gene
ly96_inf= getSampleInfo(LY96_cd)
grouptypes = getGroupTypes('C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv')
ly96_SDRs = calculateSDR(LY96_cd,ly96_inf,grouptypes)
ly96_SDRs_pops = calculateSDR(LY96_cd,ly96_inf,grouptypes,calc_SDRgroupwise=True)

#%% Load data

# All data
uniqseq = load_uniqseqs()
totdist = load_totdist()
SDRsuper = load_SDRs('C:/Users/norab/Master/Data/SDR/SDRsuper_all.csv')
SDRsub = load_SDRs('C:/Users/norab/Master/Data/SDR/SDRsub_all.csv')
SDVs = load_SDVs('E:\Master\SDV\SDV_values_all.csv')
group_SDRs = load_groupSDRs('E:\Master\SDR\SDR_values_all_groups.csv')
totdists_pops = load_totdist_pops()

#%% 
# Small subset data
group_SDRs = load_groupSDRs('C:/Users/norab/MasterDisaster/Data/SDR/SDR_values_subset.csv')
group_totdists = load_groupTotdists('C:/Users/norab/MasterDisaster/Data/meta_data/group_totdists_subset.csv')

test = group_SDRs[group_SDRs['gene'] == 'ENSG00000154589___LY96']

#%% Filter all_data: 
    
all_data1 =  all_data.drop(all_data[all_data.totdist == 0].index)
all_data2 =  all_data1.drop(all_data1[all_data1.SDR == 1].index)
all_data3 =  all_data2.drop(all_data2[all_data2.SDR > 10].index)
all_data4 = all_data3[~all_data3.isin([np.nan, np.inf, -np.inf]).any(1)]
all_data5 =  all_data4.drop(all_data4[all_data4.totdist <= 0.019].index) # At least 50 samples w/ more than 50 samples aving dists

#%% Inspect SDR, totdist, uniqseq

result = pd.merge(totdist,uniqseq,on="gene")
result = result.merge(SDRsuper, on= 'gene')
result = result.merge(SDRsub, on ='gene')

# Ratio of SDR super and sub (probably uneccesaary)
#result['ratioSupSub'] = SDRratioSupSub['ratioSupSub']

# Filter for genes with less than 30 samples having distances to any other samples
result = result[result['totdist'] > 0.0122549]

# =============================================================================
# La stå
# =============================================================================
# Full
(ggplot(result, aes('uniqseq', 'totdist', fill = 'SDRsub'))
 + geom_point(alpha = 0.7, stroke= .05, size = 3)
 + theme_classic()
 + labs(title='SUB: non-zero sample distances VS nr. unique sequences')
)

(ggplot(result, aes('uniqseq', 'totdist', fill = 'SDR'))
 + geom_point(alpha = 0.7, stroke= .05, size = 3)
 + theme_classic()
 + labs(title='SUPER: non-zero sample distances VS nr. unique sequences')
)

# Xlimited

(ggplot(result, aes('uniqseq', 'totdist', fill = 'SDRsub'))
 + geom_point(alpha = 0.7, stroke= .05, size = 4)
 + theme_classic()
 + labs(title='SUB: non-zero sample distances VS nr. unique sequences')
 + xlim(0,30)
)

(ggplot(result, aes('uniqseq', 'totdist', fill = 'SDR'))
 + geom_point(alpha = 0.7, stroke= .05, size = 4)
 + theme_classic()
 + labs(title='SUPER: non-zero sample distances VS nr. unique sequences')
 + xlim(0,30)
)


(ggplot(result, aes('uniqseq', 'totdist', fill = 'SDR'))
 + geom_point(alpha = 0.7, stroke= .01, size = 5 )
 + theme_classic()
 + labs(title='SUPER: totdist vs uniqseq')
 + xlim(0,30)
)


(ggplot(result, aes('uniqseq', 'totdist', fill = 'ratioSupSub'))
 + geom_point(alpha = 0.7, stroke= .01, size = 5 )
 + theme_classic()
 + labs(title='SUPER: totdist vs uniqseq')
 + xlim(0,100)
)

SDRratioSupSub = pd.DataFrame(result['SDR'] / result['SDRsub'])
SDRratioSubSup = pd.DataFrame(result['SDRsub'] / result['SDR'], column)
SDRratioSupSub.columns = ['ratioSupSub']


#%% Plot

# =============================================================================
# Pairplots
# =============================================================================

# Unique seq vs totdist
plt.figure()
result.plot(x='uniqseq', y='totdist', style='o')
plt.xlabel("unique sequences")
plt.ylabel("total non-zero distances in tree")
plt.show()

# Unique seq vs SDR super
plt.figure()
result.plot(x='SDR_super', y='uniqseq', style='o', color = "red")
plt.xlabel("SDR super")
plt.ylabel("unique sequences")
plt.title("SDR vs number of unique sequences in gene. No samples removed.")
plt.show()

plt.figure()
group_SDR_totdists.plot(x='percentNonZeroForPop', y='SDR', style='o')
plt.xlabel("nonZero, %")
plt.ylabel("SDR")
plt.show()


# All against all pairplots
sns.pairplot(result2)
sns.pairplot(result3)
sns.pairplot(result4)
sns.pairplot(result5)
sns.pairplot(result6)

# All in one plot
result.plot()
result4.plot()

# =============================================================================
# Barchart
# =============================================================================

result['uniqseq'].plot.hist(bins=100, figsize=(8,6))
result['SDR_super'].plot.hist(bins=100, figsize=(8,6))

result2.plot.hist(bins=100, alpha = 0.7)


# =============================================================================
# Scatterplots
# =============================================================================

result.plot.scatter(x='uniqueseq', y='SDR_super')


g = sns.catplot(x="gene", y="SDR", hue="pop", data=result)
g.set_xticklabels(rotation=90)

g = sns.catplot(x="gene", y="SDR", hue="pop", size = "log_uniqseq", jitter = False, data=result)
g.set_xticklabels(rotation=90)

g = sns.relplot(x="gene", y="SDR", hue="pop", size="uniqseq", sizes = (1,300),data=result)
g.set_xticklabels(rotation=90)


# Pairplot of SDR groups

sns.pairplot(result, hue='pop')

# Youre onto something here! Do this with all. look for patterns.
#%% analysis of single trees / genes

# =============================================================================
# Heatmap
# =============================================================================

# VAV2 ; totdist = 0,737422 , uniqseq = 244, SDR super = 0.3339, SDR sub = 0.3061, also higher variance than others

sns.heatmap(VAV2_cd)


#%% Explor group SDRs

# Question to anser:
    # Are there any popultiaons with stronger clustering for multiple genes than others?

#group_SDRs = group_SDRs[~group_SDRs.isin([np.nan, np.inf, -np.inf]).any(1)]
#group_SDRs = group_SDRs.drop(group_SDRs[group_SDRs.SDR ==1].index)

# Extract data from all_data5 (filtered): 
    
group_SDR = all_data5[['SDR', 'gene', 'pop', 'percentNonZeroForPop']]
group_SDRs_pivot = group_SDR.pivot_table(values='SDR', index=group_SDR['pop'], columns='gene', aggfunc='first')
sdr_cov = group_SDRs_pivot.cov()

sns.heatmap(sdr_cov)

g = sns.catplot(x="gene", y="SDR", hue="pop", data=group_SDR)
g.set_xticklabels(rotation=90)

#%% Make super matrix with all values - this exports the matrix "all values" 

# Join group totdists and group sdrs
group_SDRs_totdists1 = pd.merge(group_SDRs, group_totdists,  how='left', left_on=['pop','gene'], right_on = ['pop','gene'])
group_SDRs_totdists2 = pd.merge(group_SDRs_totdists1,totdist,how='outer',left_on='gene',right_index=True, validate='many_to_many')
group_SDRs_totdists3 = pd.merge(group_SDRs_totdists2,SDRs,how='outer',left_on='gene',right_index=True, validate='many_to_many')
group_SDRs_totdists4 = pd.merge(group_SDRs_totdists3,SDVs,how='outer',left_on='gene',right_index=True, validate='many_to_many')

# Export
with open('C:/Users/norab/MasterDisaster/Data/real_tree_data/allTreeInfo.csv', 'a') as f:
    group_SDRs_totdists4.to_csv(f, mode = 'w', header=f.tell()==0)



#%% toGO filter

"""
Filter genes with SDR < 0.7
"""
superSet = load_groupTypeSets('C:/Users/norab/MasterDisaster/Data/real_tree_data/phydist_population_classes.tsv')[0]
SDRfilter07 = all_data5[['SDR_super', 'SDR_sub', 'SDR', 'gene', 'pop']] 
SDRfilter07 = SDRfilter07.drop_duplicates('gene')
#SDRfilter07 = SDRfilter07 = SDRfilter07[SDRfilter07['pop'].isin(superSet)]
SDRfilter07 = SDRfilter07.drop(SDRfilter07[SDRfilter07.SDR_super >= 0.7].index)
SDRfilter07.sort_values(by=['SDR_super'], inplace=True)

sns.scatterplot(x = 'gene', y = 'SDR_super', hue='pop', data = SDRfilter07)

#  Extract gene list

geneUniverse = pd.DataFrame(all_data['gene'])
geneUniverse= pd.DataFrame(geneUniverse.drop_duplicates('gene'))

#  SDRSuper 0.7 treshold 
genes07 = SDRfilter07['gene']

# Export 

genes07.to_csv('C:/Users/norab/MasterDisaster/Data/go_enrichment/superSDRgenes07.csv', index = False)
geneUniverse.to_csv('C:/Users/norab/MasterDisaster/Data/go_enrichment/geneUniverse.csv', index = False)






