# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 08:28:05 2021

@author: norab

Inspect uniqseq count map
"""

from inspections.inspectionHelpers import *
from plotnine import ggplot, aes, geom_point, geom_line, labs
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce 



# =============================================================================
# Load stuff 
# =============================================================================
uniqseqs = load_uniqseqs()
totdist = load_totdist()
SDRs = load_SDRs()
SDVs = load_SDVs()
SDRnull = load_SDRnull()
uniqseq_counts = load_uniqseq_map()

# =============================================================================
# Merge dfs
# =============================================================================

data = SDRs.merge(uniqseq_counts_subset, on = 'gene')
data_filter1 = test[test['uniqseq_count'] >= 30]        # Filter out uniqseqs with less than 20 samples having the unique sequence
data_filter2 = data_filter1.groupby('gene').filter(lambda x: len(x) > 1)# Filter out single standing uniqseq groups
data_filter3 = data_filter2.groupby('gene').filter(lambda x: len(x) < 20)# Filter out single standing uniqseq groups
data_filter2
data2 = SDRs.merge(data_filter2, on='gene')
data3 = SDRs.merge(data_filter3, on='gene')
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


