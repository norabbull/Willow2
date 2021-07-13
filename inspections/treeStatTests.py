# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 15:06:09 2021

@author: norab
"""

from inspections.inspectionHelpers import *
from plotnine import *        # OBS: must change to all needed stuff only
from scipy.stats import kstest, norm
import numpy as np
import matplotlib.pyplot as plt



# =============================================================================
# Load stuff 
# =============================================================================
uniqseqs = load_uniq_seqs()
totdist = load_tot_dist()
SDRs = load_SDRs()
SDRsuper = load_SDRs(level='super')
SDRsub = load_SDRs(level='sub')
SDVs = load_SDVs()
SDRnull = load_SDRnull()
SDRnull47 = load_SDRnull47()


# =============================================================================
# Inspect and check SDR normality 
# =============================================================================

# 1. Visualization: Density, boxplot and QQ

##########################  Boxplot + violin ##################################

shift = 0.1

def alt_sign(x):
    "Alternate +1/-1 if x is even/odd"
    return (-1) ** x

m1 = aes(x=stage('level', after_scale='x+shift*alt_sign(x)'))              # shift outward
m2 = aes(x=stage('level', after_scale='x-shift*alt_sign(x)'), group='gene')  # shift inward

(ggplot(SDRs, aes('level', 'SDR', fill = 'level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.65, show_legend = False)
 + geom_boxplot(width = shift, alpha=0.7, size = 0.65, show_legend = False)
 + scale_fill_manual(values=['dodgerblue', 'darkorange'])
 + labs(title='SDR values for 8782 genes (all)')
)

"""
    Comment: 
        From first inspection, seems to be a lot of outliers in SDR values. 
        Long tails. 
"""

##############################  Density  ######################################


# All SDRs
(ggplot(SDRs, aes(x='SDR', fill='level')) 
 + geom_density(adjust = 1/2, alpha=0.5)
 + scale_fill_manual(values=['dodgerblue', 'darkorange'])
 + labs(title="SDR density, 8782 genes (all)")
 )

"""
    Comment: 
        Looks relativly normally distributed, but a little suspicious about the 
        long tails and high mid top. 
"""
###############################  QQ plot  #####################################

import statsmodels.api as sm
from scipy.stats import norm
import pylab

my_data = norm.rvs(size=1000)
sm.qqplot(my_data, line='45')
pylab.show()


# 2. Test for noramlity with Shapiro Wilk test IF looks suspicious. 


# Data
SDRsuperData = SDRsuper['SDR'].to_numpy()
SDRsubData = SDRsub['SDR'].to_numpy()
my_data = norm.rvs(size=1000)

ks_statistic, p_value = kstest(SDRsuperData, 'norm'); print(ks_statistic, p_value)
ks_statistic, p_value = kstest(SDRsubData, 'norm'); print(ks_statistic, p_value)
ks_statistic, p_value = kstest(my_data, 'norm'); print(ks_statistic, p_value)

"""
    Comment: 
        Interpret as follows: 
            
            - If observed data follows perfect normal distribution, the value
            of the KS statistic will be 0.
            
            P-value is used to decide if difference is large enough to reject 
            null hypothesis: 
            - If the p-values of the KS test is larger than 0.05, we assume 
            a normal distribution.
            - If the p-value of the KS test is smaller than 0.05, we do not
            assume a normal distribution. 
            
        
        Result: 
            - SDRsuper: ks_statistic = 0.74, p-val = 0.0
            - SDRsub:   ks_statistic = 0.74, p-val = 0.0
            
        The Kolmogorov Smirnov test says the SDR values are not normally
        distributed at all. No doubt.
        NB! Check if this test is appropriate more throughly. 
        
"""



# =============================================================================
# Test distributions against eachother - z test ot t test
# =============================================================================
