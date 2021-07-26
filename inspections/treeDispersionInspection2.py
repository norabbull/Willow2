# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 13:27:05 2021
@author: norab
Module to inspect data. 
Loads data form files. 
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import re
from plotnine import ggplot, aes, geom_line, geom_violin
from plotnine import *


"""
Inspections: 
    
    nullSDR - shuffled data


"""


#%% nullSDR - shuffled data

# =============================================================================
# DAXX sample info
# =============================================================================


DAXX_samples1 = pd.read_csv('C:/Users/norab/Master/save_info_SDR.csv')
DAXX_samples2 = pd.read_csv('C:/Users/norab/Master/save_info2_SDR.csv')

len(DAXX_samples1['SDR'].unique())
len(DAXX_samples2['SDR'].unique())
