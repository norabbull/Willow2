# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:39:57 2021

@author: norab
"""

import sys
import yaml
from treeRun import treeRun
import timeit

def main(configFilepath=None):
    """
    Function: 
        - call function to run program
        - all information given through 1st argument (path to configuration file)
    """
    print("Job initiated")
    # Read configfile argument
    # configFilepath = sys.argv[1]   # uncomment this to run from terminal
    
    configFilepath = configFilepath.strip()
    
    with open(configFilepath, 'r') as c:
        configFile = yaml.safe_load(c)

    # Run
    treeRun(configFile)
    
    
if __name__ == '__main__':
    
    # Test job: 
    configFilepath = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/test_GDR_10genes/job_input/main_config.yml'
    
    # GDR calculation job:
    #configFilepath = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/GDRcalculation/job_input/main_config.yml'
    
    # nullGDR calculation job:
    #configFilepath = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/GDRcalculation/job_input/main_config.yml'
        
    # non-zero phydists job: 
    #configFilepath = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/GDRcalculation/job_input/main_config.yml'
    
    
    #starttime = timeit.default_timer()
    main(configFilepath)
    #endtime = timeit.default_timer()

        
    
    
    #%% Create non-zero distance values file 

"""
    Non-zero phydist values    
    Create file containing total amount of non-zero values in each matrix

"""


# folder with all phydist files 
phydist_files = 'C:/Users/norab/Master/thesis_data/test_data/phydists_10samples'

# calculate percentage of non-zero values in phydist matrices that is non-zero
test = calc_nonZero_phydists(phydist_files)    # function defined above

# save to file 
saveto = 'C:/Users/norab/Master/thesis_data/meta_data/test_nonzero_phydists.csv'
test.to_csv(saveto, index = False, header = True)
