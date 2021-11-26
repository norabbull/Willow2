# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:39:57 2021

@author: norab
"""

import sys
import yaml
from treeRun import treeRun
 
def main(configFilepath):
    """
    Function: 
        - call function to run program
        - all information given through 1st argument (path to configuration file)
    """
    
    # Read configfile argument
    configFilepath = sys.argv[1]
    
    configFilepath = configFilepath.strip()
    
    with open(configFilepath, 'r') as c:
        configFile = yaml.safe_load(c)

    # Run
    treeRun(configFile)
    
    
if __name__ == '__main__':
    
    #configFilepath = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/test_GDR_10genes/job_input/main_config.yml'
    print("Job initiated")
    main()
    
    
    
