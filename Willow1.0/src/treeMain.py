# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:39:57 2021

@author: norab
"""

import sys
import yaml
from src.treeRun import treeRun
 
def main():
    """
    Function: 
        - call function to run program
        - all information given through 1st argument (path to configuration file)
    """
    
    # Read input configfile
    configFilepath = sys.argv[1]
    configFilepath = configFilepath.strip()
    
    with open(configFilepath, 'r') as c:
        configFile = yaml.safe_load(c)

    # Run
    treeRun(configFile)
    
if __name__ == '__main__':
    testpath = 'C:/Users/norab/Master/data/job_simulation/job_input/main_config_simulation.yml'
    sys.argv.append(testpath)
    main()
