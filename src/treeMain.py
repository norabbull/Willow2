# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:39:57 2021

@author: norab
"""

#import treeLogger
import sys
import yaml
from src.treeRun import RunStuff
#import logging

#logger = treeLogger.get_logger(__name__)
 
def main():
    
    # Read input configfile
    configFilepath = sys.argv[1]
    configFilepath = configFilepath.strip()
    
    with open(configFilepath, 'r') as c:
        configFile = yaml.safe_load(c)
    
    # Create logger

    #logger = treeLogger.get_logger(__name__)
    
    # Run
    RunStuff(configFile)
    
    # Log
    #logger.info('Session finished successfully: {0}'.format(datetime.now().strftime("%d.%m.%Y_%H.%M")))


if __name__ == '__main__':
    testpath = 'C:/Users/norab/Master/data/job_simulation/job_input/main_config_simulation.yml'
    sys.argv.append(testpath)
    main()
    

    
    #%% Old stuff
     # Test on lenovo computer
    # configFilepath = 'E:/Master/jobs/job_test/job_input/main_config_calcTest_skip_select.yml'
    # configFilepath = 'C:/Users/norab/Master/data/job_simulation/job_input/main_config_simulation.yml'
    
    # # Config file argument
    # configFilepath = sys.argv[1]
    
    # # Load config arguments
    # configFilepath = configFilepath.strip()
    
    # with open(configFilepath, 'r') as c:
    #     config_file = yaml.safe_load(c)
        
    # log = (config_file.get('output_folder') + config_file.get('output_log')).replace('datetime', datetime.now().strftime("%d.%m.%Y_%H.%M"))
    # file_handler = logging.FileHandler(filename=log)
    # stdout_handler = logging.StreamHandler(sys.stdout)
    # handlers = [file_handler, stdout_handler]
    
    # logging.basicConfig(
    #     level=logging.DEBUG, 
    #     format='[%(asctime)s] {%(name)s:%(lineno)d} %(levelname)s - %(message)s',
    #     handlers=handlers
    # )

    # logger = logging.getLogger(__name__)
    
    # run = RunStuff(config_file)
    # run.main()
