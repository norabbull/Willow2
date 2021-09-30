# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 13:13:50 2021
@author: norab
"""


import unittest
from datetime import datetime
from src.treeRun import RunStuff 

#%% 

class TestRunStuff(unittest.TestCase):
   
    def __init__():
        pass
        #RunStuff.__init__(self)
    
        
    def test_make_filelist():
        
        input_files = 'E:\\Master\\test_runs\\cophenetic_dists_selection'
        config = 'E:\\Master\\jobs\\job_test\\configs\\main_config_calcTest.yml'
        job = RunStuff(config)
        
        # Make filelist from directory
        filelist = job.make_filelist(input_files)        
   
        # Make filelist from saved filelist
        input_files = 'E:\\Master\\jobs\\job_calcSDRnull\\job_output\\unprocessed_genes_30.06.2021_13.36.csv'
        filelist = job.make_filelist(input_files)

     
    def test_runCalcSDR():
        
        input_files = 'E:\\Master\\test_runs\\cophenetic_dists_selection'
        pop_info = 'E:\\Master\\otherTreeData\\phydist_population_classes.tsv'
        SDRsuper = 'E:\\Master\\test_runs\\SDRsuper_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
        SDRsub = 'E:\\Master\\test_runs\\SDRsub_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))  # dd/mm/YY H:M
        save_unprocessed = 'E:\\Master\\test_runs\\unprocessed_files.csv'
        
        run = RunStuff()
        run.run_calcSDR(input_files = input_files, 
                             pop_info = pop_info,
                             SDRsuper_output = SDRsuper,
                             SDRsub_output = SDRsub,
                             unprocessed_files = save_unprocessed)

    def test_runCalcSingleSDR():
        
        input_files = 'E:\\Master\\test_runs\\cophenetic_dists_selection'
        pop_info = 'E:\\Master\\otherTreeData\\phydist_population_classes.tsv'
        SSDR_output_dir = 'E:\\Master\\test_runs\\singleSDR_test\\'
        unprocessed_files = 'E:\\Master\\test_runs\\singleSDR_test\\unprocessed_testRun17.06.21_16.00.csv'
        
        run = RunStuff()
        run.run_calcSingleSDRs(
            input_files = input_files, 
            pop_info = pop_info,
            SSDR_output_dir = SSDR_output_dir,
            unprocessed_files = unprocessed_files)
        
    
    def test_runCalcSDRnullDist():
        
        config_file = 'E:/Master/test_runs/nullSDR_test/test_25.06.21/test_makeNullDistSDR.yml'
        run = RunStuff(config_file)
        run.main()
        
        # Test 2
        config_file = 'E:/Master/test_runs/nullSDR_test/test_21.07.21/configs/main_config_test_nullSDR_47genes.yml'
        run = RunStuff(config_file)
        run.main()
        
        # Test 3
        config_file = 'E:/Master/test_runs/nullSDR_test/test_21.07.21/configs/main_config_test_nullSDR_DAXX.yml'
        run = RunStuff(config_file)
        run.main()
        
        
    def test_arbitrary_stuff():
                ## REMOVE
        test_gene = 'C:/Users/norab/Master/data/job_simulation/cd_simulation/sim_mat5.csv'
        pop_info = 'C:/Users/norab/Master/data/job_simulation/job_input/simulation_group_classes.tsv'
        config = 'C:/Users/norab/Master/data/job_simulation/job_input/main_config_simulation.yml'
        tree = RunStuff(config)
        tree.setup(test_gene, pop_info)
        gn = tree.getGeneName()
        cd = tree.getDistMat()
        pop_names = tree.getPopInfo()
        si = tree.getSampleInfo()
        tree.getPopInfo()
        ##
    
        

    
        
#%% Runs

if __name__ == '__main__':
    
    #TestRunStuff.test_runCalcSDR()
    
    TestRunStuff.test_runCalcSingleSDR()