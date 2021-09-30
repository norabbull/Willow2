# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 16:38:44 2021

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
    
      
    def test_runCalcSDR():
        
        input_files = 'C:\\Users\\norab\\Master\\data\\job_simulation\\cd_simulation'
        subtype_info = 'C:\\Users\\norab\\Master\\data\\job_simulation\\job_input\\simulation_subtype_info.tsv'
        SDR = 'C:\\Users\\norab\\Master\\data\\job_simulation\\job_output\\test_run_30.09\\SDR_runDate{0}.csv'.format(datetime.now().strftime("%d.%m.%Y_%H.%M"))
        save_unprocessed = 'C:\\Users\\norab\\Master\\data\\job_simulation\\job_output\\test_run_30.09'
        
        run = RunStuff()
        