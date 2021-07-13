# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 09:08:48 2021

@author: norab
"""


import pandas as pd
import re
import os
from os.path import isfile, join
from copy import deepcopy
from random import shuffle, seed


def make_filelist(input_files):

    if isinstance(input_files, str):     
        files = [join(input_files, f) for f in os.listdir(input_files) 
                     if isfile(join(input_files, f))]
    return files

def load_uniq_seqs(file_path = "C:/Users/norab/Master/Data/meta_data/9381_uniqseqs.txt"):
    """
    Input: None
    Funtion: Load data. Path is specified in function for functoinality.
    Retun: Pd dataframe with values and gene names as indices.
    """
    
    uniqseq = pd.read_csv(file_path, delimiter = ":", header = 0, index_col = None) 
    uniqseq.rename(index=lambda s: re.sub('_HUMAN__uniq.*', '', s), 
                   columns = {uniqseq.columns[0]: "uniqseq"}, inplace = True)
    uniqseq.reset_index(inplace=True)
    uniqseq.rename(columns = {"index": "gene"}, inplace=True)
    return uniqseq

def load_tot_dist(totdist_dir = "C:/Users/norab/Master/Data/meta_data/totalDistancesRefined.txt"):
    """
    Input: None
    Funtion: Load data. Path is specified in function for functoinality.
    Retun: Pd dataframe with values and gene names as indices.
    """
    
    totdist = pd.read_csv(totdist_dir, "r", delimiter = ",", index_col = 0)
    totdist.rename(index=lambda s: re.sub('_HUMAN__full.*', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    totdist.rename(index=lambda s: re.sub('^.*trees/', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    totdist.reset_index(inplace=True)
    totdist.rename(columns={"index":"gene"}, inplace = True)
    
    return totdist


def load_SDRs(file_path = 'C:/Users/norab/Master/Data/SDR/SDR_all.csv'):
    """ 
    Input: path
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """
    
    return pd.read_csv(file_path, header=0, index_col=False)
    
    # SDRs.rename(index=lambda s: re.sub('^.*.*ENS', 'ENS', s), inplace = True)
    # SDRs.columns = ['level', 'gene', 'SDR']

    

def load_SDVs(file_path = 'C:/Users/norab/Master/Data/SDV/SDV_all.csv'):
    """
    Input: None
    Return: pd DataFrame containing SDV values for super and sub popultions for each tree/gene
    """

    return pd.read_csv(file_path, header=0, index_col=False)
    #SDVs.rename(index=lambda s: re.sub('^.*.*ENS', 'ENS', s), inplace = True)
    

def load_singleSDR(gene_name, load_all = False):
    """ 
    Input: None
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """
    if load_all: 
        file_path = deepcopy(gene_name)
        gene_name = re.sub('^.*SSDR_', '', gene_name)
        gene_name = re.sub('.csv', '', gene_name)
    else: 
        file_path = 'C:/Users/norab/Master/Data/singleSDRs/SSDR_' + gene_name + '.csv'

        
    SSDRs = pd.read_csv(file_path, header=0, index_col=0)   
    SSDRs.columns = ['pop', f'SSDR_{gene_name}']
    return SSDRs

def load_allSingleSDRs(file_path = 'C:/Users/norab/Master/Data/singleSDRs/', num_rand_trees = 50):
    
    file_list = make_filelist(file_path)
    
    seed(num_rand_trees)
    shuffle(file_list)
    SSDR_all = load_singleSDR(file_list.pop(), load_all = True)
    
    n = 0
    for f in file_list: 
        SSDR_current = load_singleSDR(f, load_all = True)
        SSDR_all = SSDR_all.merge(SSDR_current, on = "pop")
        if n >= num_rand_trees: 
            break
        n +=1

    return SSDR_all
    
    
def load_SDRnull(file_path = 'C:/Users/norab/Master/Data/SDRnull/all/SDRnull_all_07.07.2021.csv'):
    """ 
    Input: Filepath
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """    
    return pd.read_csv(file_path, header=0, index_col=None)
    

    
def load_SDRnull47(file_path = 'C:/Users/norab/Master/Data/SDRnull/all/SDRnull47_all_12.07.2021.csv'):
    """ 
    Input: Filepath
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    """    
    return pd.read_csv(file_path, header=0, index_col=None)


def load_cd_mat(path):
    """
    Input:
        Name: String -  name of gene to find on the form ENSG00000071205___RHG10
        file_dir: Directory with file  of files 
    Function. 
        Read cophenetic distance matrix of a given gene
    Return: 
        
        Distance matrix
    """
    
    return pd.read_csv(path, index_col = 0)
