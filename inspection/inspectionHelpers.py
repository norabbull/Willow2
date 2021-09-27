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

def load_uniqseqs(file_path = 'C:/Users/norab/Master/Data/meta_data/9381_uniqseqs.txt'):
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

def load_totdist(path= "C:/Users/norab/Master/Data/meta_data/totalDistancesRefined.txt"):
    """
    Input: None
    Funtion: Load data. Path is specified in function for functoinality.
    Retun: Pd dataframe with values and gene names as indices.
    """
    
    totdist = pd.read_csv(path, "r", delimiter = ",", index_col = 0)
    totdist.rename(index=lambda s: re.sub('_HUMAN__full.*', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    totdist.rename(index=lambda s: re.sub('^.*trees/', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    totdist.reset_index(inplace=True)
    totdist.rename(columns={"index":"gene"}, inplace = True)
    
    return totdist

def load_totdist_pops(path = 'C:/Users/norab/Master/data/other_measures/totdist_pops_all.csv', get_info = False):
    desc = """Informtion on total distance within each defined population
                for all genes.  """
    if get_info: 
        print(desc)
    
    return pd.read_csv(path)
    
def load_SDRs(file_path = 'C:/Users/norab/Master/Data/SDR/SDR_all.csv', level='All'):
    """ 
    Input: path
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """
    SDRs = pd.read_csv(file_path, header=0, index_col=False)
    
    if level == 'super':
        return SDRs[SDRs['level'] == 'super']
    elif level == 'sub':
        return SDRs[SDRs['level'] == 'sub']
    else:
        return SDRs 


def load_SDVs(file_path = 'C:/Users/norab/Master/Data/SDV/SDV_all.csv', get_info = False):
    
    if get_info: 
        print("""pd DataFrame containing SDV values for super and sub popultions 
              for each tree/gene""")
        
    SDVs = pd.read_csv(file_path, header=0, index_col=False)
    
    return SDVs.drop_duplicates(subset=['gene', 'SDV'])

    
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
    
    
def load_SDRnull(file_path = 'C:/Users/norab/Master/Data/SDRnull/all/SDRnull_allGenes_22.07.21.csv',
                 level = "All",
                 gene_set = "All"):
    """ 
    Input: Filepath
    Return: pd DataFrame containing SDR values for super and sub popultions for each tree/gene
    
    """ 
    
    if gene_set == "47genes":
        file_path = "C:/Users/norab/Master/Data/SDRnull/all/SDRnullAll_47genes.csv"
    elif gene_set == "93genes":
        file_path = "C:/Users/norab/Master/Data/SDRnull/all/SDRnullAll_93genes.csv"
        
    SDRs = pd.read_csv(file_path, header=0, index_col=False)
    
    if 'super' in level:
        return SDRs[SDRs['level'] == 'nullSuper']
    elif 'sub' in level:
        return SDRs[SDRs['level'] == 'nullSub']
    else:
        return SDRs 
        
    


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

def load_uniqseq_map(folder_path = 'E:/Master/Data/other/uniqseq_maps/maps/', genes = None):
    """
    returns df with number of each unique sequence for each gene. 
    list of genes should be passed, otherwise all genes in specified 
    folder are included.  

    """
    files = make_filelist(folder_path)
    
    if genes:
        
        files = [f for f in files for g in genes if g in f]
         
    all_genes = pd.DataFrame(columns = ['gene', 'uniqseq_count'])
    for f in files:
        try: 
            file = pd.read_csv(f, header = None, sep = '\s', index_col = 0, engine='python')
            file.columns = ['uniqseq_count']
            file['uniqseq_count'] = file['uniqseq_count'].apply(lambda x: len(x.split(',')))
            file.sort_values(by='uniqseq_count', inplace = True, ascending = False)
            gene_name = re.sub('_HUMAN.*$','', file.index[0])
            file.insert(0, 'gene', gene_name)
            all_genes = all_genes.append(file, ignore_index=True)
        except: 
            print("Something weird with: ", f)
            
    return all_genes
    

def get_gene_entry(gene, df):

    return df[df['gene'].str.contains(gene)]
    
if __name__ == '__main__':
    genes = ['ENSG00000166347___CYB5', 'ENSG00000185946___RNPC3', 'ENSG00000160049___DFFA', 'ENSG00000143278___F13B', 'ENSG00000185101___ANO9']
    test = load_uniqseq_map('E:/Master/Data/other/uniqseq_maps/maps/')
    

#%%
# Test shit

# test = make_filelist('E:/Master/Data/other/uniqseq_maps/maps/')
# test = load_uniqseq_map('E:/Master/Data/other/uniqseq_maps/maps/')



# f_working = pd.read_csv('E:/Master/Data/other/uniqseq_maps/maps/ENSG00000136518___ACL6A_HUMAN__uniq_samplemap.tsv', sep = '\s', header = None, index_col = 0)
# f_noWork = pd.read_csv('E:/Master/Data/other/uniqseq_maps/maps/ENSG00000205922___ONEC3_HUMAN__uniq_samplemap.tsv', header = None,sep = '\s', index_col=0)

# f = pd.read_csv('E:/Master/Data/other/uniqseq_maps/maps/ENSG00000136518___ACL6A_HUMAN__uniq_samplemap.tsv', sep = '\s', index_col = 0, header= None)
# f.columns = ['uniqseq_count']
# f['uniqseq_count']= f['uniqseq_count'].apply(lambda x: len(x.split(',')))
# gene_name = re.sub('_HUMAN.*$','', f.index[0])
# f.insert(0, 'gene', gene_name)

# all_genes = pd.DataFrame(columns = ['gene', 'uniqseq_count'])
# all_genes = all_genes.append(f, ignore_index=True)
    
    