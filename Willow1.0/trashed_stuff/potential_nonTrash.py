# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 10:05:01 2021

@author: norab
"""

# IN treeMetrics: 
    
# def calcNonZerosForGroups(self, groupType='all', percent = True):
    #     """
    #     Input: 
    #         dist_mat : DataFrame with distances for a tree
    #         sample_info: a dict with numeric index values as keys and sample
    #         infor as value, stored on format [gene_name, supergroup, subgroup]
    #         group_info: Defines groups in popType. 
    #     Function:
    #         COunt number of 
        
    #     Returns:
    #         Amount of non-zero distances for each subgroup
            
    #     TO DO: 
    #         create group_info and groups with functions within this function.
    #         Have to fix imports and project structure first. 
    #     """        
    #     # find total sample cells in matrix
    #     sampleCount = pd.DataFrame(self.countGroupSamples(self.sample_info, groupType=groupType))
    #     sampleCount['totalGroupComp'] = (sampleCount['sampleCount'] * self.dist_mat.shape[0]) - sampleCount['sampleCount']    
        
    #     # Make dict counter for groups --- TEST!!!
    #     if groupType=='super':
    #         sup_sums = {group : 0 for group in self.group_info[0]}   # Placeholder for sum values
    #         sup_sample_info = dict([val[0], val[1]] for val in self.sample_info.values())
            
    #     elif groupType=='sub':
    #         sub_sums = {group : 0 for group in self.group_info[1]}   # Placeholder for sum values
    #         sub_sample_info = dict([val[0], val[2]] for val in self.sample_info.values())
            
    #     elif groupType=='all':
    #         sup_sums = {group : 0 for group in self.group_info[0]}   # Placeholder for sum values
    #         sup_sample_info = dict([val[0], val[1]] for val in self.sample_info.values())
            
    #         sub_sums = {group : 0 for group in self.group_info[1]}   # Placeholder for sum values
    #         sub_sample_info = dict([val[0], val[2]] for val in self.sample_info.values())
            
    #     #     sums1 = {group : 0 for group in group_info[0]}   # Placeholder for sum values
    #     #     sums2 = {group : 0 for group in group_info[1]}
    #     #     sums = dict(sums1, **sums2)
            
    #     #     sample_info1 = dict([val[0], val[1]] for val in sample_info.values())
    #     #     sample_info2 = dict([val[0], val[2]] for val in sample_info.values())
            
    #     #     sample_info = {**sample_info1, **sample_info2}
            
    #     #     return sample_info1
    #     else: 
    #         raise ValueError("Invalid calssification type. ")
            
    #     # map nonZero and sample_info together into matrix
        
    #     nonZeros = self.calcNonZerosForSamples(self.dist_mat, percent = True)  
    #     if groupType == 'all':
    #         nonZeros['supGroup'] = nonZeros.index.map(sup_sample_info)
    #         nonZeros['subGroup'] = nonZeros.index.map(sub_sample_info)
    #         for key, val in nonZeros.iterrows():
    #             sup_sums[val[2]] += val[0]
    #             sub_sums[val[3]] += val[0]
        
    #         sums = {**sup_sums,**sub_sums}
            
    #     elif groupType == 'super':
    #         raise ValueError("This option is not yet available")
    #     else: 
    #         raise ValueError("This option is not yet available")
            
    #     # Divide non-zero values for group on total group compariosns    
    #     sampleCount['groupSums'] = sampleCount.index.map(sums)
    #     sampleCount['percentNonZeroForGroup'] = sampleCount.groupSums / sampleCount.totalGroupComp
    #     sampleCount['gene'] = self.getGeneName()
        
    #     return sampleCount
    
    
#%% In treeMetrics


def calcNonZerosForSamples(self, percent = True):
    """
    Input:
        dist_mat : DataFrame with distances for a tree

    Returns:
        totDist: dataframe with total amount of pairdistances being non-zero as 
        percent and count per row. 
    """
     
    nonZero_row = pd.DataFrame((self.dist_mat != 0).astype(int).sum(axis=1), columns = ['nonZero_count'])
    nonZero_row['percent'] = round(nonZero_row['nonZero_count'] / (nonZero_row.shape[0] - 1), 4)
    
    return nonZero_row



#%% calcSDR


    def run_calcSDR_old(self, random=False):
        """
        random: If groups should be random or not. Used to calculate null distribution.
        skip_genes: list of gene names to be skiped from calculation. 

        """
        
        try:
            
            input_cd_folder = self.config.get('input_cd_folder').strip()
            group_info = self.input_folder + self.config.get('input_group_info').strip()
            
            output_super = self.output_folder + self.config.get('output_super').strip()
            output_sub = self.output_folder + self.config.get('output_sub').strip()
            output_unprocessed = self.output_folder + self.config.get('output_unprocessed').strip()
            
            file_list = self.make_filelist(input_cd_folder)
            skip_genes = self.input_folder + self.config.get('skip_genes').strip()
            select_genes = self.input_folder + self.config.get('select_genes').strip()
            
            
            if '.csv' in select_genes:
                select_genes= pd.read_csv(select_genes, header=None)
                select_genes.columns = ['gene']
                select_genes= list(select_genes['gene'])
                
                keep_files = []
                for file in file_list:
                    if any(gene in file for gene in select_genes):
                        keep_files.append(file)
                        file_list = keep_files
            
            if '.csv' in skip_genes:
                skip_genes= pd.read_csv(skip_genes, header=None)
                skip_genes.columns = ['gene']
                skip_genes= list(skip_genes['gene'])
                
                remove_files = []
                for file in file_list:
                    if any(gene in file for gene in skip_genes):
                        remove_files.append(file)
        
                file_list = [f for f in file_list if f not in remove_files]
            #logger.info(file_list)
            
            ind = 1
            ind_len = len(file_list)
            
            if random: 
                num_iter = int(self.config.get('num_rand_trees'))
            else:
                num_iter = 1

            for i in range(num_iter):
                for cd_file in file_list:
                    try:         
                        
                        #logger.info("File processed: {0}".format(cd_file))
                        #logger.info("File number: {0} / {1}".format(ind, ind_len))
                        #logger.info("Iter round: {0}".format(i))
                        ind +=1
    
                        tree = treeMetrics()
                        tree.setup(cd_file.strip(), group_info)
                        
                        if random: 
                            tree.shuffleSampleInfo()
    
                        tree.calcSDR()
                        
                        supSDR = tree.getSDRsuper()
                        subSDR = tree.getSDRsub()
                        
                        supSDR = [tree.getGeneName(), supSDR]
                        subSDR = [tree.getGeneName(), subSDR]
                    
                        with open(output_super, 'a', newline='') as f:   # write to file    
                            writer = csv.writer(f)
                            writer.writerow(supSDR)
                            
                        with open(output_sub, 'a', newline='') as f:   # write to file    
                            writer = csv.writer(f)
                            writer.writerow(subSDR)
        
                    except Exception: 
                       
                        #logger.exception("File disrupted:", cd_file)
                        file = str(cd_file)
                        
                        with open(output_unprocessed, 'a') as f: 
                            f.write(file)
                            #open text file
                        pass
                
        except Exception as e: 
            #logger.exception(e)
            print("error with runCalc")
            print(e)
            