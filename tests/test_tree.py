# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 08:53:43 2021

@author: norab
"""

import pandas as pd
from copy import deepcopy

from src.treeInformation import treeInfo
from src.treeMetrics import treeMetrics

class TestTreeFunctions:   
    
    def __init__(self):
        
        pass
        
    def test_setGeneName(self):
        """
        OK.
        """
        test_gene1 = 'E:/Master/cophenetic_dists/ENSG00000188324___OR6C6___CopD.csv'
        test_gene2 = 'E:/Master/cophenetic_dists/ENSG00000000938___FGR___CopD.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree1 = treeInfo()
        tree1.setup(test_gene1, pop_info)
        tree2 = treeInfo()
        tree2.setup(test_gene2, pop_info)
        
        assert tree1.getGeneName() == "ENSG00000188324___OR6C6"
        assert tree2.getGeneName() == "ENSG00000000938___FGR"
    
    def test_setPopInfo(self):
        """
        OK.
        """
        pop_info_file = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
        tree = treeInfo()
        tree.setPopInfo(pop_info_file)
        pop_info = tree.getPopInfo()
        
        sup = pop_info[0]
        sub = pop_info[1]
        
        assert sup == {'EUR', 'EAS', 'SAS', 'AFR', 'AMR'}
        assert sub == {'MXL', 'PUR', 'TSI', 'PEL', 'PJL', 'MSL', 
                       'CHB', 'ASW', 'ESN', 'STU', 'IBS', 'BEB', 
                       'ACB', 'YRI', 'ITU', 'GWD', 'CHS', 'CDX', 
                       'GBR', 'KHV', 'GIH', 'FIN', 'LWK', 'JPT', 
                       'CLM', 'CEU'}
    
    def test_setSampleInfo(self):
        """
        OK.
        """
        test_gene = 'E:/Master/cophenetic_dists/ENSG00000188324___OR6C6___CopD.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
        tree = treeInfo()
        tree.setup(test_gene, pop_info)
        
        # dist_mat = tree.getDistMat()   -- used for visual inspection
        sample_info = tree.getSampleInfo()
        
        assert isinstance(sample_info, dict)
        assert isinstance(sample_info[2], list)
        assert sample_info[0] == ['EUR___GBR___HG00261', 'EUR', 'GBR']
        assert sample_info[200] == ['EAS___CHS___HG00409', 'EAS', 'CHS']

    def test_getSampleInfo(self):
        """
        OK.
        """
        test_gene_small = 'C:\\Users\\norab\\Master\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'

        tree = treeInfo()
        tree.setup(test_gene_small, pop_info)
        # dist_mat = tree.getDistMat() -- used for visual inspection
        info = tree.getSampleInfo()
        
        assert info[2][1] =='EUR'
        assert info[2][2] =='GBR'
        assert list(info.items())[-1][-1][-1] == 'ESN'

    def test_makePsuedoPops(self):
        """
        OK
        """
        test_gene = 'C:\\Users\\norab\\Master\\Data\\real_tree_data\\dist_mat_subset\\ENSG00000001167___NFYA___CopD.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
        
        num_pops_tests = [6,17]
        for num_pops in num_pops_tests:
          
            tree = treeInfo()
            tree.setup(test_gene, pop_info)
            info = tree.getSampleInfo()
            dist_mat = tree.getDistMat()
            tree.makePsuedoPops(num_pops)
            
            num_samples = len(dist_mat) / num_pops
            
            all_groups = []
            for ind, el in info.items():
                all_groups.append(el[-1])
            
            all_groups = pd.Series(all_groups)
            all_groups.value_counts()
            
            assert len(all_groups) == len(dist_mat)
            assert sum(all_groups.value_counts()) == len(dist_mat)
            
            num_samples = int(round(num_samples))
            for val in all_groups.value_counts():
                assert val in [num_samples, num_samples+1, num_samples-1]
                       
    def test_calcPopDists(self):
        """
        OK.
        
        """
        test_gene_small = 'C:\\Users\\norab\\Master\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeMetrics()
        tree.setup(test_gene_small, pop_info)
        tree.calcPopDists() 
        pop_dists = tree.getPopDists()
        dist_mat = tree.getDistMat()
        
        for i in range(1, len(dist_mat)):
            print(i)
        
        AFR_test = pop_dists['supWith']['AFR']; AFR_test[0] = round(AFR_test[0], 6)
        EUR_test = pop_dists['supBet']['EUR']; EUR_test[0] = round(EUR_test[0], 6)
        LWK_test = pop_dists['subWith']['LWK']; LWK_test[0] = round(LWK_test[0], 6)
        CLM_test = pop_dists['subBet']['CLM']; CLM_test[0] = round(CLM_test[0], 6)
        
        # Manual calculations retrieved from visual inspection of dist_mat
        supWithSum_AFR = round(0.00236 + 0.00354 + 0.0059 + 0.00472 + 0.00118 + 
                               0.00354 + 0.00236 + 0.00236 + 0.00118 + 0.00118, 6)
        supBetSum_EUR = round(4*0.00236 + 4*0 + 4*(0.00118 + 0.00354 + 2*0.00236), 6)
        subWithSum_LWK = round(0.00236 + 0.0059 + 0.00354, 6)
        subBetSum_CLM = round(0.00472 + 5*0.00236 + 2*0.00118 + 0, 6)

        assert AFR_test == [supWithSum_AFR, 10]  # 5 AFR samples, 10 comparisons
        assert EUR_test == [supBetSum_EUR, 24]   # 4 EUR samples, 24 comparisons
        assert LWK_test == [subWithSum_LWK, 3]   # 3 LWK samples, 3 comparisons
        assert CLM_test == [subBetSum_CLM, 9]    # 1 CLM sample, 9 comparisons 

    def test_calcMeanTypeDists(self):
        """
        OK.

        """
        test_gene_small = 'C:\\Users\\norab\\Master\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        test_gene_IDS = 'E:\\Master\\cophenetic_dists\\ENSG00000010404___IDS___CopD.csv'         # Zero distance
        test_gene_DAXX = 'E:\\Master\\cophenetic_dists\\ENSG00000227046___DAXX___CopD.csv'         # Zero distance
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeMetrics()
        tree.setup(test_gene_small, pop_info)
        # dist_mat = tree.getDistMat()     -- used for visual inspection
        # pop_dists = tree.getPopDists()   -- used for visual inspection
        meanTypeDists = tree.getMeanTypeDists()
        
        # Values read from dist_mat
        subWith = round((0+0.0118)/(6+3), 8)
        subBet = round((0.01416+0.01888+0.01888+0.0472+0.0472) / (9+9+9+24+21), 8) 
        supWith = round((0.02832+0)/(10+6), 8)
        supBet = round((0.0472+0.01888+0.0472) / (25+9+24), 8)
        
        assert supWith == meanTypeDists['supWith']
        assert supBet == meanTypeDists['supBet']
        assert subWith == meanTypeDists['subWith']
        assert subBet == meanTypeDists['subBet']
        
        tree2 = treeMetrics()
        tree2.setup(test_gene_IDS, pop_info)
        meanTypeDists = tree2.getMeanPopDists()
        SDR = tree2.calcSDR()
        
        tree3 = treeMetrics()
        tree3.setup(test_gene_DAXX, pop_info)
        meanTypeDists = tree3.getMeanPopDists()
        meanPopDists = tree3.getMeanPopDists()
        SDR = tree3.calcSDR(random=True)
        SDR 
        
    def test_calcMeanPopDists(self):
        """
        OK.
        """
        test_gene_small = 'C:\\Users\\norab\\Master\\Data\\real_tree_data\\dist_mat_test\\MIO_5x5.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'    
        tree = treeMetrics()
        tree.setup(test_gene_small, pop_info)
        # dist_mat = tree.getDistMat()    -- used for visual inspection
        meanPopDists = tree.getMeanPopDists()
        
        # Values read from dist_mat
        supWithAMR = round(0.003 / 1, 6)
        supBetEUR = round((2*0.00108 + 4*0.00054) / 6, 6) 
        subWithGBR = round(0.005 / 1, 6)
        subBetSTU = round((2*0.00108 + 2*0.00054) / 4, 6)
        
        assert supWithAMR == meanPopDists['supWith']['AMR']
        assert supBetEUR == meanPopDists['supBet']['EUR']
        assert subWithGBR == meanPopDists['subWith']['GBR']
        assert subBetSTU == meanPopDists['subBet']['STU']
        
    def test_calcSDR(self):
        """
        OK.

        """
        test_gene_small = 'C:\\Users\\norab\\Master\\Data\\real_tree_data\\dist_mat_test\\MIO_5x5.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeMetrics()
        tree.setup(test_gene_small, pop_info)
        meanTypeDists = tree.getMeanTypeDists()
        tree.calcSDR()
        treeSubSDR = tree.SDRsub
        treeSupSDR = tree.SDRsuper
        
        # Manual calc: 
        # Mean within superpop dist: 0
        # Mean between superpop dist: 
        mWsup = meanTypeDists['supWith']
        mBsup = meanTypeDists['subBet']
        mWsub = meanTypeDists['subWith']
        mBsub = meanTypeDists['subBet']
        
        subSDR = round(mWsub / mBsub, 6)
        supSDR = round(mWsup / mBsup, 6)        
        
        assert subSDR == treeSubSDR
        assert supSDR == treeSupSDR
        
    def test_calcSingleSDRs(self):
        """
        OK.
        """
        test_gene_small = 'C:\\Users\\norab\\Master\\Data\\real_tree_data\\dist_mat_test\\MIO_5x5.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
        
        tree = treeMetrics()
        tree.setup(test_gene_small, pop_info)
        ## dist_mat = tree.getDistMat()        # used for visual inspection
        
        tree.calcSingleSDRs()
        singleSuperSDR = tree.getSingleSuperSDR()
        singleSubSDR = tree.getSingleSubSDR()

        # Manual calc: 
            
        # Super, AMR
        mW = 0.003 / 1
        mB = (6 * 0.00054)/6
        SDR_AMR = round(mW / mB, 3)
        assert SDR_AMR == singleSuperSDR['AMR']
        
        # Sub, GBR
        mW = 0.005 / 1
        mB = ((2*0.00108) + (4*0.00054)) / 6
        SDR_GBR = round(mW / mB, 3)
        assert SDR_GBR == singleSubSDR['GBR']        

    def test_calcSDV(self):       
        """
        OK.
        """
        test_gene = 'E:/Master/cophenetic_dists/ENSG00000175336___APOF___CopD.csv'
        
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeMetrics()
        tree.setup(test_gene, pop_info)
        tree.calcSDV()
        
        # Manual calc:
        # singleSuperSDR = tree.getSingleSuperSDR()   -- used for visual inspection
        # singleSubSDR = tree.getSingleSubSDR()       -- used for visual inspection
            
        # Super
        mean = (0.422 + 0.79 + 0.202 + 0.455 + 0.138) / 5
        n = 5-1     # ddof = 1
        vals = (0.422-mean)**2 + (0.79-mean)**2 + (0.202-mean)**2 + (0.455-mean)**2 + (0.138-mean)**2
        SDVsuper_man = round(vals / n, 4)         
        
        assert tree.getSDVsuper() == SDVsuper_man 
        
        # Sub
        mean = (0.351+0.344+0.181+0.626+0.231+0.501+0.291+0.547+0.5+0.851+0.657+
                0.539+0.285+0.269+0.101+0.512+0.394+0.163+0.618+0.184+0.12+
                0.202+0.399+0.867+0.178+0.995) / 26
        n = 26-1
        vals = ((0.351-mean)**2+(0.344-mean)**2+(0.181-mean)**2+(0.626-mean)**2+
        (0.231-mean)**2+(0.501-mean)**2+(0.291-mean)**2+(0.547-mean)**2+
        (0.5-mean)**2+(0.851-mean)**2+(0.657-mean)**2+(0.539-mean)**2+
        (0.285-mean)**2+(0.269-mean)**2+(0.101-mean)**2+(0.512-mean)**2+
        (0.394-mean)**2+(0.163-mean)**2+(0.618-mean)**2+(0.184-mean)**2+
        (0.12-mean)**2+(0.202-mean)**2+(0.399-mean)**2+(0.867-mean)**2+
        (0.178-mean)**2+(0.995-mean)**2)
        SDVsub_man = round(vals / n, 4)
        
        assert tree.getSDVsub() == SDVsub_man         
        
    def test_shuffleSampleInfo(self):
        test_gene_small = 'C:\\Users\\norab\\Master\\Data\\real_tree_data\\dist_mat_test\\FGR_10x10.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
    
        tree = treeInfo()
        tree.setup(test_gene_small, pop_info)
        # dist_mat = tree.getDistMat() -- used for visual inspection
        info1 = tree.getSampleInfo()
        info1 = info1.deepcopy()
        tree.shuffleSampleInfo()  
        info2 = tree.getSampleInfo()
        
        assert info1 != info2
    
    def test_calcSDRrandom(self):
        test_gene = 'E:/Master/cophenetic_dists/ENSG00000175336___APOF___CopD.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
        
        tree = treeMetrics()
        tree.setup(test_gene, pop_info)
        sample_info = tree.getSampleInfo()        
        sample_info_orig = deepcopy(sample_info)
        tree.calcSDR()
        subSDRorig = tree.getSDRsub()
        superSDRorig = tree.getSDRsuper()

        # New tree with shuffled shit
        tree2 = treeMetrics()
        tree2.setup(test_gene, pop_info)
        tree2.shuffleSampleInfo()
        tree2.calcSDR()
        subSDRtrand = tree2.getSDRsub()
        superSDRrand = tree2.getSDRsuper()
    

    def test_nullSDR_DAXX(self):
        
        gene = 'E:/Master/cophenetic_dists/ENSG00000227046___DAXX___CopD.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
        tree = treeMetrics()
        tree.setup(gene, pop_info)
        
        sample_info = tree.getSampleInfo()
        sample_info_before = deepcopy(sample_info)
        
        tree.shuffleSampleInfo()
        sample_info_shuffled = tree.getSampleInfo()
        
        meanTypeDist = tree.getMeanTypeDists()
        popDists = tree.getPopDists()
        tree.calcSDR()
        tree.getSDRsuper()        
        # Extract number of samples in each generated defined pop
        # There will always be
        d = pd.DataFrame.from_dict(sample_info).transpose()
        d.rename(columns={0:'sample', 1:'super', 2:'sub'}, inplace = True)
        k = d.super.value_counts()
        
# =============================================================================
#         Manual calculation of SDR from DAXX gene
# =============================================================================
        # Distance from 1 AFR sample to all others: 0.00534
        
        # When in AFR:
            # mean within:
        # 1. Total distance from non-zero sample to all other within-samples: 
            0.00534 * 670 = 3.5778
        # 2. Total number of comparisons made within:
            
            #(671!)/(2! * (671-2)!) (calculated online: https://www.calculatorsoup.com/calculators/discretemathematics/combinations.php)
            
            224785 # AFR
            135981 # EUR
            132355  # EAS
            120786 # SAS
            60378  # AMR
            
            # Total number of within pairwise distances: 
            224785 + 135981 + 132355 + 120786 + 60378 = 674285
            
            # Mean between distance: 
            2548-671 = 1877         # Num non-afr samples
            # Num pairwise distances between AFR and other pops
            1877*671 = 1259467
            # Total between distance between the AFR and other pops: 
            1877 * 0.00534 = 10.02318
            
            # Same for other pops: 
            2548-522 = 2026 # EUR
            2026*522 = 1057572
            522 * 0.00534 = 2.78748
            
            2548-515 = 2033 # EAS
            2033*515 = 1046995
            515 * 0.00534 = 2.7501
            
            2548-492 = 2056 # SAS
            2056*429 = 882024
            429 * 0.00534 = 2.29086
            
            2548-348 = 2200 # AMR
            2200*348 = 765600
            348 * 0.00534 = 1.85832
            
            # Total between distance:
            10.02318 + 2.78748 + 2.7501 + 2.29086 + 1.85832 = 19.70994
            # Total between comp = 
            1259467 + 1057572 + 1046995 + 882024 + 765600 = 5011658
        
        # 3. SDR:
            (3.5778 / 674285) / (19.70994 / 5011658)

            # Note:
            #   For between-distance calculation, all apirwise distance values
            #   are saved twice; once for each of the populations. 
            #   The number of pairwise distance comparisons is also saved twice,
            #   which means that the final mean distance will be correct as the 
            #   distance is divided by a double amount of pairwise comparisons. 
            #   it is necessary to save it this way for other types of calculations.
            
# =============================================================================
#             America
# =============================================================================
# 1. Total distance from non-zero sample to all other within-samples: 
            0.00534 * 348 = 1.85832
        # 2. Total number of comparisons made within:
            
            #(671!)/(2! * (671-2)!) (calculated online: https://www.calculatorsoup.com/calculators/discretemathematics/combinations.php)
            
            224785 # AFR
            135981 # EUR
            132355  # EAS
            120786 # SAS
            60378  # AMR
            
            # Total number of within pairwise distances: 
            224785 + 135981 + 132355 + 120786 + 60378 = 674285
            
            # Mean between distance: 
            2548-671 = 1877         # Num non-afr samples
            # Num pairwise distances between AFR and other pops
            1877*671 = 1259467
            # Total between distance between the AFR and other pops: 
            671 * 0.00534 = 3.58314
            
            # Same for other pops: 
            2548-522 = 2026 # EUR
            2026*522 = 1057572
            522 * 0.00534 = 2.78748
            
            2548-515 = 2033 # EAS
            2033*515 = 1046995
            515 * 0.00534 = 2.7501
            
            2548-492 = 2056 # SAS
            2056*429 = 882024
            429 * 0.00534 = 2.29086
            
            2548-348 = 2200 # AMR
            2200*348 = 765600
            2200 * 0.00534 = 11.748
            
            # Total between distance:
            3.58314 + 2.78748 + 2.7501 + 2.29086 + 11.748 = 23.15958
            # Total between comp = 
            1259467 + 1057572 + 1046995 + 882024 + 765600 = 5011658
        
        # 3. SDR:
            (1.85832/ 674285) / (23.15958 / 5011658)
        # AFR (3.5778 / 674285) / (19.70994 / 5011658)
        
# =============================================================================
#             Shuffled inspection:
# =============================================================================
            
            # Funny samle: AFR___LWK___NA19437
            sdrs = []
            for _ in range(20):
                
                tree.shuffleSampleInfo()
                sample_info_shuf = tree.getSampleInfo()
                if sample_info_shuf[2547][1] == 'EUR':
                    tree.calcSDR()
                    sdrs.append(tree.getSDRsuper())
                    meanTypeDist_shuf = tree.getMeanTypeDists()
                    meanPopDist_shuf = tree.getMeanPopDists()
                    
                    # print("Pop dists: ", tree.getPopDists())
                    print("Sample: ", sample_info_shuf[2547])
                    # print("MeanPopWith:", meanPopDist_shuf['supWith'])
                    # print("MeanPopBet:", meanPopDist_shuf['supBet'])
                    # print("MeanTypeWith:", meanTypeDist_shuf['supWith'])
                    # print("MeanTypeBet:", meanTypeDist_shuf['supBet'])
                    print("SDR: ", tree.getSDRsuper())
                    
    def test_DAXX_sampleinfo(self):
        
        gene = 'C:/Users/norab/OneDrive/Skrivebord/ENSG00000227046___DAXX___CopD.csv'
        pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'     
        save_info2 = pd.DataFrame(columns=['SDR', 'AFR', 'EUR', 'EAS', 'AMR', 'SAS','sample_in'])
        
        for i in range(100):
            print("Tree nr: ", i)
            
            tree = treeMetrics()
            tree.setup(gene, pop_info)
        
            sample_info = tree.getSampleInfo()
            tree.calcSDR()
            df = pd.DataFrame.from_dict(sample_info).transpose()
            df.rename(columns={0:'sample', 1:'super', 2:'sub'}, inplace = True)
            counts = df.super.value_counts().to_frame()
            counts.rename(columns={'super': f'normal_{i}'}, inplace=True)
            counts.at['SDR', f'normal_{i}'] = tree.getSDRsuper()
            counts.at['sample_in',f'normal_{i}'] = sample_info[2547][1]
            
            save_info2 = save_info2.append(counts.transpose())
            
            for k in range(10):
                print("Round nr: ", k)
                tree.shuffleSampleInfo()
                tree.calcSDR()
                
                sample_info= tree.getSampleInfo()
                df = pd.DataFrame.from_dict(sample_info).transpose()
                df.rename(columns={0:'sample', 1:'super', 2:'sub'}, inplace = True)
                counts = df.super.value_counts().to_frame()
                counts.rename(columns={'super': f'shuf_{i}_{k}'}, inplace=True)
                counts.at['SDR', f'shuf_{i}_{k}'] = tree.getSDRsuper()
                counts.at['sample_in', f'shuf_{i}_{k}'] = sample_info[2547][1]
                save_info2 = save_info2.append(counts.transpose())
            
        save_info2.to_csv('C:/Users/norab/Master/save_info2_SDR.csv')
        
        
        
        def test_setUnigseqSampleMap(self):
            """
            
            
            """
            # Initate tree
            gene = 'C:/Users/norab/Master/Data/real_tree_data/dist_mat_subset/ENSG00000001167___NFYA___CopD.csv'
            pop_info = 'C:/Users/norab/Master/Data/real_tree_data/phydist_population_classes.tsv'
            map_NFYA = 'E:/Master/Data/other/uniqseq_maps/maps/ENSG00000001167___NFYA_HUMAN__uniq_samplemap.tsv'
            tree = treeInfo()
            tree.setup(gene, pop_info)
            tree.setUniqseqMap(map_NFYA)
            uniqseq_map = tree.getUniqseqMap()     
            uniqseq_count = tree.getUniqseqCount()
        
        
            
            
# if __name__ == "__main__":
    
#     test = TestTreeFunctions()
#     test.test_setGeneName()
#     test.test_setPopInfo()
#     test.test_setSampleInfo()
#     test.test_getSampleInfo()
#     test.test_makePsuedoPops()
#     test.test_calcPopDists()
#     test.test_calcMeanTypeDists()
#     test.test_calcMeanPopDists()
#     test.test_calcSDR()
#     test.test_calcSDV()
#     test.test_calcSingleSDRs()
    