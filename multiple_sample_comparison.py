"""
Process RNAseq data for multiple samples
"""

import pandas as pd 
import numpy as np 
from scipy import stats
import scipy.stats
from scipy.stats import ttest_ind
from scipy.stats import ttest_ind_from_stats
import statsmodels
from statsmodels.stats._knockoff import RegressionFDR
from statsmodels.stats import multitest

df = pd.read_csv('/home/ec2-user/environment/RPKM_cqn_thr_annotatedv2_BCC0029.csv')

#isolate each group
data_BXPC3= [df["BXPC3_COL_a"], df["BXPC3_COL_b"], df["1BxPC3_5"], df["2BxPC3_5"], df["3BxPC3_7.5"], df["4BxPC3_7.5"]]
headers = ["BXPC3_1", "BXPC3_2", "BXPC31_3", "BXPC31_4", "BXPC31_5", "BXPC31_6"]
dfBXPC3 = pd.concat(data_BXPC3, axis=1, keys=headers)
dfBXPC3 = dfBXPC3.replace(0,1)

data_PANC1= [df["PANC1_COL_a"], df["PANC1_COL_b"], df["5Panc1_5"], df["6Panc1_5"], df["7Panc1_7.5"], df["8Panc1_7.5"]]
headers = ["PANC1_1", "PANC1_2", "PANC1_3", "PANC1_4", "PANC1_5", "PANC1_6"]
dfPANC1 = pd.concat(data_PANC1, axis=1, keys=headers)
dfPANC1 = dfPANC1.replace(0,1)

#calculate log 2 fold change
df['mean_log2fc_BXPC3'] = dfBXPC3.mean(axis=1)
df['mean_log2fc_PANC1'] = dfPANC1.mean(axis=1)

df['log2fc_BXPC3_PANC1'] = df['mean_log2fc_BXPC3'].subtract(df['mean_log2fc_PANC1'])

tstat, pvalue = scipy.stats.ttest_ind(dfBXPC3, dfPANC1,axis=1)
df['tstat'] = tstat 
df['pvalue'] = pvalue 

#Bonferroni correction for p value
df['pvalue_adjusted'] = df['pvalue'].multiply(len(df))
dfcorrection = df['pvalue_adjusted'] 
dfcorrection.values[dfcorrection > 1] = 1
df['pvalue_adjusted'] = dfcorrection

df.to_csv("processed.csv") 