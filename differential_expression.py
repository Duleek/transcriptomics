"""
Turn text files to dataframes for calculations
"""
import pandas as pd 
import numpy as np 
from scipy import stats
import scipy.stats
from scipy.stats import ttest_ind
from scipy.stats import ttest_ind_from_stats

dfmutated = pd.read_csv('/home/ec2-user/environment/Mutated.txt', delimiter = "\t")
dfnonmutated = pd.read_csv('/home/ec2-user/environment/Nonmutated.txt', delimiter = "\t")

dfcombined = pd.concat([dfmutated, dfnonmutated], axis=1)

#isolate dataframes for comparison from 1 sample 
#replace all 0 values with 1 to prevent errors from using the logarithm function

#sample/group 1
dfreads1 = dfmutated.read_count
dfreads1 = dfreads1.replace(0,1)
dfrpm1 = dfmutated.reads_per_million_miRNA_mapped
dfrpm1 = dfrpm1.replace(0,1)

#sample/group 2
dfreads2 = dfnonmutated.read_count
dfreads2 = dfreads2.replace(0,1)
dfrpm2 = dfnonmutated.reads_per_million_miRNA_mapped
dfrpm2 = dfrpm2.replace(0,1)

#compare group 1 (treatment) to group 2 (control)
#create new columns for the difference in expression between samples
dfcombined["read_count_difference"] = dfreads1.subtract(dfreads2)
dfcombined["reads_per_million_difference"] = dfrpm1.subtract(dfrpm2)

#create new columns for log 2-fold change between sample/group 1 and 2
dfcombined['log2fc_reads'] = np.log2(dfreads1) - np.log2(dfreads2)
dfcombined['log2fc_rpm'] = np.log2(dfrpm1) - np.log2(dfrpm2)


#T-test to calulate p-values

#set sample sizes for groups 1 and 2
n1 = 3
n2 = 3


df_mreads = pd.merge(dfreads1, dfreads2, left_index=True, right_index=True)
df_mrpm = pd.merge(dfrpm1, dfrpm2, left_index=True, right_index=True)

#get means and variances on multisample data
#NOTE: There should be a minimum of 3 samples in each group for a t-test. 
#NOTE: In the example comparison with 1 sample in each group, variance cannot be calulated so the variance was arbitrarily set and n was set to 3 

for ind in df_mreads.index: 
    #print(df_m['read_count_x'][ind], df_m['read_count_y'][ind]) 
    meanreads1 = df_mreads['read_count_x'][ind]
    meanreads2 = df_mreads['read_count_y'][ind]
    sdreads1 = np.sqrt(df_mreads['read_count_x'])
    sdreads2 = np.sqrt(df_mreads['read_count_y'])
    tstat_reads, pvalue_reads = ttest_ind_from_stats(meanreads1, sdreads1, n1, meanreads2, sdreads2, n2)
    
for ind in df_mrpm.index: 
    meanrpm1 = df_mrpm['reads_per_million_miRNA_mapped_x'][ind]
    meanrpm2 = df_mrpm['reads_per_million_miRNA_mapped_y'][ind]
    sdrpm1 = np.sqrt(df_mrpm['reads_per_million_miRNA_mapped_x'])
    sdrpm2 = np.sqrt(df_mrpm['reads_per_million_miRNA_mapped_y'])
    tstat_rpm, pvalue_rpm = ttest_ind_from_stats(meanrpm1, sdrpm1, n1, meanrpm2, sdrpm2, n2)
    
#add t stat and p values to main datafram
dfcombined['tstat_reads'] = tstat_reads 
dfcombined['pvalue_reads'] = pvalue_reads 
dfcombined['tstat_rpm'] = tstat_rpm 
dfcombined['pvalue_rpm'] = pvalue_rpm 

#add adjusted p value column
length = len(dfcombined)
dfcombined['padjusted_reads'] = dfcombined['pvalue_reads'].multiply(length)
dfcombined['padjusted_rpm'] = dfcombined['pvalue_rpm'].multiply(length)

#output file for further analysis/graphing
dfcombined.to_csv("output.csv") 
