#!/usr/bin/env python3

import pandas as pd
import numpy as np

from itertools import product
from scipy.stats import fisher_exact

from pathos.multiprocessing import ProcessingPool

from src import logger




def cv(df, axis):
	return  

def compute_fishers_exact(s1, s2): 
	contingency_table = pd.crosstab(s1, s2)
	return fisher_exact(contingency_table)

def pairwise_fishers_exact(df1, df2, n_cpus=24): 


	if df1.shape[1] * df2.shape[1] > 12: 
		indices = pd.DataFrame(product(df1.columns, df2.columns))
		with ProcessingPool(n_cpus) as pool: 
			args = ((df1[feat1], df2[feat2]) for feat1,feat2 in product(df1.columns, df2.columns))
			results = pool.map(lambda args: compute_fishers_exact(*args), args)
			results_df = pd.DataFrame(results, columns=["OR", "pval"])

		results_df = pd.concat([results_df, indices], axis=1)
		results_df = results_df.pivot(index=0, columns=1).rename_axis(index=df1.columns.name, columns=[None, df2.columns.name])
		return results_df
	else: 
		odds_ratios = pd.DataFrame(index=df1.columns, columns=df2.columns)
		pvals = pd.DataFrame(index=df1.columns, columns=df2.columns)

		for feat1, feat2 in product(df1.columns, df2.columns): 
			odds_ratios.loc[feat1, feat2], pvals.loc[feat1, feat2] = compute_fishers_exact(df1[feat1], df2[feat2])

		return pd.concat({"OR": odds_ratios, "pval": pvals}, axis=1)



#--------------------------------------------------------------------------------------------------#
def grouped_ttest(df, grouping): 
	assert grouping.dtype == bool
	df1 = df.loc[:,grouping]
	df2 = df.loc[:,~grouping]
	return pd.DataFrame(stats.ttest_ind(df1, df2, axis=1), columns=df.index, index=["t", "pval"]).T

