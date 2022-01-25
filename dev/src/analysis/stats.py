#!/usr/bin/env python3

from itertools import product
from pathos import multiprocessing

import numpy as np
import pandas as pd
from scipy import stats, linalg

import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)

from src import logger


def cv(df, axis):
	return  

def fast_crosstab(s1, s2): 
	df = pd.concat({'s1': s1, 's2': s2}, axis=1)
	counts = (df['s1'] * 2 + df['s2']).value_counts().sort_index()
	return counts.values.reshape(2,2)

def compute_fishers_exact(s1, s2, fast=True):
	if fast:  
		table = fast_crosstab(s1, s2)
		return stats.fisher_exact(table)
	else: 
		try:
			contingency_table = pd.crosstab(s1, s2)
			return stats.fisher_exact(contingency_table)
		except ValueError: 
			return np.nan, np.nan

def _compute_fishers_exacts(df_pairs): 
	results = []
	for s1, s2 in df_pairs: 
		results.append(compute_fishers_exact(s1, s2))
	return results

#----------------------------------------------------------------------------------------------------#
# Functions to perform pairwise feature comparisons among two boolean dfs with same indices (i.e. peaks)
#  - df1: (tf binding x peaks)
#  - df2: (genomic feature x peaks)
#----------------------------------------------------------------------------------------------------#
def pairwise_fishers_exact(df1, df2=None, n_cpus=-1): 

	if n_cpus != -1: 
		indices = pd.DataFrame(product(df1.columns, df2.columns))
		with multiprocessing.ProcessingPool(n_cpus) as pool: 
			args = ((df1[feat1], df2[feat2]) for feat1,feat2 in product(df1.columns, df2.columns))
			results = pool.map(lambda args: compute_fishers_exact(*args), args)
			results_df = pd.DataFrame(results, columns=["OR", "pval"])

		results_df = pd.concat([results_df, indices], axis=1)
		results_df = results_df.pivot(index=0, columns=1).rename_axis(index=df1.columns.name, columns=[None, df2.columns.name])
		# return results_df.astype(float).sort_index(level=1, axis=1)
	else: 
		odds_ratios = pd.DataFrame(index=df1.columns, columns=df2.columns)
		pvals = pd.DataFrame(index=df1.columns, columns=df2.columns)

		for feat1, feat2 in product(df1.columns, df2.columns): 
			odds_ratios.loc[feat1, feat2], pvals.loc[feat1, feat2] = compute_fishers_exact(df1[feat1], df2[feat2])

		results_df = pd.concat({"OR": odds_ratios, "pval": pvals}, axis=1).astype(float)

	results_df = results_df.astype(float).sort_index(axis=1, level=[0,1], ascending=[True, True])
	results_df.OR = np.log2(results_df.OR)
	results_df.pval = -np.log10(results_df.pval)
	results_df = results_df.rename(columns={"OR": "LOD", "pval": "-log10_pval"}, level=0)
	return results_df



def pairwise_feature_counts(df1, df2): 
	"""
	Returns df where index corresponds to df1 variables (columns) and columns correspond
	to df2 variables. 
	"""
	assert df1.shape[0] == df2.shape[0]

	df1_vals = df1.values.astype(int)
	df2_vals = df2.values.astype(int)

	counts = df1_vals.T.dot(df2_vals)
	return pd.DataFrame(counts, index=df1.columns, columns=df2.columns)

def pairwise_contingency_table(df1, df2): 
	return pd.concat({
	    "False": pd.concat({False: pairwise_feature_counts(~df1, ~df2), True: pairwise_feature_counts(~df1, df2)}, axis=1),
	    "True": pd.concat({True: pairwise_feature_counts(df1, df2), False: pairwise_feature_counts(df1, ~df2)}, axis=1), 
	}, axis=0).swaplevel(0,1, axis=1).swaplevel(0,1, axis=0)

# .reindex(columns=[False, True], index=[False, True])
#--------------------------------------------------------------------------------------------------#
def grouped_ttest(df, grouping): 
	assert grouping.dtype == bool
	df1 = df.loc[:,grouping]
	df2 = df.loc[:,~grouping]
	return pd.DataFrame(stats.ttest_ind(df1, df2, axis=1), columns=df.index, index=["t", "pval"]).T

