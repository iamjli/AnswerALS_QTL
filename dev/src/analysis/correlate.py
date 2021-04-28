#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy.stats import t, multitest

from src import logger
from src.analysis import utils


__all__ = [
	"pearson_by_row", "pearson_by_pairs", 
	# "pairwise_cramers_v", "pairwise_theils_u", "pairwise_correlation_ratio",
]

#----------------------------------------------------------------------------------------------------#
# Row-wise correlations (columns are identical)
#----------------------------------------------------------------------------------------------------#
def pearson_by_row(df1, df2, dof=None, fdr_method="bonferroni"): 
	"""Computes pearson statistics for corresponding rows."""

	# Ensure input dataframes have the same indices
	assert df1.columns.equals(df2.columns), "columns are not the same"
	if not df1.index.equals(df2.index): df1, df2 = utils.unify_dfs(df1, df2)

	if dof is None: dof = df1.shape[1] - 2

	vals1, vals2 = df1.values.T, df2.values.T

	vals_norm1 = vals1 - vals1.mean(axis=0)
	vals_norm2 = vals2 - vals2.mean(axis=0)

	covariance = (vals_norm1 * vals_norm2).sum(axis=0)
	variance1 = np.power(vals_norm1, 2).sum(axis=0)
	variance2 = np.power(vals_norm2, 2).sum(axis=0)

	with np.errstate(divide="ignore", invalid="ignore"):
		r = covariance / np.sqrt(variance1 * variance2)
		pvals = get_pvalue_from_corr_coef(r, dof)
		fdr = multitest.multipletests(pvals, method=fdr_method)[1]

	return pd.DataFrame({"pearson_r": r, "pval": pvals, "fdr": fdr}, index=df1.index)


def pearson_by_pairs(df1, df2, pairs_df, dof=None, fdr_method="bonferroni"): 

	assert df1.columns.equals(df2.columns), "columns are not the same"

	df1 = df1.reindex(pairs_df.iloc[:,0]).set_index(pairs_df.index)
	df2 = df2.reindex(pairs_df.iloc[:,1]).set_index(pairs_df.index)

	results = pearson_by_row(df1, df2, dof=dof, fdr_method=fdr_method)
	return pd.concat([pairs_df, results], axis=1)

	if dof is None: dof = df1.shape[1] - 2


def get_pvalue_from_corr_coef(coeffs, dof): 
	tstat = np.abs(coeffs * dof / np.sqrt(1 - np.power(coeffs, 2)))
	return t.sf(tstat, dof) * 2


#----------------------------------------------------------------------------------------------------#
# Column-wise correlations
# 
# Used to compare feature matrices. Description and implementations here:
# 	https://towardsdatascience.com/the-search-for-categorical-correlation-a1cf7f1888c9
#	https://github.com/shakedzy/dython/blob/master/dython/nominal.py
#----------------------------------------------------------------------------------------------------#
# from collections import Counter

# _REPLACE = 'replace'
# _DROP = 'drop'
# _DROP_SAMPLES = 'drop_samples'
# _DROP_FEATURES = 'drop_features'
# _SKIP = 'skip'
# _DEFAULT_REPLACE_VALUE = 0.0


# def _inf_nan_str(x):
#     if np.isnan(x):
#         return 'NaN'
#     elif abs(x) == np.inf:
#         return 'inf'
#     else:
#         return ''

# def conditional_entropy(x, y, nan_strategy=_REPLACE, nan_replace_value=_DEFAULT_REPLACE_VALUE, log_base: float = math.e):
#     """
#     Calculates the conditional entropy of x given y: S(x|y)
#     Wikipedia: https://en.wikipedia.org/wiki/Conditional_entropy
#     Parameters:
#     -----------
#     x : list / NumPy ndarray / Pandas Series
#         A sequence of measurements
#     y : list / NumPy ndarray / Pandas Series
#         A sequence of measurements
#     nan_strategy : string, default = 'replace'
#         How to handle missing values: can be either 'drop' to remove samples
#         with missing values, or 'replace' to replace all missing values with
#         the nan_replace_value. Missing values are None and np.nan.
#     nan_replace_value : any, default = 0.0
#         The value used to replace missing values with. Only applicable when
#         nan_strategy is set to 'replace'.
#     log_base: float, default = e
#         specifying base for calculating entropy. Default is base e.
#     Returns:
#     --------
#     float
#     """
#     if nan_strategy == _REPLACE:
#         x, y = replace_nan_with_value(x, y, nan_replace_value)
#     elif nan_strategy == _DROP:
#         x, y = remove_incomplete_samples(x, y)
#     y_counter = Counter(y)
#     xy_counter = Counter(list(zip(x, y)))
#     total_occurrences = sum(y_counter.values())
#     entropy = 0.0
#     for xy in xy_counter.keys():
#         p_xy = xy_counter[xy] / total_occurrences
#         p_y = y_counter[xy[1]] / total_occurrences
#         entropy += p_xy * math.log(p_y / p_xy, log_base)
#     return entropy

# def cramers_v(x, y, bias_correction=True, nan_strategy=_REPLACE, nan_replace_value=_DEFAULT_REPLACE_VALUE):
#     """
#     Calculates Cramer's V statistic for categorical-categorical association.
#     This is a symmetric coefficient: V(x,y) = V(y,x)
#     Original function taken from: https://stackoverflow.com/a/46498792/5863503
#     Wikipedia: https://en.wikipedia.org/wiki/Cram%C3%A9r%27s_V
#     Parameters:
#     -----------
#     x : list / NumPy ndarray / Pandas Series
#         A sequence of categorical measurements
#     y : list / NumPy ndarray / Pandas Series
#         A sequence of categorical measurements
#     bias_correction : Boolean, default = True
#         Use bias correction from Bergsma and Wicher,
#         Journal of the Korean Statistical Society 42 (2013): 323-328.
#     nan_strategy : string, default = 'replace'
#         How to handle missing values: can be either 'drop' to remove samples
#         with missing values, or 'replace' to replace all missing values with
#         the nan_replace_value. Missing values are None and np.nan.
#     nan_replace_value : any, default = 0.0
#         The value used to replace missing values with. Only applicable when
#         nan_strategy is set to 'replace'.
#     Returns:
#     --------
#     float in the range of [0,1]
#     """
#     if nan_strategy == _REPLACE:
#         x, y = replace_nan_with_value(x, y, nan_replace_value)
#     elif nan_strategy == _DROP:
#         x, y = remove_incomplete_samples(x, y)
#     confusion_matrix = pd.crosstab(x, y)
#     chi2 = ss.chi2_contingency(confusion_matrix)[0]
#     n = confusion_matrix.sum().sum()
#     phi2 = chi2 / n
#     r, k = confusion_matrix.shape
#     if bias_correction:
#         phi2corr = max(0, phi2 - ((k - 1) * (r - 1)) / (n - 1))
#         rcorr = r - ((r - 1) ** 2) / (n - 1)
#         kcorr = k - ((k - 1) ** 2) / (n - 1)
#         if min((kcorr - 1), (rcorr - 1)) == 0:
#             warnings.warn(
#                 "Unable to calculate Cramer's V using bias correction. Consider using bias_correction=False",
#                 RuntimeWarning)
#             return np.nan
#         else:
#             return np.sqrt(phi2corr / min((kcorr - 1), (rcorr - 1)))
#     else:
#         return np.sqrt(phi2 / min(k - 1, r - 1))

# def theils_u(x, y, nan_strategy=_REPLACE, nan_replace_value=_DEFAULT_REPLACE_VALUE):
#     """
#     Calculates Theil's U statistic (Uncertainty coefficient) for categorical-
#     categorical association. This is the uncertainty of x given y: value is
#     on the range of [0,1] - where 0 means y provides no information about
#     x, and 1 means y provides full information about x.
#     This is an asymmetric coefficient: U(x,y) != U(y,x)
#     Wikipedia: https://en.wikipedia.org/wiki/Uncertainty_coefficient
#     Parameters:
#     -----------
#     x : list / NumPy ndarray / Pandas Series
#         A sequence of categorical measurements
#     y : list / NumPy ndarray / Pandas Series
#         A sequence of categorical measurements
#     nan_strategy : string, default = 'replace'
#         How to handle missing values: can be either 'drop' to remove samples
#         with missing values, or 'replace' to replace all missing values with
#         the nan_replace_value. Missing values are None and np.nan.
#     nan_replace_value : any, default = 0.0
#         The value used to replace missing values with. Only applicable when
#         nan_strategy is set to 'replace'.
#     Returns:
#     --------
#     float in the range of [0,1]
#     """
#     if nan_strategy == _REPLACE:
#         x, y = replace_nan_with_value(x, y, nan_replace_value)
#     elif nan_strategy == _DROP:
#         x, y = remove_incomplete_samples(x, y)
#     s_xy = conditional_entropy(x, y)
#     x_counter = Counter(x)
#     total_occurrences = sum(x_counter.values())
#     p_x = list(map(lambda n: n / total_occurrences, x_counter.values()))
#     s_x = ss.entropy(p_x)
#     if s_x == 0:
#         return 1
#     else:
#         return (s_x - s_xy) / s_x

# def correlation_ratio(categories, measurements, nan_strategy=_REPLACE, nan_replace_value=_DEFAULT_REPLACE_VALUE):
#     """
#     Calculates the Correlation Ratio (sometimes marked by the greek letter Eta)
#     for categorical-continuous association.
#     Answers the question - given a continuous value of a measurement, is it
#     possible to know which category is it associated with?
#     Value is in the range [0,1], where 0 means a category cannot be determined
#     by a continuous measurement, and 1 means a category can be determined with
#     absolute certainty.
#     Wikipedia: https://en.wikipedia.org/wiki/Correlation_ratio
#     Parameters:
#     -----------
#     categories : list / NumPy ndarray / Pandas Series
#         A sequence of categorical measurements
#     measurements : list / NumPy ndarray / Pandas Series
#         A sequence of continuous measurements
#     nan_strategy : string, default = 'replace'
#         How to handle missing values: can be either 'drop' to remove samples
#         with missing values, or 'replace' to replace all missing values with
#         the nan_replace_value. Missing values are None and np.nan.
#     nan_replace_value : any, default = 0.0
#         The value used to replace missing values with. Only applicable when
#         nan_strategy is set to 'replace'.
#     Returns:
#     --------
#     float in the range of [0,1]
#     """
#     if nan_strategy == _REPLACE:
#         categories, measurements = replace_nan_with_value(
#             categories, measurements, nan_replace_value)
#     elif nan_strategy == _DROP:
#         categories, measurements = remove_incomplete_samples(
#             categories, measurements)
#     categories = convert(categories, 'array')
#     measurements = convert(measurements, 'array')
#     fcat, _ = pd.factorize(categories)
#     cat_num = np.max(fcat) + 1
#     y_avg_array = np.zeros(cat_num)
#     n_array = np.zeros(cat_num)
#     for i in range(0, cat_num):
#         cat_measures = measurements[np.argwhere(fcat == i).flatten()]
#         n_array[i] = len(cat_measures)
#         y_avg_array[i] = np.average(cat_measures)
#     y_total_avg = np.sum(np.multiply(y_avg_array, n_array)) / np.sum(n_array)
#     numerator = np.sum(
#         np.multiply(n_array, np.power(np.subtract(y_avg_array, y_total_avg),
#                                       2)))
#     denominator = np.sum(np.power(np.subtract(measurements, y_total_avg), 2))
#     if numerator == 0:
#         eta = 0.0
#     else:
#         eta = np.sqrt(numerator / denominator)
#     return eta

# def identify_nominal_columns(dataset):
#     return identify_columns_by_type(dataset, include=['object', 'category'])


# def identify_numeric_columns(dataset):
#     return identify_columns_by_type(dataset, include=['int64', 'float64'])



# def pairwise_cramers_v(categorical_df): 

# 	feature_names = categorical_df.columns
# 	results_df = pd.DataFrame(index=feature_names, columns=feature_names)

# 	for feat in feature_names: results_df.loc[feat, feat] = 1

# 	for feat1,feat2 in combinations(feature_names, 2): 
# 		corr = cramers_v(categorical_df[feat1], categorical_df[feat2])
# 		results_df.loc[feat1, feat2] = corr
# 		results_df.loc[feat2, feat1] = corr

# 	return results_df

# def pairwise_theils_u(categorical_df): 

# 	feature_names = categorical_df.columns
# 	results_df = pd.DataFrame(index=feature_names, columns=feature_names)

# 	for feat in feature_names: results_df.loc[feat, feat] = 1

# 	for feat1,feat2 in permutations(feature_names, 2): 
# 		results_df.loc[feat1, feat2] = theils_u(categorical_df[feat1], categorical_df[feat2])

# 	return results_df

# def pairwise_correlation_ratio(categories, measurements):
# 	"""
# 	Compare categorical with numerical data
# 	"""
# 	fcat, _ = pd.factorize(categories)
# 	cat_num = np.max(fcat)+1
# 	y_avg_array = np.zeros(cat_num)
# 	n_array = np.zeros(cat_num)
# 	for i in range(0,cat_num):
# 		cat_measures = measurements[np.argwhere(fcat == i).flatten()]
# 		n_array[i] = len(cat_measures)
# 		y_avg_array[i] = np.average(cat_measures)
# 	y_total_avg = np.sum(np.multiply(y_avg_array,n_array))/np.sum(n_array)
# 	numerator = np.sum(np.multiply(n_array,np.power(np.subtract(y_avg_array,y_total_avg),2)))
# 	denominator = np.sum(np.power(np.subtract(measurements,y_total_avg),2))
# 	if numerator == 0:
# 		eta = 0.0
# 	else:
# 		eta = np.sqrt(numerator/denominator)
# 	return eta


#----------------------------------------------------------------------------------------------------#
# Old functions 
#----------------------------------------------------------------------------------------------------#
def _correlate_feature_pairs(df1, df2, pairs=None, method="pearson", fast=True): 
	"""
	pairs (pd.DataFrame): two columns specifying rows in df1 and df2

	correlate_feature_pairs(pairs.as_df()[["gene_id", "peak_id"]].head(50), rna.vals, atac.vals)
	"""
	if pairs is None: 
		common_idx = df1.index & df2.index
		if len(df1) != len(df2): 
			logger.write("Inputs are of unequal length. Reporting the intersection.")
			df1 = df1.reindex(common_idx)
			df2 = df2.reindex(common_idx)
		pairs = pd.DataFrame.from_dict({
			"feature_1": common_idx, 
			"feature_2": common_idx
		})

	col1, col2 = pairs.columns[:2]
	r_col, p_col = "{}_r".format(method), "{}_pval".format(method)

	if fast:
		from pathos.multiprocessing import ProcessingPool as Pool
		from pathos.helpers import cpu_count
		n_cpu = cpu_count()

		# Group by first column (i.e. `gene_id`)
		pairs_grouped = pairs.groupby(col1)[col2]

		def apply_corr(data, method=method):
			name,s,df = data
			corrs = df.T.corrwith(s, method=method).rename(r_col).to_frame().reset_index()
			corrs[col1] = name
			return corrs  # `gene_id`, `peak_id`, `pearson_r`

		with Pool(n_cpu) as pool:
			# Generator with name of first column, series with that column, and dataframe associated with its pairs
			data = ((name, df1.loc[name], df2.loc[group]) for name, group in pairs_grouped)
			logger.write("Parallizing with {} cpus".format(n_cpu))
			results = pool.map(apply_corr, data)
		
		results = pd.concat(results)

	else: 
		def apply_corr(s):
			df1_vals = df1.loc[s.name]  # series from df1
			df2_vals = df2.loc[s]		# dataframe of df2 values from groupby
			return df2_vals.T.corrwith(df1_vals, method=method)

		groupby_obj = pairs.groupby(pairs.columns[0])[pairs.columns[1]]
		results = groupby_obj.apply(apply_corr).rename(r_col).to_frame().reset_index()

	# add p-value
	from scipy import stats
	n = df1.shape[1]
	tt = results[r_col] / np.sqrt((1 - results[r_col]**2) / (n - 2))
	results[p_col] = stats.t.sf(np.abs(tt), n-1) * 2

	# return results[[col1, col2, r_col, p_col]]
	return results



