#!/usr/bin/env python3

import pandas as pd
import numpy as np

from src import logger, RESULTS_PATHS, BASE_DIR

# tensorqtl libraries
logger.write("Importing torch may cause other basic packages to fail.")
import torch
from tensorqtl.core import impute_mean, calculate_corr
from tensorqtl.cis import calculate_association



class Residualizer(object):
	"""
	Based on `tensorqtl.core.Residualizer` but added support for dataframes.
	"""
	def __init__(self, C):
		"""
		C: samples x covariates
		"""
		self.C = C
		if isinstance(C, pd.DataFrame): 
			if not C.index.str.startswith("NEU").all(): 
				logger.write("Warning: check that input is indexed by sample name")
			C_t = torch.tensor(C.values, dtype=torch.float32).to("cpu") 
		elif isinstance(C, torch.Tensor): 
			C_t = C
		else: 
			logger.write("Must provide as dataframe or tensor.")

		# center and orthogonalize
		self.Q_t, _ = torch.qr(C_t - C_t.mean(0))
		self.dof = C_t.shape[0] - 2 - C_t.shape[1]

	@classmethod
	def load_rna(cls, covariate_path=RESULTS_PATHS["rna_covariates"]): 
		logger.write("Loading RNA covariates from: {}".format(covariate_path.relative_to(BASE_DIR)))
		covariates_df = pd.read_csv(covariate_path, sep="\t", index_col=0)
		return cls(covariates_df.T)

	@classmethod
	def load_atac(cls, covariate_path=RESULTS_PATHS["atac_covariates"]): 
		logger.write("Loading ATAC covariates from: {}".format(covariate_path.relative_to(BASE_DIR)))
		covariates_df = pd.read_csv(covariate_path, sep="\t", index_col=0)
		return cls(covariates_df.T)

	def transform(self, M, center=True):
		"""Residualize rows of M wrt columns of C. Does not necessarily need to be normalized."""

		is_df = isinstance(M, pd.DataFrame)
		M_t = torch.tensor(M.values, dtype=torch.float).to("cpu") if is_df else M

		if center:
			# center row means
			M0_t = M_t - M_t.mean(1, keepdim=True) 
		else:
			M0_t = M_t

		# the second term is the components of M that are explainable by Q. First projects into covariate space, then projects back
		# Note that normally the projection back would be Q_inverse, but because it's orthonormal, that's equal to Q^T
		M_t_transformed = M_t - torch.mm(torch.mm(M0_t, self.Q_t), self.Q_t.t())  # keep original mean

		if is_df: 
			return pd.DataFrame(M_t_transformed.numpy(), index=M.index, columns=M.columns)
		else: 
			return M_t_transformed


def QTL_pairwise(genotypes_df, phenotypes_df, covariates_df=None, report_maf=False, return_r_matrix=False): 
	"""
	Wrapper for `tensorqtl.core.calculate_corr` and reimplementation of `tensorqtl.cis.calculate_association`.
	Sample names must be axis 0 for each input (i.e. index for pd.Series and columns for pd.DataFrame)
	"""
	if isinstance(genotypes_df, pd.Series): 
		genotypes_df = genotypes_df.to_frame().T
	if isinstance(phenotypes_df, pd.Series): 
		phenotypes_df = phenotypes_df.to_frame().T

	assert genotypes_df.columns.equals(phenotypes_df.columns)

	# Prepare variables as torch tensors
	genotypes_t  = torch.tensor(genotypes_df.values, dtype=torch.float).to("cpu")
	phenotypes_t = torch.tensor(phenotypes_df.values, dtype=torch.float).to("cpu")
	impute_mean(genotypes_t)

	if covariates_df is not None: 
		assert genotypes_df.columns.equals(covariates_df.columns)
		residualizer = Residualizer(torch.tensor(covariates_df.T.values, dtype=torch.float32).to("cpu"))
		dof = residualizer.dof
	else: 
		residualizer = None
		dof = genotypes_t.shape[1] - 2

	# Compute pairwise correlations and associated stats
	r_nominal_t, genotype_var_t, phenotype_var_t = calculate_corr(genotypes_t, phenotypes_t, residualizer=residualizer, return_var=True)
	std_ratio_t = torch.sqrt(phenotype_var_t.reshape(1,-1) / genotype_var_t.reshape(-1,1))
	r_nominal_t = r_nominal_t.squeeze()
	r2_nominal_t = r_nominal_t.double().pow(2)
	slope_t = r_nominal_t * std_ratio_t.squeeze()
	tstat_t = r_nominal_t * torch.sqrt(dof / (1 - r2_nominal_t))
	slope_se_t = (slope_t.double() / tstat_t).float()

	if return_r_matrix: 
		return pd.DataFrame(r_nominal_t.numpy(), index=genotypes_df.index, columns=phenotypes_df.index)

	# Prepare results as dataframe. 
	genotype_ids, phenotype_ids = zip(*product(genotypes_df.index, phenotypes_df.index))
	results =  pd.DataFrame({
		"variant_id": genotype_ids,
		"phenotype_id": phenotype_ids, 
		"tstat": tstat_t.flatten().numpy(), 
		"slope": slope_t.flatten().numpy(),
		"slope_se": slope_se_t.flatten().numpy(),
		"r": r_nominal_t.flatten().numpy(),
		"r2": r2_nominal_t.flatten().numpy()
	})
	results["pval_nominal"] = 2*stats.t.cdf(-np.abs(results["tstat"]), dof)

	if report_maf: 
		# calculate MAF
		n2 = 2 * genotypes_t.shape[1] 				# total number of alleles
		af_t = genotypes_t.sum(1) / n2 				# allele frequency
		ix_t = af_t <= 0.5							
		maf_t = torch.where(ix_t, af_t, 1 - af_t)	# minor allele frequency
		# calculate MA samples and counts
		m = genotypes_t > 0.5
		a = m.sum(1).int()
		b = (genotypes_t < 1.5).sum(1).int()
		ma_samples_t = torch.where(ix_t, a, b)		# number of samples with a minor allele
		a = (genotypes_t * m.float()).sum(1).int()
		ma_count_t = torch.where(ix_t, a, n2-a)		# total number of minor alleles

		results["maf"]		= maf_t.flatten().numpy()
		results["ma_samples"] = ma_samples_t.flatten().numpy()
		results["ma_count"]   = ma_count_t.flatten().numpy()

	return results.sort_values("pval_nominal")

def QTL_diff(genotypes_df, phenotypes_df, covariates_df, condition_s, report_maf=True): 

	if isinstance(genotypes_df, pd.Series): 
		genotypes_df = genotypes_df.to_frame().T
	if isinstance(phenotypes_df, pd.Series): 
		phenotypes_df = phenotypes_df.to_frame().T

	als_qtl = QTL_pairwise(
		genotypes_df  = genotypes_df.loc[:,condition_s == "ALS"], 
		phenotypes_df = phenotypes_df.loc[:,condition_s == "ALS"], 
		covariates_df = covariates_df.loc[:,condition_s == "ALS"] if covariates_df is not None else None,
		report_maf = report_maf
	).drop(columns=["ma_count"], errors="ignore")

	ctr_qtl = QTL_pairwise(
		genotypes_df  = genotypes_df.loc[:,condition_s == "CTR"], 
		phenotypes_df = phenotypes_df.loc[:,condition_s == "CTR"], 
		covariates_df = covariates_df.loc[:,condition_s == "CTR"] if covariates_df is not None else None, 
		report_maf = report_maf
	).drop(columns=["ma_count"], errors="ignore")

	merged_qtl = als_qtl.merge(ctr_qtl, on=["variant_id", "phenotype_id"], suffixes=["_ALS", "_CTR"])
	# Compute differential statistics
	merged_qtl["z"] = (merged_qtl["slope_ALS"] - merged_qtl["slope_CTR"]) / (merged_qtl["slope_se_ALS"] + merged_qtl["slope_se_CTR"])
	merged_qtl["pval"] = stats.norm.sf(abs(merged_qtl["z"]))*2
	return merged_qtl

def QTL_by_pairs(G, omic_df, pairs, covariates_df, report_maf=False, condition_s=None):
	"""Iterate over either genotype or omic (rather than performing full pairwise)."""
	phen_id = pairs.columns[0]
	# Iterate over the feature that requires fewer QTL calls
	if pairs["variant_id"].nunique() < pairs[phen_id].nunique(): 
		logger.write("Iterating over variants...")
		grouped_pairs = pairs.groupby("variant_id")
	else: 
		logger.write("Iterating over phenotypes...")
		grouped_pairs = pairs.groupby(phen_id)

	results = []
	for i,feature in enumerate(grouped_pairs.groups.keys()): 
		logger.update("{} of {} features tested".format(i,grouped_pairs.ngroups))
		df = grouped_pairs.get_group(feature)
		genotypes_df = G.get_genotypes(df["variant_id"].unique())
		phenotypes_df = omic_df.loc[df[phen_id].unique()]

		if condition_s is None: 
			res = QTL_pairwise(genotypes_df, phenotypes_df, covariates_df, report_maf)
		else: 
			res = QTL_diff(genotypes_df, phenotypes_df, covariates_df, condition_s, report_maf)

		results.append(res)
	logger.flush()
	return pd.concat(results).rename(columns={"phenotype_id": phen_id})



######################################
## 			NORMALIZATIONS			##
######################################

# def annotate_symbol(df, col="gene_id"): 
# 	return df.assign(symbol=df[col].map(ENSG))

def grouped_ttest(df, grouping): 
	assert grouping.dtype == bool
	df1 = df.loc[:,grouping]
	df2 = df.loc[:,~grouping]
	return pd.DataFrame(stats.ttest_ind(df1, df2, axis=1), columns=df.index, index=["t", "pval"]).T


def pairwise_correlation(df1, df2):
	"""Pairwise correlation among rows from 2 dataframes with the same columns."""
	n = len(df1.columns)
	v1, v2 = df1.values.T, df2.values.T
	sums = np.multiply.outer(v2.sum(0), v1.sum(0))
	stds = np.multiply.outer(v2.std(0), v1.std(0))
	return pd.DataFrame((v2.T.dot(v1) - sums / n) / stds / n, df2.index, df1.index)


def correlate_feature_pairs(df1, df2, pairs=None, method="pearson", fast=True): 
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


def correlate_omics(omic1, omic2, window=0, method="pearson"): 
	from .bed import join_regions_by_window

	pairs = join_regions_by_window(omic1.pos_pr, omic2.pos_pr, window=window)
	pairs = pairs.as_df()[[omic1.phen_name, omic2.phen_name]]  # get `gene_id` and `peak_id` columns
	pairs = pairs.head(50)
	return correlate_feature_pairs(pairs, omic1.vals, omic2.vals, method=method)



