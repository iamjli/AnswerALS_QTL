#!/usr/bin/env python3

from itertools import product

import numpy as np
import pandas as pd
from scipy import stats
import torch
from tensorqtl import core

from src import logger


def QTL_pairwise(genotypes_df, phenotypes_df, residualizer=None, report_maf=False, return_r_matrix=False): 
	"""
	Wrapper for `tensorqtl.core.calculate_corr` and reimplementation of `tensorqtl.cis.calculate_association`.
	Sample names must be axis 0 for each input (i.e. index for pd.Series and columns for pd.DataFrame)
	"""
	if isinstance(genotypes_df, pd.Series): genotypes_df = genotypes_df.to_frame().T
	if isinstance(phenotypes_df, pd.Series): phenotypes_df = phenotypes_df.to_frame().T

	assert genotypes_df.columns.equals(phenotypes_df.columns)

	# Prepare variables as torch tensors
	genotypes_t  = torch.tensor(genotypes_df.values, dtype=torch.float).to("cpu")
	phenotypes_t = torch.tensor(phenotypes_df.values, dtype=torch.float).to("cpu")
	core.impute_mean(genotypes_t)

	dof = genotypes_t.shape[1] - 2 if residualizer is None else residualizer.dof

	# Compute pairwise correlations and associated stats
	r_nominal_t, genotype_var_t, phenotype_var_t = core.calculate_corr(genotypes_t, phenotypes_t, residualizer=residualizer, return_var=True)
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


#--------------------------------------------------------------------------------------------------#

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



