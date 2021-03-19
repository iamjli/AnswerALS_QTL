#!/usr/bin/env python3
import sys

import pandas as pd
import numpy as np

import torch
from tensorqtl.core import Residualizer, impute_mean
from tensorqtl.cis import calculate_cis_nominal

from scipy import interpolate, stats

from . import CHROMS
from .bed import df_to_pr, pr_to_df



def single_QTL(phenotype, genotypes, covariates): 
	"""
	Arguments: 
		- phenotype (1-d np.array): sample
		- genotype (2-d np.array): genotype by sample
		- covariates (2-d np.array): covariate by sample
	"""
	assert np.array(phenotype).ndim == 1
	n_samples = len(phenotype)
	dof = n_samples - 2

	phenotype_t = torch.tensor(np.array(phenotype).reshape(1,n_samples), dtype=torch.float).to("cpu")
	genotypes_t = torch.tensor(np.array(genotypes).reshape(-1,n_samples), dtype=torch.float).to("cpu")
	impute_mean(genotypes_t)

	covariates_t = torch.tensor(np.array(covariates).T.reshape(n_samples,-1), dtype=torch.float).to("cpu")
	residualizer = Residualizer(covariates_t)

	res = calculate_cis_nominal(genotypes_t, phenotype_t, residualizer)

	tstat, slope, slope_se, maf, ma_samples, ma_count = [i.cpu().numpy() for i in res]
	df = pd.DataFrame({
	    'pval_nominal':2*stats.t.cdf(-np.abs(tstat), dof),
	    'slope':slope, 'slope_se':slope_se,
	    'tstat':tstat, 'maf':maf, 'ma_samples':ma_samples, 'ma_count':ma_count,
	}, index=genotypes.index)

	return df

def test_diff_QTL(phenotype, genotypes, covariates, condition): 

	als_qtl = single_QTL(
		phenotype  = phenotype.loc[condition == "ALS"], 
		genotypes  = genotypes.loc[:,condition == "ALS"], 
		covariates = covariates.loc[:,condition == "ALS"]
	)

	ctr_qtl = single_QTL(
		phenotype  = phenotype.loc[condition == "CTR"], 
		genotypes  = genotypes.loc[:,condition == "CTR"], 
		covariates = covariates.loc[:,condition == "CTR"]
	)

	merged_qtl = pd.concat([als_qtl.add_suffix("_ALS"), ctr_qtl.add_suffix("_CTR")], axis=1)
	merged_qtl["z"] = (merged_qtl["slope_ALS"] - merged_qtl["slope_CTR"]) / (merged_qtl["slope_se_ALS"] + merged_qtl["slope_se_CTR"])
	merged_qtl["pval"] = stats.norm.sf(abs(merged_qtl["z"]))*2
	return merged_qtl


	# genotype = np.array(genotype)
	# if genotype.ndim == 1: 
	# 	genotype = genotype.reshape(1,-1)
	# genotype_t = torch.tensor(genotype, dtype=torch.float).to("cpu")

	# phenotype = np.array(phenotype).reshape(1,-1)
	# phenotype_t = torch.tensor(phenotype, dtype=torch.float).to("cpu")

	# covariates = np.array(covariates).T
	# covariates_t = torch.tensor(covariates, dtype=torch.float).to("cpu")

	# residualizer = Residualizer(covariates_t)


class QTL:

	def __init__(self, results_dir, omic, fdr=0.05, initialize=False): 
		# NOTE 3/19/21: tss distances may be off
		# in original run, the start position was used for TSS for genes on negative strand
		# after this was corrected, the signs for calculated distances are flipped for those on negative strand

		assert (omic == "rna") or (omic == "atac")

		self.results_dir = results_dir
		self.omic = omic
		self._fdr = fdr

		if omic == "rna": 
			self.phenotype_id = "gene_id"
			self.distance_id = "tss_distance"
		else: 
			self.phenotype_id = "peak_id"
			self.distance_id = "peak_distance"

		self._cis 		 = None  # all variant-phenotype pairs
		self._leads 	 = None  # top associations along with FDR thresholds
		self._ind 		 = None  # cis independent results
		self._sig 		 = None  # significant cis results
		self._phenotypes = None  # unique phenotypes
		self._variants   = None  # unique variants

		if initialize: 
			self.sig
			self.leads

	@property
	def fdr(self):
		return self._fdr

	@fdr.setter
	def fdr(self, fdr): 
		self._fdr = fdr
		self._leads = None
		self._sig = None

	def _reformat_tensorqtl_outputs(self, df): 

		# rename index if indexed by `phenotype_id` (i.e. leads and ind results)
		if df.index.name == "phenotype_id":
			df.index.name = self.phenotype_id

		# rename columns
		df = df.rename(columns={ "phenotype_id": self.phenotype_id, "tss_distance": self.distance_id })

		# add `chr` prefix if ATAC
		if self.omic == "atac": 
			if df.index.name == "peak_id": 
				df.index = "chr" + df.index
			else: 
				df["peak_id"] = "chr" + df["peak_id"]

		# correct TSS distance

		# sort by pvalue
		df = df.sort_values("pval_nominal")

		# return QTLResults object
		return QTLResults(df)

	@property
	def cis(self):
		if self._cis is None: 
			self._cis = import_cis_nominal_results(self.results_dir, columns=["phenotype_id", "variant_id", "tss_distance", "pval_nominal", "slope"], filter_id=True)
			self._cis = self._reformat_tensorqtl_outputs(self._cis)
			# self._cis = self._cis.rename(columns={ "phenotype_id": self.phenotype_id, "tss_distance": self.distance_id })
			# if self.omic == "atac": 
			# 	self._cis[self.phenotype_id] = "chr" + self._cis[self.phenotype_id]
			# self._cis = QTLResults(self._cis)
		return self._cis

	@property
	def leads(self):
		# to be used for `cis_independent` mode.
		if self._leads is None: 
			self._leads = import_cis_results(self.results_dir)
			self._leads = calculate_qvalues(self._leads, fdr=self.fdr)
			self._leads = self._reformat_tensorqtl_outputs(self._leads)

			# phenotype_id needs to be treated differently since it's the index
			# self._leads.index.name = self.phenotype_id
			# self._leads = self._leads.rename(columns={ "tss_distance": self.distance_id })
			# if self.omic == "atac": 
			# 	self._leads.index = "chr" + self._leads.index
			# self._leads = QTLResults(self._leads)
			# self._leads = self._leads.sort_pval
		return self._leads

	@property
	def ind(self):
		if self._ind is None: 
			self._ind = import_cis_independent_results(self.results_dir)
			self._ind = self._reformat_tensorqtl_outputs(self._ind)
			# self._ind = self._ind.rename(columns={ "phenotype_id": self.phenotype_id, "tss_distance": self.distance_id })
			# if self.omic == "atac": 
			# 	self._ind[self.phenotype_id] = "chr" + self._ind[self.phenotype_id]
			# self._ind = QTLResults(self._ind)
			# self._ind = self._ind.sort_pval
		return self._ind
	
	@property
	def sig(self):
		if self._sig is None: 
			self._sig = self.cis.loc[ self.cis["pval_nominal"] <= self.cis[self.phenotype_id].map(self.leads["pval_nominal_threshold"]) ]
			self._sig = self._sig.sort_values("pval_nominal")
			# self._sig = QTLResults(self._sig)
			# self._sig = self._sig.sort_pval
		return self._sig

	@property
	def phenotypes(self):
		if self._phenotypes is None: 
			self._phenotypes = self.cis[self.phenotype_id].unique()
		return self._phenotypes

	@property
	def variants(self):
		if self._variants is None: 
			self._variants = self.cis["variant_id"].unique()
		return self._variants

	@property
	def sig_phenotypes(self):
		return self.sig[self.phenotype_id].unique()
	
	@property
	def sig_variants(self):
		return self.sig["variant_id"].unique()
	
	def info(self): 
		print("Number of QTLs:", self.cis.shape[0])
		print("Number of QTLs (FDR={}):".format(self.fdr), self.sig.shape[0])
		print("Number of unique phenotypes:", len(self.phenotypes))
		print("Number of unique phenotypes (FDR={}):".format(self.fdr), self.sig[self.phenotype_id].nunique())
		print("Number of unique variants:", len(self.variants))
		print("Number of unique variants (FDR={}):".format(self.fdr), self.sig["variant_id"].nunique())

	def write_cis_qtl(self, path=None, write=True): 
		"""
		Produces tensorqtl input for `cis_independent` mode.
			1. merge `cis` results (aka the top variants associated per genotype)
			2. append columns indicating FDR
		"""
		results = import_cis_results(self.results_dir)
		results = calculate_qvalues(results, fdr=self.fdr)

		if write: 
			if path is None: 
				path = self.results_dir / "cis_qtl.fdr_{}.txt.gz".format(str(self.fdr))
			results.to_csv(path, sep="\t", index=True, header=True, compression="gzip") 
			
		return results


## Process tensorQTL results

def import_cis_results(results_dir, chrom=None): 
	"""Imports results of tensorqtl's `cis` mode."""
	if chrom is None: 
		chroms = CHROMS
	elif type(chrom) is not list: 
		chroms = [chrom]
	else: 
		pass

	all_df = []
	for chrom in chroms: 
		try:
			cis_df = pd.read_csv("{}/{}.cis_qtl.txt.gz".format(results_dir, chrom), sep="\t", index_col=0, compression="gzip")
			all_df.append(cis_df)
		except FileNotFoundError:
			print("Could not find cis results for {}".format(chrom))
	return pd.concat(all_df)

def import_cis_nominal_results(results_dir, chrom=None, columns=None, filter_id=False): 
	"""Imports results of tensorqtl's `cis_nominal` mode."""
	if chrom is None: 
		chroms = CHROMS
	elif type(chrom) is not list: 
		chroms = [chrom]
	else: 
		pass

	all_df = []
	for i,chrom in enumerate(chroms): 
		sys.stdout.write("\rLoading `nominal cis` results for {} ".format(chrom))
		sys.stdout.flush()

		cis_nominal_df = pd.read_parquet("{}/{}.cis_qtl_pairs.{}.parquet".format(results_dir, chrom, chrom[3:]))

		if filter_id:  # remove rows where variant id is unknown
			cis_nominal_df = cis_nominal_df[cis_nominal_df["variant_id"] != "."]
		if columns is not None:  # grab specific columns
			cis_nominal_df = cis_nominal_df[columns]

		all_df.append(cis_nominal_df)
	all_df = pd.concat(all_df, ignore_index=True)
	print("{} pairs loaded".format(len(all_df)))
	return all_df

def import_cis_independent_results(results_dir, chrom=None): 
	"""Imports results of tensorqtl's `cis_independent` mode."""
	if chrom is None: 
		chroms = CHROMS
	elif type(chrom) is not list: 
		chroms = [chrom]
	else: 
		pass

	all_df = []
	for chrom in chroms: 
		cis_independent_df = pd.read_csv("{}/{}.cis_qtl.txt.gz.cis_independent_qtl.txt.gz".format(results_dir, chrom), sep="\t", compression="gzip")
		all_df.append(cis_independent_df)
	return pd.concat(all_df)

def compute_thresholds(results_dir, fdr=0.1): 
	cis_df = import_cis_results(results_dir)
	cis_df = calculate_qvalues(cis_df, fdr=fdr)
	return cis_df

def qvalue(pv, m=None, pi0=None):
	"""
	Estimates q-values from p-values
	Args
	=====
	m: number of tests. If not specified m = pv.size
	verbose: print verbose messages? (default False)
	lowmem: use memory-efficient in-place algorithm
	pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
		 For most GWAS this is not necessary, since pi0 is extremely likely to be
		 1

	Adapted from https://github.com/nfusi/qvalue/blob/master/qvalue/qvalue.py
	"""
	assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

	if m is None: m = len(pv) # set the number of tests
		
	if pi0 is None: 
		if len(pv) < 100: # skip interpolation if number of hypotheses is small
			pi0 = 1.0
		
		else: 
			lambdas = np.arange(0, 0.95, 0.01)
			pi0_arr = np.array([(pv > lam).sum()/m/(1-lam) for lam in lambdas])
		
			# fit natural cubic spline
			tck = interpolate.splrep(lambdas, pi0_arr, k=3)
			pi0 = interpolate.splev(lambdas[-1], tck)
		
			print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)
			if pi0 > 1: 
				print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
				pi0 = 1
		
	assert(pi0 >= 0 and pi0 <= 1)
	
	t = np.array(pv).copy()  # array of thresholds
	rank = t.argsort().argsort()
	fdr = pi0 * m * t / (rank + 1)  # FDR(t) at different thresholds
	qv = np.array([ fdr[t >= pval].min() for pval in pv ])  
	return qv, pi0

def calculate_qvalues(res_df, fdr=0.05, qvalue_lambda=None, logger=False):
	"""Annotate permutation results with q-values, p-value threshold"""

	if logger: 
		print('Computing q-values')
		print('  * Number of phenotypes tested: {}'.format(res_df.shape[0]))
		print('  * Correlation between Beta-approximated and empirical p-values: : {:.4f}'.format(
			stats.pearsonr(res_df['pval_perm'], res_df['pval_beta'])[0]))

	# calculate q-values
	qval, pi0 = qvalue(res_df['pval_beta'])
	if logger:
		print('  * Proportion of significant phenotypes (1-pi0): {:.2f}'.format(1 - pi0))
		print('  * QTL phenotypes @ FDR {:.2f}: {}'.format(fdr, np.sum(qval<=fdr)))

	# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
	pval_pass = res_df.loc[qval<=fdr, 'pval_beta'].sort_values() # series of significant p-values at FDR threshold
	pval_fail = res_df.loc[qval>fdr, 'pval_beta'].sort_values()  # series of non-significant p-values
	
	if len(pval_pass) > 0: # make sure there is at least one significant phenotype
		# get p-value threshold by finding the boundary corresponding to FDR
		if len(pval_fail) > 0: 
			pthreshold = (pval_pass.values[-1] + pval_fail.values[0]) / 2
		else: 
			pthreshold = pval_pass.values[-1]
		if logger:
			print('  * min p-value threshold @ FDR {}: {:.6g}'.format(fdr, pthreshold))
		pthresholds = stats.beta.ppf(pthreshold, res_df['beta_shape1'], res_df['beta_shape2'])
	else: 
		pthresholds = np.zeros(len(res_df))
		
	stats_df = pd.DataFrame({
		"qval": qval, 
		"pval_nominal_threshold": pthresholds
	}, index=res_df.index)
	
	return pd.concat([res_df, stats_df], axis=1)


## Methods for QTL results

class QTLResults(pd.DataFrame): 

	@property
	def _constructor(self):
		return QTLResults

	@property
	def omic(self):
		if ("gene_id" in self.columns) or (self.index.name == "gene_id"): 
			return "rna"
		elif ("peak_id" in self.columns) or (self.index.name == "peak_id"): 
			return "atac"
		else: 
			print("Could not determine omic type.")

	# @property
	# def annotate(self): 
	# 	assert self.omic == "rna"
	# 	if "gene_id" in self.columns: 
	# 		return self.assign(symbol=self["gene_id"].map(ENSG))
	# 	elif self.index.name == "gene_id": 
	# 		return self.assign(symbol=self.index.map(ENSG))
	# 	else: 
	# 		print("Cannot annotate with gene symbol.")

	@property
	def sort_pval(self): 
		return self.sort_values("pval_nominal")

	@property
	def sort_distance(self):
		if self.omic == "rna": 
			return self.iloc[self["tss_distance"].abs().argsort()]
		else: 
			return self.iloc[self["peak_distance"].abs().argsort()]

	# @property
	# def sort_tss_dist(self):
	# 	return self.iloc[self["tss_distance"].abs().argsort()]

	# @property
	# def sort_peak_dist(self):
	# 	return self.iloc[self["peak_distance"].abs().argsort()]

	def select(self, var=None, gene_id=None, peak_id=None, variant_id=None, symbol=None): 

		if var is not None: # Infer variable type
			assert (gene_id is None) & (peak_id is None) & (variant_id is None) & (symbol is None)
			if type(var) == str: var = [var]
			if var[0].startswith("ENSG"): 
				gene_id = var
			elif var[0].startswith("chr"): 
				peak_id = var
			elif var[0].startswith("rs"): 
				variant_id = var
			else: 
				print("Assuming these are ENSG symbols...")
				symbol = var

		if gene_id is not None: 
			if type(gene_id) == str: gene_id = [gene_id]
			if "gene_id" in self.columns:
				self = self[self["gene_id"].isin(gene_id)]
			else: 
				self = self.loc[ensg]
		if peak_id is not None: 
			if type(peak_id) == str: peak_id = [peak_id]
			if "peak_id" in self.columns: 
				self = self[self["peak_id"].isin(peak_id)]
			else: 
				self = self.loc[peak_id]
		if variant_id is not None: 
			if type(variant_id) == str: variant_id = [variant_id]
			self = self[self["variant_id"].isin(variant_id)]
		if symbol is not None: 
			if type(symbol) == str: symbol = [symbol]
			if "symbol" in self.columns:
				self = self[self["symbol"].isin(symbol)]
			else: 
				self = self.loc[ensg]

		return self

	def filter(self, pval=None, distance=None): 

		if pval is not None: 
			self = self[self["pval_nominal"] <= pval]

		if distance is not None: 
			if self.omic == "rna": 
				self = self[self["tss_distance"].abs() <= distance]
			else: 
				self = self[self["peak_distance"].abs() <= distance]

		return self

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


def correlate_feature_pairs(pairs, df1, df2, method="pearson", fast=True): 
	"""
	pairs (pd.DataFrame): two columns specifying rows in df1 and df2

	correlate_feature_pairs(pairs.as_df()[["gene_id", "peak_id"]].head(50), rna.vals, atac.vals)
	"""
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
			print("Parallizing with {} cpus".format(n_cpu))
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

	return results[[col1, col2, r_col, p_col]]


def correlate_omics(omic1, omic2, window=0, method="pearson"): 
	from .bed import join_regions_by_window

	pairs = join_regions_by_window(omic1.pos_pr, omic2.pos_pr, window=window)
	pairs = pairs.as_df()[[omic1.phen_name, omic2.phen_name]]  # get `gene_id` and `peak_id` columns
	pairs = pairs.head(50)
	return correlate_feature_pairs(pairs, omic1.vals, omic2.vals, method=method)






# def tensorqtl_post(results_dir, fdr=0.05, write=True): 
# 	"""
# 	Ingests `cis` and `cis_nominal` results. Performs the following: 
# 	 - merges `cis` results and adds columns 
# 	 	 - `qval`: q-value for each top QTL
# 	 	 - `pval_nominal_threshold`: p-value threshold for `cis_nominal` results
# 	 - merges and filters `cis_nominal results`
# 	"""

# 	# Merge all lead cis results 
# 	cis_df = import_cis_results(results_dir)

# 	# Compute FDR and nominal p-value thresholds
# 	cis_df = calculate_qvalues(cis_df, fdr=fdr)

# 	# Merge and filter all variant-gene pairs
# 	sig_qtls_df = []
# 	for chrom in CHROMS: 
# 		cis_nominal_df = import_cis_nominal_results(results_dir, chrom=chrom)

# 		# Filter for significant pairs by checking if nominal pvalue is less than the threshold for the corresponding gene
# 		sig_rows = cis_nominal_df["pval_nominal"] <= cis_nominal_df["phenotype_id"].map(cis_df["pval_nominal_threshold"])

# 		print("--", chrom)
# 		print("# total pairs:", len(cis_nominal_df))
# 		print("# significant (FDR < {}):".format(str(fdr)), sig_rows.sum())
# 		sig_qtls_df.append(cis_nominal_df.loc[sig_rows])
		
# 	sig_qtls_df = pd.concat(sig_qtls_df).reset_index(drop=True)

# 	if write: 
# 		cis_df.to_csv("{}/cis_QTL.qval.txt".format(results_dir), sep="\t", index=True, header=True)
# 		sig_qtls_df.to_csv("{}/cis_nominal_QTL.FDR_{}.txt".format(results_dir, str(fdr)), sep="\t", index=False, header=True)

# 	return cis_df, sig_qtls_df



# def import_all_snps(results_dir, chrom=None): 
# 	"""Gets all SNPs used in `cis_nominal` mode."""
# 	if chrom is None: 
# 		chroms = CHROMS
# 		print("Importing all chromosomes. Could take a while.")
# 	elif type(chrom) is not list: 
# 		chroms = [chrom]
# 	else: 
# 		pass

# 	from pathos.multiprocessing import ProcessingPool as Pool
# 	def get_snps(chrom): 
# 		cis_nominal_df = import_cis_nominal_results(results_dir, chrom=chrom)
# 		return cis_nominal_df["variant_id"].unique().tolist()

# 	with Pool() as pool:
# 		results = pool.map(get_snps, chroms)
# 	return flatten(results)
