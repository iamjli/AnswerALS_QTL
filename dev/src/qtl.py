#!/usr/bin/env python3
import sys

import pandas as pd
import numpy as np

from scipy import interpolate, stats

from . import DATA, CHROMS, logger

import torch
from tensorqtl.core import impute_mean, calculate_corr
from tensorqtl.cis import calculate_association
from itertools import product


class Residualizer(object):
	"""
	Based on `tensorqtl.core.Residualizer` but added support for dataframes.
	"""
	def __init__(self, C):
		"""
		C: samples x covariates
		"""
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

	def transform(self, M, center=True):
		"""Residualize rows of M wrt columns of C"""

		is_df = isinstance(M, pd.DataFrame)
		M_t = torch.tensor(M.values, dtype=torch.float).to("cpu") if is_df else M

		if center:
			M0_t = M_t - M_t.mean(1, keepdim=True) 
		else:
			M0_t = M_t

		M_t_transformed = M_t - torch.mm(torch.mm(M0_t, self.Q_t), self.Q_t.t())  # keep original mean

		if is_df: 
			return pd.DataFrame(M_t_transformed.numpy(), index=M.index, columns=M.columns)
		else: 
			return M_t_transformed


def QTL_pairwise(genotypes_df, phenotypes_df, covariates_df=None, report_maf=False): 
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

	# Prepare results as dataframe. 
	genotype_ids, phenotype_ids = zip(*product(genotypes_df.index, phenotypes_df.index))
	results =  pd.DataFrame({
		"variant_id": genotype_ids,
		"phenotype_id": phenotype_ids, 
		"tstat": tstat_t.flatten().numpy(), 
		"slope": slope_t.flatten().numpy(),
		"slope_se": slope_se_t.flatten().numpy(),
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

		results["maf"]        = maf_t.flatten().numpy()
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
		covariates_df = covariates_df.loc[:,condition_s == "ALS"],
		report_maf = report_maf
	).drop(columns=["ma_count"])

	ctr_qtl = QTL_pairwise(
		genotypes_df  = genotypes_df.loc[:,condition_s == "CTR"], 
		phenotypes_df = phenotypes_df.loc[:,condition_s == "CTR"], 
		covariates_df = covariates_df.loc[:,condition_s == "CTR"], 
		report_maf = report_maf
	).drop(columns=["ma_count"])

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







# def test_diff_QTL(phenotype, genotypes, covariates, condition): 

# 	als_qtl = single_QTL(
# 		phenotype  = phenotype.loc[condition == "ALS"], 
# 		genotypes  = genotypes.loc[:,condition == "ALS"], 
# 		covariates = covariates.loc[:,condition == "ALS"]
# 	)

# 	ctr_qtl = single_QTL(
# 		phenotype  = phenotype.loc[condition == "CTR"], 
# 		genotypes  = genotypes.loc[:,condition == "CTR"], 
# 		covariates = covariates.loc[:,condition == "CTR"]
# 	)

# 	merged_qtl = pd.concat([als_qtl.add_suffix("_ALS"), ctr_qtl.add_suffix("_CTR")], axis=1)
# 	merged_qtl["z"] = (merged_qtl["slope_ALS"] - merged_qtl["slope_CTR"]) / (merged_qtl["slope_se_ALS"] + merged_qtl["slope_se_CTR"])
# 	merged_qtl["pval"] = stats.norm.sf(abs(merged_qtl["z"]))*2
# 	return merged_qtl



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

		# if indexed by phenotype, reset index, then set it again later
		indexed_by_phenotype = df.index.name == "phenotype_id"
		if indexed_by_phenotype:
			df = df.reset_index()

		# rename columns
		df = df.rename(columns={ "phenotype_id": self.phenotype_id, "tss_distance": self.distance_id })

		# add `chr` prefix if ATAC
		if self.omic == "atac": 
			if not df["peak_id"].iloc[0].startswith("chr"): 
				df["peak_id"] = "chr" + df["peak_id"]

		# annotate eQTLs with gene symbol
		if self.omic == "rna": 
			df["symbol"] = df["gene_id"].map(DATA.ENSG)

		# correct TSS distance
		if self.omic == "rna": 
			logger.write("Recalculating distances")
			# adjust sign for TSS on negative strand.
			strand = df["gene_id"].map(DATA.rna_metadata["strand"])
			df.loc[strand == "-", "tss_distance"] *= -1
			df["tss_distance"] = df["tss_distance"].astype(int)

			# tss = df["gene_id"].map(rna_metadata["tss"])
			# strand = df["gene_id"].map(rna_metadata["strand"])
			# variant_pos = df["variant_id"].map(self.R.rsid["pos"])
			# # distance of variant wrt gene. Negative values indicate upstream variant
			# distance = variant_pos - tss
			# distance.loc[strand == "-"] *= -1
			# df["tss_distance"] = distance.copy()
		else: 
			# May need to correct for midpoint shift
			pass

		# sort by pvalue
		# df = df.sort_values("pval_nominal")

		if indexed_by_phenotype: 
			df = df.set_index(self.phenotype_id)

		# return QTLResults object
		return QTLResults(df)

	@property
	def cis(self):
		if self._cis is None: 
			self._sig = None  # recalculate sig based on this dataframe
			df = import_cis_nominal_results(self.results_dir, columns=["phenotype_id", "variant_id", "tss_distance", "pval_nominal", "slope"], filter_id=True)
			self._cis = self._reformat_tensorqtl_outputs(df)
		return self._cis

	@property
	def leads(self):
		# to be used for `cis_independent` mode.
		if self._leads is None: 
			leads_file = self.results_dir / "cis_qtl.fdr_{}.txt.gz".format(str(self.fdr))
			if leads_file.exists(): 
				logger.write("Loading cis leads file with qvals from:", str(leads_file))
				df = pd.read_csv(leads_file, sep="\t", index_col=0, compression="gzip")
			else: 
				logger.write("Generating cis leads with qvals...")
				df = import_cis_results(self.results_dir)
				df = calculate_qvalues(df, fdr=self.fdr)
			self._leads = self._reformat_tensorqtl_outputs(df)
		return self._leads

	@property
	def ind(self):
		if self._ind is None: 
			df = import_cis_independent_results(self.results_dir)
			self._ind = self._reformat_tensorqtl_outputs(df)
		return self._ind
	
	@property
	def sig(self):
		if self._sig is None: 
			sig_file = self.results_dir / "cis_qtl_nominal_sig.fdr_{}.txt.gz".format(str(self.fdr))
			if sig_file.exists(): 
				logger.write("Loading sig QTLs from file from:", str(sig_file))
				sig = pd.read_csv(sig_file, sep="\t", compression="gzip")
			else: 
				df = self.cis.loc[ self.cis["pval_nominal"] <= self.cis[self.phenotype_id].map(self.leads["pval_nominal_threshold"]) ]
				sig = df.sort_values("pval_nominal")
			self._sig = QTLResults(sig)
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
		logger.write("Number of QTLs:", self.cis.shape[0])
		logger.write("Number of QTLs (FDR={}):".format(self.fdr), self.sig.shape[0])
		logger.write("Number of unique phenotypes:", len(self.phenotypes))
		logger.write("Number of unique phenotypes (FDR={}):".format(self.fdr), self.sig[self.phenotype_id].nunique())
		logger.write("Number of unique variants:", len(self.variants))
		logger.write("Number of unique variants (FDR={}):".format(self.fdr), self.sig["variant_id"].nunique())

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

	def write_sig_qtl(self, path=None, write=True): 
		"""
		Writes list of significant QTLs in parsed format. 
		"""
		if write: 
			if path is None: 
				path = self.results_dir / "cis_qtl_nominal_sig.fdr_{}.txt.gz".format(str(self.fdr))
			self.sig.to_csv(path, sep="\t", index=False, header=True, compression="gzip")



## Methods for QTL results

class QTLResults(pd.DataFrame): 

	_omic = None
	_phen_idx = None

	@property
	def _constructor(self):
		return QTLResults

	@property
	def omic(self):
		"""Determine omic type from column content."""
		if self._omic is None: 
			if ("gene_id" in self.columns) or (self.index.name == "gene_id"): 
				self._omic = "rna"
			elif ("peak_id" in self.columns) or (self.index.name == "peak_id"): 
				self._omic = "atac"
			else: 
				logger.write("Could not determine omic type.")
		return self._omic

	@property
	def phenotype_id(self):
		if self.omic == "rna": 
			return "gene_id"
		else: 
			return "peak_id"

	@property
	def phen_idx(self):
		if self._phen_idx is None: 
			logger.write("Building phenotype id index")
			self._phen_idx = {}
			for i,name in enumerate(self[self.phenotype_id].values): 
				if name not in self._phen_idx: 
					self._phen_idx[name] = [i,i]
				self._phen_idx[name][1] = i
		return self._phen_idx

	def fast_select(self, phenotype_id): 
		start_idx, end_idx = self.phen_idx[phenotype_id]
		return self.loc[start_idx:end_idx]

	@property
	def annotate(self): 
		if "gene_id" in self.columns: 
			return self.assign(symbol=self["gene_id"].map(DATA.ENSG))
		elif self.index.name == "gene_id": 
			return self.assign(symbol=self.index.map(DATA.ENSG))
		else: 
			logger.write("Cannot annotate with gene symbol.")

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
				logger.write("Assuming these are ENSG symbols...")
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

	return results[[col1, col2, r_col, p_col]]


def correlate_omics(omic1, omic2, window=0, method="pearson"): 
	from .bed import join_regions_by_window

	pairs = join_regions_by_window(omic1.pos_pr, omic2.pos_pr, window=window)
	pairs = pairs.as_df()[[omic1.phen_name, omic2.phen_name]]  # get `gene_id` and `peak_id` columns
	pairs = pairs.head(50)
	return correlate_feature_pairs(pairs, omic1.vals, omic2.vals, method=method)








######################################
## 		  TENSORQTL HELPERS			##
######################################

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
		sys.stdout.write("\rLoading `cis` results for {} ".format(chrom))
		sys.stdout.flush()
		try:
			cis_df = pd.read_csv("{}/{}.cis_qtl.txt.gz".format(results_dir, chrom), sep="\t", index_col=0, compression="gzip")
			all_df.append(cis_df)
		except FileNotFoundError:
			logger.write("Could not find cis results for {}".format(chrom))
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
	logger.write("{} pairs loaded".format(len(all_df)))
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
		sys.stdout.write("\rLoading `cis independent` results for {} ".format(chrom))
		sys.stdout.flush()
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
	verbose: log verbose messages? (default False)
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
		
			logger.write("qvalues pi0=%.3f, estimated proportion of null features " % pi0)
			if pi0 > 1: 
				logger.write("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
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
		logger.write('Computing q-values')
		logger.write('  * Number of phenotypes tested: {}'.format(res_df.shape[0]))
		logger.write('  * Correlation between Beta-approximated and empirical p-values: : {:.4f}'.format(
			stats.pearsonr(res_df['pval_perm'], res_df['pval_beta'])[0]))

	# calculate q-values
	qval, pi0 = qvalue(res_df['pval_beta'])
	if logger:
		logger.write('  * Proportion of significant phenotypes (1-pi0): {:.2f}'.format(1 - pi0))
		logger.write('  * QTL phenotypes @ FDR {:.2f}: {}'.format(fdr, np.sum(qval<=fdr)))

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
			logger.write('  * min p-value threshold @ FDR {}: {:.6g}'.format(fdr, pthreshold))
		# return pthreshold, res_df['beta_shape1'], res_df['beta_shape2']
		pthresholds = np.array([ stats.beta.ppf(pthreshold, row["beta_shape1"], row["beta_shape2"]) for _,row in res_df.iterrows() ])
		# pthresholds = stats.beta.ppf(pthreshold, res_df['beta_shape1'], res_df['beta_shape2'])
	else: 
		pthresholds = np.zeros(len(res_df))
		
	stats_df = pd.DataFrame({
		"qval": qval, 
		"pval_nominal_threshold": pthresholds
	}, index=res_df.index)
	
	return pd.concat([res_df, stats_df], axis=1)


