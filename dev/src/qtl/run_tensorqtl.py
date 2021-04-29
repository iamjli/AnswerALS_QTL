#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

from src import base_dir, logger
from src.load import hg38


_chroms = hg38.chroms

class TensorQTLManager: 

	def __init__(self, results_dir, plink_prefix, phenotype_path, covariates_path, params): 

		self.results_dir = Path(results_dir).resolve()
		self.plink_prefix = Path(plink_prefix).resolve()
		self.phenotype_path = Path(phenotype_path).resolve()
		self.covariates_path = Path(covariates_path).resolve()
		self.params = params

		self.config_file = self.results_dir / "config.txt"

		self._validate()

	#----------------------------------------------------------------------------------------------------#
	# Initializers
	#----------------------------------------------------------------------------------------------------#
	def _validate(self): 
		assert self.phenotype_path.is_file() and self.phenotype_path.name.endswith("bed.gz")
		assert self.covariates_path.exists()

	def initialize(self, overwrite_config=False): 
		"""Sets project directory structure."""
		self.results_dir.mkdir(exist_ok=True)

		for mode in ["cis", "cis_nominal"]: 
			output_dir = self.results_dir / f"{mode}_output"
			output_dir.mkdir(exist_ok=True)

		self.split_phenotypes_by_chrom()
		self.save_config(overwrite_config)

	def split_phenotypes_by_chrom(self): 
		"""Split phenotype file into chromosome chunks."""
		output_dir = self.results_dir / "phenotypes_by_chrom"
		output_dir.mkdir(exist_ok=True)

		# Get paths of phenotypes split by chrom
		phenotype_handle = self.phenotype_path.name.replace(".bed.gz", "")
		phenotype_chrom_paths = [output_dir / f"{phenotype_handle}.{chrom}.bed.gz" for chrom in _chroms]

		# If files have not already been generated
		if not all([path.is_file() for path in phenotype_chrom_paths]):

			logger.write("Splitting chromosome files...")
			phenotypes_df = pd.read_csv(self.phenotype_path, sep="\t", compression="gzip", low_memory=False, dtype={'#chr':str, '#Chr':str})

			# Iterate through chromosomes and output phenotypes
			for chrom,output_path in zip(_chroms, phenotype_chrom_paths):
				phenotypes_by_chrom_df = phenotypes_df[phenotypes_df["#chr"] == chrom[3:]]
				phenotypes_by_chrom_df.to_csv(output_path, sep="\t", index=False, compression="gzip")

	#----------------------------------------------------------------------------------------------------#
	# Manage tensorqtl runs
	#----------------------------------------------------------------------------------------------------#
	def get_tensorqtl_cmd(self, mode, chrom): 
		"""Returns tensorqtl command."""
		phenotype_handle = self.phenotype_path.name.replace(".bed.gz", "")
		input_paths = {
			"results_prefix": f"{self.results_dir}/{mode}_output/{chrom}",
			"genotype_path": f"{self.plink_prefix}.{chrom}",
			"phenotype_path": f"{self.results_dir}/phenotypes_by_chrom/{phenotype_handle}.{chrom}.bed.gz",
			"covariates_path": f"{self.covariates_path}",
		}
		return _tensorqtl_cmd(mode, **input_paths, **self.params)

	def to_run(self, mode=None): 
		"""Displays all commands left to run."""
		if mode == "cis": 
			for chrom in _chroms: 
				if not self.cis_output_path(chrom).is_file(): 
					print(self.get_tensorqtl_cmd("cis", chrom))
		elif mode == "cis_nominal": 
			for chrom in _chroms: 
				if not self.cis_nominal_output_path(chrom).is_file(): 
					print(self.get_tensorqtl_cmd("cis_nominal", chrom))
		elif mode == None: 
			self.to_run(mode="cis")
			self.to_run(mode="cis_nominal")

	#----------------------------------------------------------------------------------------------------#
	# Post-process cis batch results
	#----------------------------------------------------------------------------------------------------#
	def cis_output_path(self, chrom): 
		return self.results_dir / f"cis_output/{chrom}.cis_qtl.txt.gz"

	def cis_processed_path(self, fdr):
		"""Results path to top QTL for each phenotype.""" 
		return self.results_dir / "cis_qtl.fdr_{}.txt.gz".format(str(fdr))

	def _load_tensorqtl_cis_output(self, path): 
		return pd.read_csv(path, sep="\t", index_col=0, compression="gzip")

	def process_cis_results(self, fdr): 
		"""
		Generates file with FDR thresholds.
		"""
		cis_processed_path = self.cis_processed_path(fdr)
		if cis_processed_path.is_file(): 
			logger.write("cis file already exists. Exiting.")
			return 

		# 1. Load and merge all individual chromosome results
		cis_results = [self._load_tensorqtl_cis_output(self.cis_output_path(chrom)) for chrom in _chroms]
		cis_results = pd.concat(cis_results)

		# 2. Compute pval thresholds per phenotype
		cis_qval_df = calculate_qvalues(cis_results, fdr)

		# 3. Write results
		cis_qval_df.to_csv(cis_processed_path, sep="\t", index=True, header=True, compression="gzip") 
		logger.write("cis results written to:", str(cis_processed_path))

	#----------------------------------------------------------------------------------------------------#
	# Post-process cis nominal batch results
	#----------------------------------------------------------------------------------------------------#
	def cis_nominal_output_path(self, chrom): 
		return self.results_dir / f"cis_nominal_output/{chrom}.cis_qtl_pairs.{chrom[3:]}.parquet" 

	def cis_nominal_processed_path(self, fdr): 
		"""Results path to all sig QTLs using thresholds determined in leads file."""
		return self.results_dir / "cis_qtl_nominal_sig.fdr_{}.txt.gz".format(str(fdr))

	def _load_tensorqtl_cis_nominal_output(self, path, columns=None):
		if columns is None: columns = ["phenotype_id", "variant_id", "tss_distance", "pval_nominal", "slope"]
		return pd.read_parquet(path, columns=columns)

	def process_cis_nominal_results(self, fdr): 
		"""
		Generates file with all significant QTLs.
		"""
		cis_nominal_processed_path = self.cis_nominal_processed_path(fdr)
		if cis_nominal_processed_path.is_file(): 
			logger.write("cis nominal file already exists. Exiting.")
			return 

		# 1. Import cis nominal thresholds
		thresholds = self.load_cis_results(fdr)["pval_nominal_threshold"]

		# 2. Load each cis nominal file and apply threshold filter
		cis_nominal_results = []
		for chrom in _chroms: 
			logger.update(f"Reading cis nominal results for {chrom}...")
			path = self.cis_nominal_output_path(chrom)
			results_df = self._load_tensorqtl_cis_nominal_output(path)
			results_df = results_df.loc[results_df["pval_nominal"] <= results_df["phenotype_id"].map(thresholds)]
			cis_nominal_results.append(results_df)
		cis_nominal_results = pd.concat(cis_nominal_results, ignore_index=True)

		# 3. Write results
		cis_nominal_results.to_csv(cis_nominal_processed_path, sep="\t", index=True, header=True, compression="gzip") 
		logger.write("cis results written to:", str(cis_nominal_processed_path))

	def pickle_cis_nominal(self): 
		"""
		Loading all cis_nominal results takes a long time to load, so in the rare occasion that we need them, 
		we can reduce loading time by loading them as a pickle. 
		"""
		pass

	#----------------------------------------------------------------------------------------------------#
	# Save/load session
	#----------------------------------------------------------------------------------------------------#
	def save_config(self, overwrite=False):
		if self.config_file.is_file() and not overwrite:
			logger.write("Config file already exists")
		else: 
			with open(self.config_file, "w") as f: 
				f.write(f"results_dir\t{self.results_dir}\n")
				f.write(f"plink_prefix\t{self.plink_prefix}\n")
				f.write(f"phenotype_path\t{self.phenotype_path}\n")
				f.write(f"covariates_path\t{self.covariates_path}\n")
				for key,val in self.params.items(): 
					f.write(f"{key}\t{val}\n")
			logger.write("Config file saved to:", self.config_file)

	@classmethod
	def load_config(cls, config_file=None, dirname=None): 
		"""Load from config"""
		if config_file is None: 
			config_file = base_dir / "tensorqtl_runs" / dirname / "config.txt"

		tensorqtl_paths, params = _parse_tensorqtl_config_file(config_file)
		return cls(**tensorqtl_paths, params=params)

	#----------------------------------------------------------------------------------------------------#
	# Load fully processed results
	#----------------------------------------------------------------------------------------------------#
	def load_cis_results(self, fdr): 
		path = self.cis_processed_path(fdr)
		return pd.read_csv(path, sep="\t", index_col=0, compression="gzip")

	def load_cis_nominal_results(self, fdr): 
		path = self.cis_nominal_processed_path(fdr)
		return pd.read_csv(path, sep="\t", index_col=0, compression="gzip")

	def load_cis_nominal_results_unfiltered(self): 
		cis_nominal_output_dir = self.results_dir / "cis_nominal_output"
		cis_nominal_results = []
		for chrom in _chroms: 
			logger.update(f"Reading cis nominal results for {chrom}...")
			path = self.cis_nominal_output_path(chrom)
			results_df = self._load_tensorqtl_cis_nominal_output(path, chrom)
			cis_nominal_results.append(results_df)
		return pd.concat(cis_nominal_results, ignore_index=True)


def _tensorqtl_cmd(mode, **kwargs): 
	"""
	tensorqtl commands 
	"""
	cmd = (
		"time python3 -m tensorqtl "
		"{genotype_path} {phenotype_path} {results_prefix} "
		"--covariates {covariates_path} "
		"--maf_threshold {maf_threshold} "
		"--window {window} "
		"--search {search} "
		"--mode {mode} "
	).format(mode=mode, **kwargs)
	
	if mode == "cis_independent": 
		cmd += "--cis_output {cis_output_path}".format(kwargs["cis_output_path"])
	
	return cmd

def _parse_tensorqtl_config_file(config_file): 
	"""
	Returns paths and params used to initialize a tensorqtl run.
	"""
	assert Path(config_file).is_file()

	with open(config_file, "r") as f: 
		params = { line.split("\t")[0]:line.rstrip("\n").split("\t")[1] for line in f.readlines() }

	tensorqtl_paths = {
		"results_dir": params.pop("results_dir"),
		"plink_prefix": params.pop("plink_prefix"),
		"phenotype_path": params.pop("phenotype_path"),
		"covariates_path": params.pop("covariates_path"),
	}

	return tensorqtl_paths, params
	
#----------------------------------------------------------------------------------------------------#
# Post-processing: qvalue functions for computing cis nominal pval thresholds
#----------------------------------------------------------------------------------------------------#
def calculate_qvalues(res_df, fdr=0.05, qvalue_lambda=None, logger=False):
	"""Annotate permutation results with q-values, p-value threshold"""
	from scipy import interpolate, stats

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
	from scipy import interpolate

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
