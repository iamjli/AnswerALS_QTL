#!/usr/bin/env python3

import numpy as np
import pandas as pd

from pathlib import Path

from src import logger, hg38


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


def parse_tensorqtl_config_file(config_file): 
	"""
	Returns paths and params used to initialize a tensorqtl run.
	"""
	assert Path(config_file).is_file()

	with open(config_file, "r") as f: 
		params = { line.split("\t")[0]:line.rstrip("\n").split("\t")[1] for line in f.readlines() }

	tensorqtl_paths = {
		"results_dir": params.pop("results_dir"),
		"plink_prefix": params.pop("plink_prefix"),
		"phenotype_prefix": params.pop("phenotype_prefix"),
		"covariates_path": params.pop("covariates_path"),
	}

	return tensorqtl_paths, params


class TensorQTLRun: 
	
	def __init__(self, results_dir, plink_prefix, phenotype_prefix, covariates_path, params): 

		self.results_dir = Path(results_dir)
		self.plink_prefix = Path(plink_prefix)
		self.phenotype_prefix = Path(phenotype_prefix)
		self.covariates_path = Path(covariates_path)
		self.params = params

		self.results_dir.mkdir(exist_ok=True)
		self.config_file = self.results_dir / "config.txt"

		# Set input and output paths for tensorqtl by chrom batch
		self.chroms = hg38.chroms
		self._batch_input_paths, self._batch_output_paths = dict(), dict()
		for chrom in self.chroms: 
			self._batch_input_paths[chrom], self._batch_output_paths[chrom] = self._get_batch_tensorqtl_paths(chrom)

	@classmethod
	def load_config(cls, config_file): 
		"""Load from config"""
		tensorqtl_paths, params = parse_tensorqtl_config_file(config_file)
		return cls(**tensorqtl_paths, params=params)

	def save_config(self, overwrite=False): 

		if self.config_file.is_file() and not overwrite:
			logger.write("Config file already exists")
		else: 
			with open(self.config_file, "w") as f: 
				f.write("{}\t{}\n".format("results_dir", self.results_dir))
				f.write("{}\t{}\n".format("plink_prefix", self.plink_prefix))
				f.write("{}\t{}\n".format("phenotype_prefix", self.phenotype_prefix))
				f.write("{}\t{}\n".format("covariates_path", self.covariates_path))
				for key,val in self.params.items(): 
					f.write("{}\t{}\n".format(key, val))
			logger.write("Config file saved to:", self.config_file)


	#----------------------------------------------------------------------------------------------------#
	# Run tensorQTL in batches
	#----------------------------------------------------------------------------------------------------#
	def get_tensorqtl_cmd(self, mode, chrom): 
		"""Returns tensorqtl command."""
		return _tensorqtl_cmd(mode, **self._batch_input_paths[chrom], **self.params)

	def to_run(self, mode=None): 
		"""Displays all commands left to run."""
		if mode == None: 
			self.to_run(mode="cis")
			self.to_run(mode="cis_nominal")
		else:
			for chrom in self.chroms: 
				if not self._batch_output_paths[chrom][f"{mode}_results"].is_file():
					print(self.get_tensorqtl_cmd(mode, chrom))

	def _get_batch_tensorqtl_paths(self, chrom):
		"""Generates input and output paths for tensorqtl run by chrom batch."""
		input_paths = {
			"results_prefix": f"{self.results_dir}/{chrom}",
			"genotype_path": f"{self.plink_prefix}.{chrom}",
			"phenotype_path": f"{self.phenotype_prefix}.{chrom}.bed.gz",
			"covariates_path": f"{self.covariates_path}",
		}

		results_prefix = input_paths["results_prefix"]
		output_paths = {
			"cis_log": Path(f"{results_prefix}.tensorQTL.cis.log"),
			"cis_results": Path(f"{results_prefix}.cis_qtl.txt.gz"),
			"cis_nominal_log": Path(f"{results_prefix}.tensorQTL.cis_nominal.log"),
			"cis_nominal_results": Path(f"{results_prefix}.cis_qtl_pairs.{chrom[3:]}.parquet"),
		}
		
		return input_paths, output_paths

	#----------------------------------------------------------------------------------------------------#
	# Post-process cis batch results
	#----------------------------------------------------------------------------------------------------#
	def cis_processed_path(self, fdr):
		"""Results path to top QTL for each phenotype.""" 
		return self.results_dir / "cis_qtl.fdr_{}.txt.gz".format(str(fdr))

	def _import_cis_batch_result(self, path): 
		return pd.read_csv(path, sep="\t", index_col=0, compression="gzip")

	def process_cis_results(self, fdr): 
		"""Generates file with FDR thresholds."""
		cis_processed_path = self.cis_processed_path(fdr)
		if cis_processed_path.is_file(): 
			logger.write("cis file already exists. Exiting.")
			return 

		# 1. Load and merge all individual chromosome results
		cis_chrom_dfs = []
		for chrom in self.chroms: 
			cis_batch_path = self._batch_output_paths[chrom]["cis_results"]
			cis_chrom_dfs.append(self._import_cis_batch_result(cis_batch_path))
		cis_df = pd.concat(cis_chrom_dfs)

		# 2. Compute pval thresholds per phenotype
		cis_qval_df = calculate_qvalues(cis_df, fdr)

		# 3. Write results
		cis_qval_df.to_csv(cis_processed_path, sep="\t", index=True, header=True, compression="gzip") 
		logger.write("cis results written to:", str(cis_processed_path))

	#----------------------------------------------------------------------------------------------------#
	# Post-process cis nominal batch results
	#----------------------------------------------------------------------------------------------------#
	def cis_nominal_processed_path(self, fdr): 
		"""Results path to all sig QTLs using thresholds determined in leads file."""
		return self.results_dir / "cis_qtl_nominal_sig.fdr_{}.txt.gz".format(str(fdr))
		# if chrom is None: 
		# 	return self.results_dir / "cis_qtl_nominal_sig.fdr_{}.txt.gz".format(str(fdr))
		# else: 
		# 	return self.results_dir / "cis_qtl_nominal_sig.fdr_{}.{}.txt.gz".format(str(fdr), chrom)

	def _import_cis_nominal_batch_results(self, path, columns=None): 

		if columns is None: 
			columns = ["phenotype_id", "variant_id", "tss_distance", "pval_nominal", "slope"]
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
		thresholds = self.load_cis_results()["pval_nominal_threshold"]

		# 2. Load each cis nominal file and apply threshold filter
		cis_nominal_sig_dfs = []
		for chrom in self.chroms: 
			logger.update(f"Reading cis nominal results for {chrom}...")
			cis_nominal_batch_path = self._batch_output_paths[chrom]["cis_results"]
			cis_nominal_batch_df = self._import_cis_nominal_batch_results(cis_nominal_batch_path)

			cis_nominal_batch_sig_df = cis_nominal_batch_df.loc[cis_nominal_batch_df["pval_nominal"] <= cis_nominal_batch_df.map(thresholds)]
			cis_nominal_sig_dfs.append(cis_nominal_batch_sig_df)
		cis_nominal_sig_df = pd.concat(cis_nominal_sig_dfs, ignore_index=True)

		# 3. Write results
		cis_nominal_sig_df.to_csv(cis_nominal_processed_path, sep="\t", index=True, header=True, compression="gzip") 
		logger.write("cis results written to:", str(cis_nominal_processed_path))

	def pickle_cis_nominal(self): 
		"""
		Loading all cis_nominal results takes a long time to load, so in the rare occasion that we need them, 
		we can reduce loading time by loading them as a pickle. 
		"""
		pass

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
		cis_nominal_dfs = []
		for chrom in self.chroms: 
			logger.update(f"Reading cis nominal results for {chrom}...")
			cis_nominal_batch_path = self._batch_output_paths[chrom]["cis_results"]
			cis_nominal_dfs.append(self._import_cis_nominal_batch_results(cis_nominal_batch_path))
		return pd.concat(cis_nominal_dfs, ignore_index=True)


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
