#!/usr/bin/env python3
import sys

import pandas as pd
import numpy as np

from scipy import interpolate, stats

from pathos.multiprocessing import ProcessingPool
from itertools import product
from pathlib import Path

from src import data, BASE_DIR, RESULTS_PATHS, CHROMS, logger



class QTL:

	def __init__(self, results_dir, omic, fdr=0.05, initialize=False): 
		# NOTE 3/19/21: tss distances may be off
		# in original run, the start position was used for TSS for genes on negative strand
		# after this was corrected, the signs for calculated distances are flipped for those on negative strand

		assert (omic == "rna") or (omic == "atac")

		self.results_dir = results_dir
		self.omic = omic
		self._fdr = fdr

		self._config_file = self.results_dir / "config.txt"
		self._tensorqtl_run = TensorQTLRun.load_config(self._config_file)
		assert self._tensorqtl_run is not None

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

	@classmethod
	def load_rna(cls, results_dir=RESULTS_PATHS["rna_results_dir"], **kwargs): 
		return cls(results_dir=results_dir, omic="rna", **kwargs)

	@classmethod
	def load_atac(cls, results_dir=RESULTS_PATHS["atac_results_dir"], **kwargs): 
		return cls(results_dir=results_dir, omic="atac", **kwargs)

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
			df["symbol"] = df["gene_id"].map(data.ENSG)

		# correct TSS distance
		if self.omic == "rna": 
			logger.write("Recalculating distances")
			# adjust sign for TSS on negative strand.
			strand = df["gene_id"].map(data.rna_metadata["strand"])
			df.loc[strand == "-", "tss_distance"] *= -1
			df["tss_distance"] = df["tss_distance"].astype(int)
		else: 
			# May need to correct for midpoint shift
			pass

		if indexed_by_phenotype: 
			df = df.set_index(self.phenotype_id)

		# return QTLResults object
		return QTLResults(df)

	@property
	def cis(self):
		if self._cis is None: 
			_cis = self._tensorqtl_run.cis_nominal_qtl_results_unfiltered()
			self._cis = self._reformat_tensorqtl_outputs(_cis)
		return self._cis

	@property
	def leads(self):
		# to be used for `cis_independent` mode.
		if self._leads is None: 
			_leads = self._tensorqtl_run.cis_qtl_results(self.fdr, write=False)
			self._leads = self._reformat_tensorqtl_outputs(_leads)
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
			_sig = self._tensorqtl_run.cis_nominal_qtl_results(fdr=self.fdr, write=False)
			self._sig = self._reformat_tensorqtl_outputs(_sig)
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

	# def write_cis_qtl(self, path=None, write=True): 
	# 	"""
	# 	Produces tensorqtl input for `cis_independent` mode.
	# 		1. merge `cis` results (aka the top variants associated per genotype)
	# 		2. append columns indicating FDR
	# 	"""
	# 	results = import_cis_results(self.results_dir)
	# 	results = calculate_qvalues(results, fdr=self.fdr)

	# 	if write: 
	# 		if path is None: 
	# 			path = self.results_dir / "cis_qtl.fdr_{}.txt.gz".format(str(self.fdr))
	# 		results.to_csv(path, sep="\t", index=True, header=True, compression="gzip") 

	# 	return results

	# def write_sig_qtl(self, path=None, write=True): 
	# 	"""
	# 	Writes list of significant QTLs in parsed format. 
	# 	"""
	# 	if write: 
	# 		if path is None: 
	# 			path = self.results_dir / "cis_qtl_nominal_sig.fdr_{}.txt.gz".format(str(self.fdr))
	# 		self.sig.to_csv(path, sep="\t", index=False, header=True, compression="gzip")



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
			return self.assign(symbol=self["gene_id"].map(data.ENSG))
		elif self.index.name == "gene_id": 
			return self.assign(symbol=self.index.map(data.ENSG))
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


######################################
## 		  TENSORQTL HELPERS			##
######################################


def tensorqtl_cmd(mode, **kwargs): 
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

class TensorQTLRun: 
	
	def __init__(self, results_dir, plink_prefix, phenotype_prefix, covariates_path, params): 
		
		self.results_dir = Path(results_dir)
		self.plink_prefix = Path(plink_prefix)
		self.phenotype_prefix = Path(phenotype_prefix)
		self.covariates_path = Path(covariates_path)
		self.params = params

		self.results_dir.mkdir(exist_ok=True)
		
		self.CHROMS = ["chr"+str(i) for i in range(1,23)] + ["chrX", "chrY"]
		
		self.input_paths, self.output_paths = dict(), dict()
		for chrom in self.CHROMS: 
			self.input_paths[chrom], self.output_paths[chrom] = self.get_tensorqtl_paths(chrom)

		self.config_file = self.results_dir / "config.txt"

	@classmethod
	def load_config(cls, config_file):
		if Path(config_file).exists():
			with open(config_file, "r") as f: 
				params = { line.split("\t")[0]:line.rstrip("\n").split("\t")[1] for line in f.readlines() }
			results_dir = params.pop("results_dir")
			plink_prefix = params.pop("plink_prefix")
			phenotype_prefix = params.pop("phenotype_prefix")
			covariates_path = params.pop("covariates_path")
			return cls(results_dir, plink_prefix, phenotype_prefix, covariates_path, params)
		else: 
			logger.write("Could not find config file.")
			return None

	def get_tensorqtl_paths(self, chrom):
		
		input_paths = {
			"results_prefix":  "{}/{}".format(self.results_dir, chrom), 
			"genotype_path":   "{}.{}".format(self.plink_prefix, chrom),
			"phenotype_path":  "{}.{}.bed.gz".format(self.phenotype_prefix, chrom), 
			"covariates_path": "{}".format(self.covariates_path),
		}
		
		output_paths = {
			"cis_log":			    Path("{results_prefix}.tensorQTL.cis.log".format(**input_paths)), 
			"cis_results":		    Path("{results_prefix}.cis_qtl.txt.gz".format(**input_paths)), 
			"cis_nominal_log":	    Path("{results_prefix}.tensorQTL.cis_nominal.log".format(**input_paths)),
			"cis_nominal_results":  Path("{results_prefix}.cis_qtl_pairs.{chrom_tag}.parquet".format(chrom_tag=chrom[3:], **input_paths)),
		}
		
		return input_paths, output_paths
			
	def cis_output_path(self, fdr): 
		return self.results_dir / "cis_qtl.fdr_{}.txt.gz".format(str(fdr))

	def cis_nominal_output_path(self, fdr, chrom=None): 
		if chrom is None: 
			return self.results_dir / "cis_qtl_nominal_sig.fdr_{}.txt.gz".format(str(fdr))
		else: 
			return self.results_dir / "cis_qtl_nominal_sig.fdr_{}.{}.txt.gz".format(str(fdr), chrom)
	
	def get_tensorqtl_cmd(self, mode, chrom): 
		return tensorqtl_cmd(mode, **self.input_paths[chrom], **self.params)
		
	def get_result_path(self, mode, chrom): 
		output_paths = self.output_paths[chrom]
		if mode == "cis": 
			return output_paths["cis_results"]
		elif mode == "cis_nominal": 
			return output_paths["cis_nominal_results"]
		else: 
			pass
		
	def to_run(self, mode=None): 
		for chrom in self.CHROMS: 
			if not self.get_result_path(mode, chrom).exists(): 
				print(self.get_tensorqtl_cmd(mode, chrom))
		
		if mode == None: 
			self.to_run(mode="cis")
			self.to_run(mode="cis_nominal")
			
	def finished(self, mode):
		for chrom in self.CHROMS: 
			if not self.get_result_path(chrom, mode).exists(): 
				return False
		else: 
			return True

	@staticmethod
	def _import_cis_results(path): 
		return pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
	
	@staticmethod
	def import_cis_results(paths, concat=True, pathos=False): 

		if pathos: 
			with ProcessingPool(8) as pool: 
				results = pool.map(TensorQTLRun._import_cis_results, paths)
		else: 
			results = [TensorQTLRun._import_cis_results(path) for path in paths]
		if concat: 
			return pd.concat(results)
		else: 
			return results

	@staticmethod
	def _import_cis_nominal_results(path, thresholds=None): 

		cis_nominal = pd.read_parquet(path, columns=["phenotype_id", "variant_id", "tss_distance", "pval_nominal", "slope"])
		if thresholds is not None: 
			cis_nominal = cis_nominal.loc[cis_nominal["pval_nominal"] <= cis_nominal["phenotype_id"].map(thresholds)]

		return cis_nominal

	@staticmethod
	def import_cis_nominal_results(paths=None, thresholds=None, concat=True, pathos=False): 

		logger.write("Loading all cis results...")
		if pathos: 
			cis_nominal_paths = [(path, thresholds) for path in paths]
			with ProcessingPool(8) as pool: 
				results = pool.map(lambda args: TensorQTLRun._import_cis_nominal_results(*args), cis_nominal_paths)
			logger.write("Finished loading individual parquets")
		else: 
			results = [TensorQTLRun._import_cis_nominal_results(path, thresholds) for path in paths]

		if concat:
			logger.write("Merging individual parquets")
			return pd.concat(results, ignore_index=True)
		else: 
			return results
	
	def cis_qtl_results(self, fdr, write=False): 

		path = self.cis_output_path(fdr)

		if path.exists(): 
			logger.write("Loading cis results from existing file.")
			results = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
		else: 
			cis_paths = (output_paths["cis_results"] for _,output_paths in self.output_paths.items())
			results = TensorQTLRun.import_cis_results(cis_paths)
			results = calculate_qvalues(results, fdr)
			if write:
				results.to_csv(path, sep="\t", index=True, header=True, compression="gzip") 

		return results

	def cis_nominal_qtl_results(self, fdr, write=False): 

		path = self.cis_nominal_output_path(fdr)

		if path.exists(): 
			logger.write("Loading cis_nominal results from existing file.")
			results = pd.read_csv(path, sep="\t", compression="gzip")
		else: 
			cis_nominal_paths = [output_paths["cis_nominal_results"] for _,output_paths in self.output_paths.items()]
			thresholds = self.cis_qtl_results(fdr, write=False)
			results = TensorQTLRun.import_cis_nominal_results(cis_nominal_paths, thresholds["pval_nominal_threshold"], concat=True)
			if write: 
				results.to_csv(path, sep="\t", index=False, header=True, compression="gzip") 

		return results

	def cis_nominal_qtl_results_unfiltered(self): 
		cis_nominal_paths = (output_paths["cis_nominal_results"] for _,output_paths in self.output_paths.items())
		results = TensorQTLRun.import_cis_nominal_results(cis_nominal_paths, concat=True)
		return results

	def save_config(self, overwrite=False): 

		path = self.config_file

		if (not path.exists()) | overwrite: 
			with open(path, "w") as f: 
				f.write("{}\t{}\n".format("results_dir", self.results_dir))
				f.write("{}\t{}\n".format("plink_prefix", self.plink_prefix))
				f.write("{}\t{}\n".format("phenotype_prefix", self.phenotype_prefix))
				f.write("{}\t{}\n".format("covariates_path", self.covariates_path))
				for key,val in self.params.items(): 
					f.write("{}\t{}\n".format(key, val))
			logger.write("Config file saved to:", path)


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

# def compute_thresholds(results_dir, fdr=0.1): 
# 	cis_df = import_cis_results(results_dir)
# 	cis_df = calculate_qvalues(cis_df, fdr=fdr)
# 	return cis_df

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










