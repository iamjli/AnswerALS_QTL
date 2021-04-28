#!/usr/bin/env python3

import numpy as np
import pandas as pd

from src import base_dir, logger
from src.load import aals, hg38
from src.qtl import preprocess, residualize, run_tensorqtl
from src.utils.regions import Regions


_omic_tag = {"atac": "atac", "rna": "rna", "erna": "rna", "srna": "rna"}  

class FeatureData: 
	"""Accessor for feature counts data and QTLs."""

	def __init__(self, feature, counts_prefix=None, qtl_dir=None, qtl_fdr=0.05): 

		self._feature = feature
		self._omic = _omic_tag[self._feature]

		if counts_prefix: self.initialize_counts(counts_prefix)
		if qtl_dir: self.initialize_qtl(qtl_dir, qtl_fdr)

	def initialize_counts(self, counts_prefix): 

		# Set counts data paths
		self._counts_prefix = counts_prefix
		self._counts_data_paths = {
			"metadata": base_dir / f"tensorqtl_runs/phenotypes_210422/{counts_prefix}.metadata.txt.gz", 
			"counts": base_dir / f"tensorqtl_runs/phenotypes_210422/{counts_prefix}.counts.txt.gz", 
			"tpm": base_dir / f"tensorqtl_runs/phenotypes_210422/{counts_prefix}.tpm.txt.gz", 
			"tmm": base_dir / f"tensorqtl_runs/phenotypes_210422/{counts_prefix}.tmm.txt.gz", 
			"gtex": base_dir / f"tensorqtl_runs/phenotypes_210422/{counts_prefix}.gtex.txt.gz", 
			"tmm_factors": base_dir / f"tensorqtl_runs/phenotypes_210422/{counts_prefix}.tmm_factors.txt.gz", 
		}
		for key,path in self._counts_data_paths.items(): 
			assert path.is_file(), f"{path} does not exist"

		# Initialize counts property attributes
		self._regions, self._lengths = [None] * 2
		self._counts, self._tpm, self._tmm, self._gtex = [None] * 4
		self._mask = None
		self._tmm_factors = None

	def initialize_qtl(self, qtl_dir, qtl_fdr): 

		# Set QTL data paths
		self._qtl_dir, self._qtl_fdr = qtl_dir, qtl_fdr
		self._qtl_config_path = base_dir / f"tensorqtl_runs/{qtl_dir}/config.txt"
		assert self._qtl_config_path.is_file(), f"{self._qtl_config_path} does not exist"

		# Initialize QTL property attributes
		self._qtl_accessor = run_tensorqtl.TensorQTLManager.load_config(self._qtl_config_path)
		self._qtl_sig, self._qtl_leads, self._qtl_all = [None] * 3

		self._covariates_path = self._qtl_accessor.covariates_path
		assert self._covariates_path.is_file(), f"{self._covariates_path} does not exist"
		self._residualizer = None

	#----------------------------------------------------------------------------------------------------#
	# Access counts data
	#----------------------------------------------------------------------------------------------------#
	@property
	def regions(self):
		if self._regions is None: self._regions, self._lengths = preprocess.load_counts_metadata(self._counts_data_paths["metadata"])
		return self._regions

	@property
	def lengths(self):
		if self._lengths is None: self._regions, self._lengths = preprocess.load_counts_metadata(self._counts_data_paths["metadata"])
		return self._lengths
	
	@property
	def counts(self):
		if self._counts is None: self._counts = pd.read_csv(self._counts_data_paths["counts"], sep="\t", index_col=0)
		return self._counts

	@property
	def tpm(self):
		if self._tpm is None: self._tpm = pd.read_csv(self._counts_data_paths["tpm"], sep="\t", index_col=0)
		return self._tpm

	@property
	def tmm(self):
		if self._tmm is None: self._tmm = pd.read_csv(self._counts_data_paths["tmm"], sep="\t", index_col=0)
		return self._tmm

	@property
	def gtex(self):
		if self._gtex is None: self._gtex = pd.read_csv(self._counts_data_paths["gtex"], sep="\t", index_col=0)
		return self._gtex

	@property
	def mask(self):
		if self._mask is None: 
			self._mask = self.counts.index.isin(self.gtex.index)
		return self._mask
	
	@property
	def tmm_factors(self):
		if self._tmm_factors is None: self._tmm_factors = pd.read_csv(self._counts_data_paths["tmm_factors"], sep="\t", index_col=0)["tmm"]
		return self._tmm_factors
	
	@property
	def residualizer(self):
		if self._residualizer is None: self._residualizer = residualize.Residualizer.load_covariates(path=self._covariates_path)
		return self._residualizer

	def residualize(self, df): 
		return self.residualizer.transform(df)

	#----------------------------------------------------------------------------------------------------#
	# Access QTL results
	#----------------------------------------------------------------------------------------------------#
	@property
	def qtl_sig(self):
		if self._qtl_sig is None: 
			self._qtl_sig = self._qtl_accessor.load_cis_nominal_results(self._qtl_fdr)
			self._reformat_qtls(self._qtl_sig)
		return self._qtl_sig

	@property
	def qtl_leads(self):
		if self._qtl_leads is None: 
			self._qtl_leads = self._qtl_accessor.load_cis_results(self._qtl_fdr)
			self._reformat_qtls(self._qtl_leads)
		return self._qtl_leads

	@property
	def qtl_all(self):
		if self._qtl_all is None: 
			self._qtl_all = self._qtl_accessor.load_cis_nominal_results_unfiltered()
			self._reformat_qtls(self._qtl_all)
		return self._qtl_all
	
	def _reformat_qtls(self, qtl_results): 
		"""Reformat column names to reflect feature."""
		if self._feature == "rna" or self._feature == "srna": 
			phenotype_label, distance_label = "gene_id", "distance"
		elif self._feature == "atac" or self._feature == "erna": 
			phenotype_label, distance_label = "peak_id", "distance"
		else: 
			raise ValueError

		qtl_results.rename(columns={"phenotype_id": phenotype_label, "tss_distance": distance_label}, inplace=True)
		if qtl_results.index.name == "phenotype_id": qtl_results.index.name = phenotype_label

