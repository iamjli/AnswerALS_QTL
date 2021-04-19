#!/usr/bin/env python3

import numpy as np
import pandas as pd

from src import logger, base_dir
from src import aals, hg38
from src.counts import CountsData
from src.qtl.run_tensorqtl import TensorQTLManager
from src.qtl.residualize import Residualizer




# Defaults which may be overidden

# Mainly for normalization
_default_counts_prefix = {
	"atac": "atac_pysam_200bp",
	"rna":  "rna",
	"erna": "erna_pysam_1kb",
	"srna": "",
}
_default_normalize_prefix = {
	"atac": "atac_pysam_200bp",
	"rna":  "rna",
	"erna": "rna",
	"srna": "rna",
}
_default_qtl_dirname = {
	"atac": "1_eQTL_within_1mb_210418", 
	"rna":  "2_caQTL_within_10kb_peak_200bp_diffBind_210418", 
	"erna": "", 
	"srna": "", 
}


class FeatureData: 
	"""Accessor for molecular feature-related data."""

	# Immutable tags which specify counts normalization, residualizer, prefix
	_omic_tag = {"atac": "atac", "rna": "rna", "erna": "rna", "srna": "rna"}  
	# _genomic_feature_tag = {"atac": "peak", "rna": "gene", "erna": "peak", "srna": "gene"}

	def __init__(self, feature, **kwargs): 

		self._feature = feature
		self._omic = self._omic_tag[self._feature]

		self._set_default_data_params(**kwargs)
		self._reset_properties()


	def _set_default_data_params(self, counts_prefix=None, normalize_prefix=None, qtl_dirname=None): 

		# dataset to use for raw counts themselves
		self._counts_prefix = counts_prefix if counts_prefix is not None else _default_counts_prefix[self._feature]

		# dataset to use for normalization (TMM, residual). Essentially just the source omic
		self._normlize_prefix = normalize_prefix if normalize_prefix is not None else _default_normalize_prefix[self._feature]

		# QTL
		self._qtl_dirname = qtl_dirname if qtl_dirname is not None else _default_qtl_dirname[self._feature]


	def _reset_properties(self): 

		# Data accessors
		self._counts_accessor = None
		self._qtl_accessor = None

		# self._depth_per_million = None
		self._residualizer = None


	#----------------------------------------------------------------------------------------------------#
	# Dynamically load Accessors instances.
	#----------------------------------------------------------------------------------------------------#
	@property
	def counts_accessor(self):
		if self._counts_accessor is None: 
			logger.write("Loading counts data...")
			self._counts_accessor = CountsData.load_data(prefix=self._counts_prefix)
		return self._counts_accessor

	@property
	def qtl_accessor(self):
		if self._qtl_accessor is None: 
			logger.write("Loading qtl accessor...")
			self._qtl_accessor = TensorQTLManager.load_config(dirname=self._qtl_dirname)
		return self._qtl_accessor
	

	#----------------------------------------------------------------------------------------------------#
	# Access feature metadata
	#----------------------------------------------------------------------------------------------------#
	@property
	def regions(self):
		return self.counts_accessor.regions
	

	#----------------------------------------------------------------------------------------------------#
	# Access counts data
	#----------------------------------------------------------------------------------------------------#
	@property
	def counts(self):
		return self.counts_accessor.counts

	@property
	def tmm(self):
		return self.counts_accessor.tmm

	@property
	def gtex(self):
		return self.counts_accessor.gtex

	#----------------------------------------------------------------------------------------------------#
	# Access normalizers
	#----------------------------------------------------------------------------------------------------#
	@property
	def tmm_norm_factors(self):
		return self.counts_accessor._tmm_factors
	
	@property
	def depth_per_million(self):
		return self.omic_accessor.depth_per_million
	
	@property
	def residualizer(self):
		if self._residualizer is None:
			logger.write("Loading residualizer...")
			self._residualizer = Residualizer.load_covariates(self._normlize_prefix)
		return self._residualizer








# class OmicsAccessor: 
# 	"""Accessor for RNA/ATAC files. Probably overkill but whatever."""

# 	def __init__(self, omic_tag, counts_prefix): 

# 		self.omic_tag = omic_tag
# 		self.counts_prefix = counts_prefix

# 		self._depth_per_million = None
# 		self._residualizer = None

# 	@property
# 	def depth_per_million(self): 
# 		if self._depth_per_million is None: 
# 			from src.query.bam import get_depth_per_million
# 			paths = aals.bam_paths[f"{self.omic_tag}_bam"]
# 			self._depth_per_million = get_depth_per_million(paths)
# 		return self._depth_per_million

# 	@property
# 	def residualizer(self):
# 		if self._residualizer is None: 
# 			self._residualizer = Residualizer.load_covariates(self.counts_prefix)
# 		return self._residualizer
	
	




#----------------------------------------------------------------------------------------------------#
# Load counts data
#----------------------------------------------------------------------------------------------------#






#----------------------------------------------------------------------------------------------------#
# Instantiate bam queries
#----------------------------------------------------------------------------------------------------#




#----------------------------------------------------------------------------------------------------#
# Instantiate VCF queries
#----------------------------------------------------------------------------------------------------#




