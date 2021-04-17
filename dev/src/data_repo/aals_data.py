#!/usr/bin/env python3

import pandas as pd

from pathlib import Path

from src import base_dir, logger


_aals_data_paths = {
	"metadata": base_dir / "tensorqtl_runs/harmonized_metadata.210409.txt", 
	"bams": base_dir / "tensorqtl_runs/harmonized_data_paths.local.filtered.210409.txt", 
}

class AALSData: 
	"""Data associated with our samples such as sample paths and metadata."""

	def __init__(self): 

		self.paths = _aals_data_paths

		self._metadata = None
		self._bam_paths = None
		self._ALSC_metadata = None

		self._sample_names = None

		assert (self.sample_names == self.bam_paths.index).all()

	@property
	def sample_names(self):
		if self._sample_names is None: 
			self._sample_names = self.metadata.index.copy()
		return self._sample_names
	
	@property
	def metadata(self):
		"""Harmonized metadata"""
		if self._metadata is None: 
			self._metadata = _load_metadata(self.paths["metadata"])
		return self._metadata
	
	@property
	def bam_paths(self):
		if self._bam_paths is None: 
			self._bam_paths = _load_bam_paths(self.paths["bams"])
		return self._bam_paths

	@property
	def ALSC_metadata(self):
		if self._ALSC_metadata is None: 
			self._ALSC_metadata = _load_bam_paths(self.paths["bams"])
		return self._ALSC_metadata


#----------------------------------------------------------------------------------------------------#
# Load data 
#----------------------------------------------------------------------------------------------------#
def _load_metadata(path):
	"""Loads harmonized metadata."""
	return pd.read_csv(path, sep="\t", index_col=0)

def _load_bam_paths(path): 
	"""Loads bam paths for RNA and ATAC used in tensorqtl."""
	return pd.read_csv(path, sep="\t", index_col=0)

def _load_ALS_Consortium_metadata(path): 
	"""Loads ALS Consortium metadata."""
	return pd.read_excel(path, sheet_name=0, engine="openpyxl")


aals = AALSData()