#!/usr/bin/env python3

from pathlib import Path

import pandas as pd

from src import base_dir, logger


# _aals_data_paths = {
# 	"metadata": base_dir / "tensorqtl_runs/harmonized_metadata.210418.txt", 
# 	"bams": base_dir / "tensorqtl_runs/harmonized_data_paths.local.210418.txt", 
# 	"vcf": base_dir / "tensorqtl_runs/genomes_210409/biallelic_known_snps.harmonized.VQSR_filtered_99.rsID.GT_only.vcf.gz", 
# 	"aals_metadata": base_dir / "data/metadata/ALS Consortium DNA Metadata 20201015 .xlsx",
# 	"gene_coords": base_dir / "tensorqtl_runs/phenotypes_210422/rna.metadata.txt.gz",
# }

_aals_data_paths = {
	"metadata": base_dir / "tensorqtl_runs/harmonized_metadata.210418.txt", 
	"bams": base_dir / "tensorqtl_runs/harmonized_data_paths.local.210418.txt", 
	"vcf": base_dir / "tensorqtl_runs/beagle5_210512/filt_anno/biallelic_known_snps.harmonized.VQSR_filtered_99.rsID.GT_only.vcf.gz", 
	"aals_metadata": base_dir / "data/metadata/ALS Consortium DNA Metadata 20201015 .xlsx",
	"gene_coords": base_dir / "tensorqtl_runs/phenotypes_210422/rna.metadata.txt.gz",
}

class AALSData: 
	"""Data associated with our samples such as sample paths and metadata."""

	def __init__(self, paths): 

		self.paths = paths

		self._metadata = None
		self._bam_paths = None
		self._ALSC_metadata = None

		self._sample_names = None

		self._peak_TF_clusters = None
		self._gene_coords = None

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
	def atac_bams(self):
		return self.bam_paths["atac_bam"]

	@property
	def rna_bams(self):
		return self.bam_paths["rna_bam"]

	@property
	def ALSC_metadata(self):
		if self._ALSC_metadata is None: 
			self._ALSC_metadata = _load_ALS_Consortium_metadata(self.paths["aals_metadata"])
		return self._ALSC_metadata

	@property
	def gene_coords(self):
		if self._gene_coords is None: 
			self._gene_coords = _load_gene_coords(self.paths["gene_coords"])
		return self._gene_coords
	

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

def _load_gene_coords(path): 
	"""Loads ENSG coords processed by `src.qtl.ProcessCounts`."""
	from src.qtl import preprocess
	regions, _ = preprocess.load_counts_metadata(path)
	return regions

#----------------------------------------------------------------------------------------------------#
aals = AALSData(_aals_data_paths)