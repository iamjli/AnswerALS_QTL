#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

from src import logger
from src.load import aals
from src.qtl import normalize
from src.utils.regions import Regions


class ProcessCounts: 
	"""
	Methods for uniform processing of omics data.
	"""
	chroms = [f"chr{i}" for i in range(1,23)] + ["chrX", "chrY"]

	def __init__(self, counts, regions, lengths, sample_frac_threshold=0.2, count_threshold=6, tpm_threshold=0.1): 
		
		self.counts = counts.copy()
		self.regions = Regions(regions.copy())
		self.lengths = lengths.copy()

		self.sample_names = self.counts.columns
		self.validate()

		# Library norm factors 
		self._cpm_factors = None
		self._tmm_factors = None

		# Counts matrices
		self._tpm = None
		self._tmm = None
		self._gtex = None

		# Thresholding parameters
		self.sample_frac_threshold = sample_frac_threshold
		self.count_threshold = count_threshold
		self.tpm_threshold = tpm_threshold
		self._mask = None

	def validate(self): 
		assert self.counts.index.isin(self.regions.index).all()
		assert self.sample_names.equals(aals.sample_names)

	def initialize(self): 
		_ = self.tpm, self.tmm, self.gtex, self.mask

	@property
	def mask(self):
		"""Filter low read counts."""
		if self._mask is None: 
			n_threshold = self.counts.shape[1] * self.sample_frac_threshold
			self._mask = (
				((self.counts >= self.count_threshold).sum(axis=1) >= n_threshold) & 
				((self.tpm >= self.tpm_threshold).sum(axis=1) >= n_threshold) & 
				(self.regions["chrom"].isin(self.chroms))
			)
		return self._mask

	#----------------------------------------------------------------------------------------------------#
	# Library size normalization factors
	#----------------------------------------------------------------------------------------------------#
	@property
	def cpm_factors(self):
		if self._cpm_factors is None:
			cpm = self.counts.sum(axis=0) / 1e6
			self._cpm_factors = cpm.rename("cpm")
		return self._cpm_factors

	@property
	def tmm_factors(self):
		if self._tmm_factors is None: 
			logger.write("Computing TMM factors...")
			# Abusing nomenclature. Technically, effective library size is TMM * CPM
			tmm_factors = normalize.edgeR_calcNormFactors(self.counts) * self.cpm_factors
			self._tmm_factors = pd.Series(data=tmm_factors, index=self.sample_names, name="tmm")
		return self._tmm_factors

	#----------------------------------------------------------------------------------------------------#
	# Normalized counts matrices
	#----------------------------------------------------------------------------------------------------#
	@property
	def tpm(self):
		if self._tpm is None: 
			self._tpm = normalize.counts_to_tpm(self.counts, self.lengths)
			rpk = self.counts.div(self.lengths / 1000, axis=0)
			self._tpm = rpk.div(rpk.sum() / 1e6, axis=1)
		return self._tpm
	
	@property
	def tmm(self):
		"""edgeR normalization: normalize by library size (counts per million) and TMM factors."""
		if self._tmm is None: 
			self._tmm = self.counts / self.tmm_factors
		return self._tmm

	@property
	def gtex(self):
		if self._gtex is None: 
			logger.write("Computing GTEX matrix...")
			self._gtex = normalize.inverse_normal_transform(self.tmm[self.mask])
		return self._gtex

	#----------------------------------------------------------------------------------------------------#
	# Save/load
	#----------------------------------------------------------------------------------------------------#
	def save_data(self, prefix, output_dir, overwrite=False): 

		metadata_path = Path(output_dir) / f"{prefix}.metadata.txt.gz"
		tmm_factors_path = Path(output_dir) / f"{prefix}.tmm_factors.txt.gz"
		counts_path = Path(output_dir) / f"{prefix}.counts.txt.gz"
		tpm_path = Path(output_dir) / f"{prefix}.tpm.txt.gz"
		tmm_path = Path(output_dir) / f"{prefix}.tmm.txt.gz"
		gtex_path = Path(output_dir) / f"{prefix}.gtex.txt.gz"

		if not overwrite: 
			assert not metadata_path.is_file() and not tmm_factors_path.is_file() and not counts_path.is_file()\
				and not tpm_path.is_file() and not tmm_path.is_file() and not gtex_path.is_file()

		combined_metadata = pd.concat({"regions": self.regions.df, "lengths": self.lengths}, axis=1)
		combined_metadata.to_csv(metadata_path, sep="\t", index=True, header=True, compression="gzip")

		self.tmm_factors.to_csv(tmm_factors_path, sep="\t", index=True, header=True, compression="gzip")
		self.counts.to_csv(counts_path, sep="\t", index=True, header=True, compression="gzip")
		self.tpm.to_csv(tpm_path, sep="\t", index=True, header=True, float_format="%.6f", compression="gzip")
		self.tmm.to_csv(tmm_path, sep="\t", index=True, header=True, float_format="%.6f", compression="gzip")
		self.gtex.to_csv(gtex_path, sep="\t", index=True, header=True, float_format="%.6f", compression="gzip")

	def set_norm_factor(self, prefix, output_dir=None): 

		tmm_factors_path = Path(output_dir) / f"{prefix}.tmm_factors.txt.gz"
		self._tmm_factors = pd.read_csv(tmm_factors_path, sep="\t", index_col=0)["tmm"]


def load_counts_metadata(path): 
	"""Load counts regions and lengths generated by `ProcessCounts`."""
	metadata = pd.read_csv(path, sep="\t", index_col=0, header=[0,1])
	regions = Regions(metadata["regions"])
	lengths = metadata["lengths"].iloc[:,0]
	return regions, lengths


#----------------------------------------------------------------------------------------------------#
# Process RNA counts data
#----------------------------------------------------------------------------------------------------#
def load_rna_counts_file(counts_path, name=None): 
	"""Loads one counts file."""
	expression = pd.read_csv(counts_path, sep='\t', skiprows=1, index_col=0)
	expression.index.rename("gene_id", inplace=True)

	# Parse metadata
	metadata = pd.DataFrame({
		"chrom": "chr"+expression["Chr"].str.split(";").str[0].astype(str), 
		"start": expression["Start"].str.split(";").str[0].astype(int), 
		"end": expression["End"].str.split(";").str[-1].astype(int), 
		"strand": expression["Strand"].str.split(';').str[0].astype(str),
		"length": expression["Length"].astype(int)
	})

	# Process counts
	counts = expression.iloc[:,5]
	if name: counts.rename(name, inplace=True)
	return metadata, counts

def get_rna_raw_counts(counts_paths, output_path):
	"""Merge raw RNA counts from series of paths."""
	if Path(output_path).is_file(): 
		return pd.read_csv(output_path, sep="\t", index_col=0, compression="gzip")
	else:
		rna_counts_dic = {} 
		for i,(guid,path) in enumerate(counts_paths.items()): 
			_,counts = load_rna_counts_file(path, guid)
			rna_counts_dic[guid] = counts
			# expression = pd.read_csv(path, sep='\t', skiprows=1, index_col=0)
			# expression.index.rename("gene_id", inplace=True)
			# rna_counts_dic[guid] = expression.iloc[:, 5].rename(guid)
			logger.update("Read {} counts files".format(i))

		# save raw counts matrix
		rna_raw_counts = pd.DataFrame.from_dict(rna_counts_dic)
		rna_raw_counts.to_csv(output_path, sep="\t", index=True, header=True, compression="gzip")
		return rna_raw_counts

#----------------------------------------------------------------------------------------------------#
# Process ATAC counts data
#----------------------------------------------------------------------------------------------------#
def load_atac_diffbind(counts_path): 

	atac = pd.read_csv(counts_path, sep="\t")
	regions, counts = atac.iloc[:,:3], atac.iloc[:,3:]

	regions.columns = ["chrom", "start", "end"]
	regions.index = regions["chrom"] + ":" + regions["start"].astype(str) + "-" + regions["end"].astype(str)
	regions.index.name = "peak_id"

	counts.index = regions.index.copy()
	return regions, counts

#----------------------------------------------------------------------------------------------------#
# Prepare phenotype files for tensorqtl and PEER correction
#----------------------------------------------------------------------------------------------------#
def write_phenotype_files(counts_data_obj, prefix, output_dir): 
	"""
	Output phenotypes in tensorqtl format. Produces 2 files: 
	 - tensorqtl phenotype input
	 - phenotype file that can be read by PEER
	"""
	regions, gtex = pd.DataFrame(counts_data_obj.regions).copy(), counts_data_obj.gtex.copy()
	regions = regions.reindex(gtex.index)

	# phenotype file for tensorqtl
	tensorqtl_phenotype_path = Path(output_dir) / f"{prefix}_gtex.bed.gz"
	if tensorqtl_phenotype_path.is_file(): 
		logger.write("Existing phenotype file found. Skipping...")
	else: 
		logger.write("Writing phenotype file for tensorqtl")
		tensorqtl_phenotype_df = _get_counts_in_tensorqtl_fomat(regions, gtex)
		tensorqtl_phenotype_df.to_csv(tensorqtl_phenotype_path, sep="\t", index=False, compression="gzip")

	# PEER does not accept strand column, so write the same file without it
	PEER_phenotype_path = Path(output_dir) / f"{prefix}_gtex.for_PEER.bed.gz"
	if PEER_phenotype_path.is_file(): 
		logger.write("Existing PEER file found. Skipping...")
	else: 
		logger.write("Writing phenotype file for PEER")
		PEER_phenotype_df = tensorqtl_phenotype_df.drop(columns=["strand"], errors="ignore")
		PEER_phenotype_df.to_csv(PEER_phenotype_path, sep="\t", index=False, compression="gzip")

def _get_counts_in_tensorqtl_fomat(phenotype_pos, phenotype_df): 
	"""
	Format dataframes
	"""
	assert phenotype_pos.index.equals(phenotype_df.index)

	# tensorqtl requires essentially a bed file prepended to the counts data
	tensorqtl_df = pd.concat([phenotype_pos, phenotype_df], axis=1)	
	tensorqtl_df["gene_id"] = phenotype_pos.index

	# reorder columns
	if "strand" in tensorqtl_df.columns: 
		metadata_cols = ["chrom", "start", "end", "gene_id", "strand"]
	else: 
		metadata_cols = ["chrom", "start", "end", "gene_id"]
	cols = metadata_cols + tensorqtl_df.columns.drop(metadata_cols).tolist()
	tensorqtl_df = tensorqtl_df[cols]

	# strip `chr`
	tensorqtl_df["chrom"] = tensorqtl_df["chrom"].str[3:]
	tensorqtl_df = tensorqtl_df.rename(columns={"chrom": "#chr"})

	return tensorqtl_df

def PEER_cmd(PEER_exec_path, phenotype_file, covariates_file, num_peer, output_prefix, output_dir):
	"""
	Command to execute PEER covariate correction. Be sure to use r-4.0.3
	""" 
	return f"time Rscript {PEER_exec_path} {phenotype_file} {output_prefix} {num_peer} -o {output_dir} --covariates {covariates_file}"


def get_covariates(metadata_path): 

	harmonized_metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)
	patient_covariates = pd.DataFrame.from_dict({
		"condition": harmonized_metadata["condition"].replace({"ALS": 1, "CTR": -1}),
		"sex": harmonized_metadata["Sex"].replace({"Female": 1, "Male": -1})
	}, orient="index")
	return patient_covariates

def write_covariates(covariates, output_path): 

	covariates.T.to_csv(output_path, sep="\t", index=True, header=True)