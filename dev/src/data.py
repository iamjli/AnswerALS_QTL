#!/usr/bin/env python3
import pandas as pd
import numpy as np
import pyranges as pr

from scipy import stats

from pathlib import Path

from . import BASE_DIR, CHROMS



##########################
##	   EXTERNAL DATA    ##
##########################

class External:

	def __init__(self, base_dir=BASE_DIR): 

		self.base_dir = BASE_DIR

		# sample metadata
		self._metadata  = None	# GUID (index): sample metadata
		self._bam_paths = None	# GUID (index): paths to bam files and counts
		self._ALSC_metadata = None

		# phenotype annotations
		self._rsid   = None		# rsID (index): positions and metadata 
		self._snp2tf = None		# rsID (index): TF and scores
		self._ENSG   = None		# ensg (index): gene symbol

		# GWAS metadata
		self._PE_gwas = None	# psychENCODE GWAS
		self._PM_gwas = None	# projectMinE GWAS

		# bed files
		self._PE_enh = None		# psychENCODE enhancers

		# gene lists
		self._OT = None			# Open targets ALS gene list

	@property
	def rsid(self):
		if self._rsid is None: 
			self._rsid = _load_rsIDs()
		return self._rsid

	@property
	def snp2tf(self):
		if self._snp2tf is None: 
			self._snp2tf = _load_snp2tf()
		return self._snp2tf

	@property
	def ENSG(self):
		if self._ENSG is None: self._ENSG = _load_ENSG()
		return self._ENSG

	@property
	def metadata(self):
		if self._metadata is None: self._metadata = _load_metadata()
		return self._metadata
	
	@property
	def bam_paths(self):
		if self._bam_paths is None: self._bam_paths = _load_bam_paths()
		return self._bam_paths

	@property
	def ALSC_metadata(self):
		if self._ALSC_metadata is None: 
			self._ALSC_metadata = _load_ALS_Consortium_metadata()
		return self._ALSC_metadata

	@property
	def OT(self):
		if self._OT is None: 
			self._OT = _load_opentargets()
		return self._OT

	@property
	def PE_gwas(self):
		if self._PE_gwas is None: 
			self._PE_gwas = _load_psychencode()
		return self._PE_gwas

	@property
	def PE_enh(self):
		if self._PE_enh is None: 
			self._PE_enh = _load_psychencode_enhancers()
		return self._PE_enh

	@property
	def PM_gwas(self):
		if self._PM_gwas is None: 
			self._PM_gwas = _load_project_mine()
		return self._PM_gwas



def _load_ENSG():
	"""Loads series of gene symbols indexed by ENSG."""
	path = BASE_DIR / "data/external/ENSG_to_symbol.tsv"
	ENSG = pd.read_csv(path, sep="\t", names=["gene_id", "symbol"], skiprows=1)
	return ENSG.set_index('gene_id')["symbol"]

def _load_metadata():
	"""Loads harmonized metadata."""
	path = BASE_DIR / "tensorqtl/harmonized_metadata.210313.txt" 
	return pd.read_csv(path, sep="\t", index_col=0)

def _load_bam_paths(): 
	"""Loads bam paths for RNA and ATAC used in tensorqtl."""
	path = BASE_DIR / "tensorqtl/harmonized_data_paths.210313.txt"
	return pd.read_csv(path, sep="\t", index_col=0)

def _load_ALS_Consortium_metadata(): 
	"""Loads ALS Consortium metadata."""
	path = BASE_DIR / "data/metadata/ALS Consortium DNA Metadata 20201015 .xlsx"
	return pd.read_excel(path, sheet_name=0, engine="openpyxl")

def _load_rsIDs(): 
	"""See README.md for how this was generated. Section `Generate variant list`."""
	path = BASE_DIR / "tensorqtl/genomes/snp_list_filtered.biallelic.harmonized.VQSR_filtered_99.rsID.parquet"
	return pd.read_parquet(path)

def _load_snp2tf(): 
	"""Load TF annotations for SNPs."""
	path = BASE_DIR / "data/external/SNP2TF/snp2tfbs_JASPAR_CORE_2014_vert.bed.gz"

	import gzip
	with gzip.open(path, "r") as f: 
		if collapse: 
			results = {}
			for line in f: 
				row = line.decode().rstrip("\n").split("\t")
				variant_id, tfs, scores = row[5], row[7], row[8]
				scores = ",".join([str(int(s)) for s in scores.split(",")])

				for var in variant_id.split(";"): # some rows contain two SNPs together
					if var not in results: 
						results[var] = dict(tfs=tfs, scores=scores)
					else: 
						results[var]["tfs"]    += "," + tfs
						results[var]["scores"] += "," + scores

	results_df = pd.DataFrame.from_dict(results, orient="index")
	results_df.index.name = "variant_id"
	return results_df

def _load_opentargets():
	path = BASE_DIR / "data/external/targets_associated_with_amyotrophic_lateral_sclerosis.csv"
	als = pd.read_csv(path)
	als.columns = als.columns.str.split(".").str[-1]
	als.set_index("symbol", inplace=True)
	return als

def _load_psychencode(): 
	path = BASE_DIR / "data/external/PsychENCODE/DER-08a_hg38_eQTL.significant.txt"
	PE_eqtls = pd.read_csv(path, sep='\t', usecols=["gene_id", "SNP_id", "nominal_pval", "regression_slope", "top_SNP"])
	PE_eqtls["SNP_id"] = "chr" + PE_eqtls["SNP_id"]
	PE_eqtls["gene_id"] = PE_eqtls["gene_id"].str.split(".").str[0]
	return PE_eqtls

def _load_psychencode_enhancers(): 
	path = BASE_DIR / "data/external/PsychENCODE/DER-04a_hg38lft_PEC_enhancers.bed"
	return pr.read_bed(str(path))

def _load_project_mine(): 
	path = BASE_DIR / "data/external/Summary_Statistics_GWAS_2016/als.sumstats.lmm.parquet"
	return pd.read_parquet(path)





