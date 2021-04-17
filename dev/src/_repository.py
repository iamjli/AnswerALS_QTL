#!/usr/bin/env python3
import pandas as pd
import numpy as np
import pyranges as pr

from scipy import stats

from pathlib import Path

from . import BASE_DIR, CHROMS, logger


TENSORQTL_DIR = BASE_DIR / "tensorqtl_runs"

GENOMES_TAG 	= "genomes_210409"
PHENOTYPES_TAG 	= "phenotypes_210409"
EQTL_TAG		= "_210409_rna_gtex_maf01_PEER10_sex_condition"
CAQTL_TAG		= "_210409_atac_gtex_maf01_PEER10_sex_condition_window_10000"

RESULTS_PATHS = {
	"harmonized_metadata": 	TENSORQTL_DIR / "harmonized_metadata.210409.txt", 
	"bams": 				TENSORQTL_DIR / "harmonized_data_paths.filtered.210409.txt",
	"vcf": 					TENSORQTL_DIR / GENOMES_TAG / "biallelic_known_snps.harmonized.VQSR_filtered_99.rsID.GT_only.vcf.gz",
	"rsid": 				TENSORQTL_DIR / GENOMES_TAG / "snp_list.biallelic_known_snps.harmonized.VQSR_filtered_99.rsID.GT_only.parquet", 
	"phenotype_dir":		TENSORQTL_DIR / PHENOTYPES_TAG,
	"rna_omic_dump": 		TENSORQTL_DIR / PHENOTYPES_TAG /"_rna.counts.txt.gz", 
	# "atac_omic_dump": 		TENSORQTL_DIR / PHENOTYPES_TAG /"_omics_dump.atac_counts.txt.gz", 
	# "rna_tmm_factors": 		TENSORQTL_DIR / PHENOTYPES_TAG /"_omics_dump.rna_tmm_norm_factors.txt.gz", 
	# "atac_tmm_factors": 	TENSORQTL_DIR / PHENOTYPES_TAG /"_omics_dump.atac_tmm_norm_factors.txt.gz", 
	"rna_covariates": 		TENSORQTL_DIR / PHENOTYPES_TAG / "rna_gtex.PEER_10.obs_sex_condition.PEER_covariates.txt",
	"atac_covariates": 		TENSORQTL_DIR / PHENOTYPES_TAG / "atac_gtex.PEER_10.obs_sex_condition.PEER_covariates.txt",
	"rna_results_dir": 		TENSORQTL_DIR / EQTL_TAG, 
	"atac_results_dir": 	TENSORQTL_DIR / CAQTL_TAG, 
}

repository = 1


class Repository:

	def __init__(self, results_paths=RESULTS_PATHS): 

		self.results_paths = results_paths

		# sample metadata
		self._metadata  = None	# GUID (index): sample metadata
		self._bam_paths = None	# GUID (index): paths to bam files and counts
		self._ALSC_metadata = None

		# sample counts (phenotype data)
		self._rna_metadata = None
		# self._rna_tmm_norm_factors = None
		# self._atac_tmm_norm_factors = None
		# self._tensorqtl_rna = None
		# self._tensorqtl_atac = None

		# annotations
		self._rsid   = None		# rsID (index): positions and metadata 
		self._snp2tf = None		# rsID (index): TF and scores
		self._ENSG   = None		# ensg (index): gene symbol
		self._gencode = None

		# GWAS
		self._PM_gwas = None	# projectMinE GWAS

		# eQTL 
		self._PE_eqtl = None	# psychENCODE eQTL

		# bed files
		self._PE_enh = None		# psychENCODE enhancers
		self._chrom_sizes = None

		# gene lists
		self._OT = None			# Open targets ALS gene list

	def load_stored_vars(self, stored_obj): 
		for key,val in stored_obj.__dict__.items(): 
			if self.__dict__[key] is not None: 
				self.__dict__[key] = val

	def set_variable(self, var_name, val): 
		if self.__dict__["_"+var_name] is None: 
			logger.write("{} has already been assigned a value.".format(var_name))
		else: 
			self.__dict__["_"+var_name] = val

	@property
	def rsid(self):
		if self._rsid is None: 
			self._rsid = _load_rsid()
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
	def gencode(self):
		if self._gencode is None: self._gencode = _load_gencode_annos()
		return self._gencode

	@property
	def metadata(self):
		if self._metadata is None: self._metadata = _load_metadata()
		return self._metadata
	
	@property
	def bam_paths(self):
		if self._bam_paths is None: self._bam_paths = _load_bam_paths()
		return self._bam_paths

	@property
	def rna_metadata(self):
		if self._rna_metadata is None: 
			df = pd.read_csv(self.results_paths["rna_omic_dump"], sep="\t", index_col=0, compression="gzip")
			df = df.iloc[:,:5]
			df["tss"] = df["start"].copy()
			df.loc[df["strand"] == "-", "tss"] = df.loc[df["strand"] == "-", "end"]
			self._rna_metadata = df
		return self._rna_metadata

	@property
	def rna_tmm_norm_factors(self):
		if self._rna_tmm_norm_factors is None: 
			self._rna_tmm_norm_factors = pd.read_csv(self.results_paths["rna_tmm_factors"], sep="\t", index_col=0, compression="gzip").iloc[:,0]
		return self._rna_tmm_norm_factors

	@property
	def atac_tmm_norm_factors(self):
		if self._atac_tmm_norm_factors is None: 
			self._atac_tmm_norm_factors = pd.read_csv(self.results_paths["atac_tmm_factors"], sep="\t", index_col=0, compression="gzip").iloc[:,0]
		return self._atac_tmm_norm_factors

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
	def PE_eqtl(self):
		if self._PE_eqtl is None: 
			self._PE_eqtl = _load_psychencode()
		return self._PE_eqtl

	@property
	def chrom_sizes(self):
		if self._chrom_sizes is None: 
			self._chrom_sizes = _load_chrom_sizes()
		return self._chrom_sizes
	

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



EXTERNAL_PATHS = {
	"chrom_sizes": 		BASE_DIR / "data/external/hg38/hg38.chrom.sizes",
	"ensg": 			BASE_DIR / "data/external/ENSG_to_symbol.tsv", 
	"als_consortium": 	BASE_DIR / "data/metadata/ALS Consortium DNA Metadata 20201015 .xlsx",
	"snp2tf": 			BASE_DIR / "data/external/SNP2TF/snp2tfbs_JASPAR_CORE_2014_vert.bed.gz",
	"open_targets": 	BASE_DIR / "data/external/targets_associated_with_amyotrophic_lateral_sclerosis.csv", 
	"pe_eqtl": 			BASE_DIR / "data/external/PsychENCODE/DER-08a_hg38_eQTL.significant.txt",
	"pe_enhancers": 	BASE_DIR / "data/external/PsychENCODE/DER-04a_hg38lft_PEC_enhancers.bed", 
	"project_mine": 	BASE_DIR / "data/external/Summary_Statistics_GWAS_2016/als.sumstats.lmm.parquet",
	"project_mine_hg38": 	BASE_DIR / "data/external/Summary_Statistics_GWAS_2016/als.sumstats.lmm.hg38.parquet", 
	"gencode_gtf": 		BASE_DIR / "data/external/gencode/gencode.v34.basic.annotation.gtf", 
}

def _load_chrom_sizes(): 
	"""Loads chromosome sizes."""
	logger.write("Loading chromosome sizes...")
	chrom_sizes = pd.read_csv(EXTERNAL_PATHS["chrom_sizes"], sep="\t", names=["chrom", "end"])
	chrom_sizes = chrom_sizes[chrom_sizes["chrom"].isin(CHROMS)]
	return chrom_sizes

def _load_ENSG():
	"""Loads series of gene symbols indexed by ENSG."""
	logger.write("Loading ENSG...")
	ENSG = pd.read_csv(EXTERNAL_PATHS["ensg"], sep="\t", names=["gene_id", "symbol"], skiprows=1)
	return ENSG.set_index('gene_id')["symbol"]

def _load_metadata():
	"""Loads harmonized metadata."""
	logger.write("Loading sample metadata...")
	return pd.read_csv(RESULTS_PATHS["harmonized_metadata"], sep="\t", index_col=0)

def _load_bam_paths(): 
	"""Loads bam paths for RNA and ATAC used in tensorqtl."""
	logger.write("Loading bam paths...")
	return pd.read_csv(RESULTS_PATHS["bams"], sep="\t", index_col=0)

def _load_ALS_Consortium_metadata(): 
	"""Loads ALS Consortium metadata."""
	logger.write("Loading ALS Consortium metadata...")
	return pd.read_excel(EXTERNAL_PATHS["als_consortium"], sheet_name=0, engine="openpyxl")

def _load_rsid(): 
	"""See README.md for how this was generated. Section `Generate variant list`."""
	logger.write("Loading rsID file...")
	rsid = pd.read_parquet(RESULTS_PATHS["rsid"])
	# initialize loc - for some reason the first loc takes forever
	rsid.loc[rsid.index[0]]
	return rsid

def _load_snp2tf(collapse=True): 
	"""Load TF annotations for SNPs."""
	logger.write("Loading snp2tf annotations...")

	import gzip
	with gzip.open(EXTERNAL_PATHS["snp2tf"], "r") as f: 
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
	logger.write("Loading OpenTargets gene list...")
	als = pd.read_csv(EXTERNAL_PATHS["open_targets"])
	als.columns = als.columns.str.split(".").str[-1]
	als.set_index("symbol", inplace=True)
	return als

def _load_psychencode(): 
	logger.write("Loading psychENCODE GWAS...")
	PE_eqtls = pd.read_csv(EXTERNAL_PATHS["pe_eqtl"], sep='\t', usecols=["gene_id", "SNP_id", "nominal_pval", "regression_slope", "top_SNP"])
	PE_eqtls["SNP_id"] = "chr" + PE_eqtls["SNP_id"]
	PE_eqtls["gene_id"] = PE_eqtls["gene_id"].str.split(".").str[0]
	return PE_eqtls

def _load_psychencode_enhancers(): 
	logger.write("Loading psychENCODE enhancers...")
	return pr.read_bed(str(EXTERNAL_PATHS["pe_enhancers"]))

def _load_project_mine(): 
	logger.write("Loading Project MinE GWAS...")
	return pd.read_parquet(EXTERNAL_PATHS["project_mine_hg38"])

def _load_gencode_annos(): 
	logger.write("Loading Gencode annotations...")

	gencode_annos = pr.read_gtf(str(EXTERNAL_PATHS["gencode_gtf"]))

	gencode_annos = pr.PyRanges(pd.concat([
	    gencode_annos.features.tss().slack(1000).as_df().assign(Feature="tss"),
	    gencode_annos.features.tes().slack(1000).as_df().assign(Feature="tes"), 
	    gencode_annos[gencode_annos.Feature == "exon"].as_df(),
	    gencode_annos.features.introns().as_df()
	]))

	return gencode_annos



def load_rna_forward_counts_at_peaks(): 
	counts_path = BASE_DIR / "deeptools/rna_multicov_at_peaks.forward.readCounts.tab"

	counts_df = pd.read_csv(counts_path, sep="\t")
	counts_df.columns = counts_df.columns.str.replace("'", "")
	counts_df.rename(columns={counts_df.columns[0]: "chrom"}, inplace=True)
	counts_df.rename(columns={col:col.split("-")[1] for col in counts_df.columns[3:]}, inplace=True)
	counts_df["chrom"] = "chr"+counts_df.astype(str)

	counts_regions_df = pd.read_csv(BASE_DIR / "deeptools/atac_regions.bed", sep="\t", names=["chrom", "start", "end", "peak_id"], index_col=3)

	counts_df = counts_regions_df.merge(counts_df, on=["chrom", "start", "end"]).set_index(counts_regions_df.index)
	counts_df = counts_df.drop(columns=["chrom", "start", "end"])

	return counts_df

def load_rna_reverse_counts_at_peaks(): 
	counts_path = BASE_DIR / "deeptools/rna_multicov_at_peaks.reverse.readCounts.tab"

	counts_df = pd.read_csv(counts_path, sep="\t")
	counts_df.columns = counts_df.columns.str.replace("'", "")
	counts_df.rename(columns={counts_df.columns[0]: "chrom"}, inplace=True)
	counts_df.rename(columns={col:col.split("-")[1] for col in counts_df.columns[3:]}, inplace=True)
	counts_df["chrom"] = "chr"+counts_df.astype(str)

	counts_regions_df = pd.read_csv(BASE_DIR / "deeptools/atac_regions.bed", sep="\t", names=["chrom", "start", "end", "peak_id"], index_col=3)

	counts_df = counts_regions_df.merge(counts_df, on=["chrom", "start", "end"]).set_index(counts_regions_df.index)
	counts_df = counts_df.drop(columns=["chrom", "start", "end"])

	return counts_df






# def _load_multiBamSummary(counts_path, regions_path): 
# 	"""Loads multiBamSummary of counts."""
# 	logger.write("Loading multiBamSummary...")


# 	counts_path = BASE_DIR / "deeptools/{}.readCounts.tab".format(prefix)
# 	counts_path = BASE_DIR / "deeptools/{}.readCounts.tab".format(prefix)
