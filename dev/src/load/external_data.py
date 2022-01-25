#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import numpy as np

from src import base_dir, logger


# _external_data_paths = {
# 	# "rsid": base_dir / "tensorqtl_runs/genomes_210409/snp_list.biallelic_known_snps.harmonized.VQSR_filtered_99.rsID.GT_only.pkl",
# 	"rsid": base_dir / "tensorqtl_runs/genomes_210409/snp_positions.biallelic_known_snps.harmonized.VQSR_filtered_99.rsID.GT_only.pickle",
# 	"ensg": base_dir / "data/external/ENSG_to_symbol.tsv", 
# 	"snp2tf": base_dir / "data/external/SNP2TF/snp2tfbs_JASPAR_CORE_2014_vert.bed.gz",
# 	"open_targets": base_dir / "data/external/targets_associated_with_amyotrophic_lateral_sclerosis.csv", 
# 	"pe_eqtl": base_dir / "data/external/PsychENCODE/DER-08a_hg38_eQTL.significant.txt",
# 	"pe_enhancers": base_dir / "data/external/PsychENCODE/DER-04a_hg38lft_PEC_enhancers.bed", 
# 	# "project_mine": base_dir / "data/external/Summary_Statistics_GWAS_2016/als.sumstats.lmm.parquet",
# 	"project_mine": base_dir / "data/external/Summary_Statistics_GWAS_2016/processed_rsid_210409.als.sumstats.lmm.pickle",
# 	"encode_tfs": base_dir / "data/external/encode_TFs_bed/combined_peaks.bed",
# }

_external_data_paths = {
	# "rsid": base_dir / "tensorqtl_runs/genomes_210409/snp_list.biallelic_known_snps.harmonized.VQSR_filtered_99.rsID.GT_only.pkl",
	"rsid": base_dir / "tensorqtl_runs/beagle5_210512/filt_anno/snp_list.biallelic_known_snps.harmonized.VQSR_filtered_99.rsID.GT_only.pickle",
	"ensg": base_dir / "data/external/ENSG_to_symbol.tsv", 
	"snp2tf": base_dir / "data/external/SNP2TF/snp2tfbs_JASPAR_CORE_2014_vert.bed.gz",
	"open_targets": base_dir / "data/external/targets_associated_with_amyotrophic_lateral_sclerosis.csv", 
	"pe_eqtl": base_dir / "data/external/PsychENCODE/DER-08a_hg38_eQTL.significant.txt",
	"pe_enhancers": base_dir / "data/external/PsychENCODE/DER-04a_hg38lft_PEC_enhancers.bed", 
	# "project_mine": base_dir / "data/external/Summary_Statistics_GWAS_2016/als.sumstats.lmm.parquet",
	"project_mine": base_dir / "210615_analysis_phased/processed_data/processed_rsid_210615.als.sumstats.lmm.pickle",
	"encode_tfs": base_dir / "data/external/encode_TFs_bed/combined_peaks.bed",
}

class ExternalData: 
	"""Externally downloaded data that has been preprocessed."""

	def __init__(self, paths): 

		self.paths = paths

		self._rsid = None
		self._ensg = None
		self._snp2tf = None
		self._open_targets = None
		self._PE_eqtl = None
		self._PE_enh = None
		self._PM_gwas = None
		self._encode_tfs = None
	
	@property
	def rsid(self):
		if self._rsid is None: 
			self._rsid = _load_rsid_pickle(self.paths["rsid"])
			# self._rsid = _load_rsid_parquet(self.paths["rsid"])
		return self._rsid

	@property
	def ensg(self):
		if self._ensg is None: 
			self._ensg = _load_ensg(self.paths["ensg"])
		return self._ensg

	@property
	def snp2tf(self):
		if self._snp2tf is None: 
			self._snp2tf = _load_snp2tf(self.paths["snp2tf"])
		return self._snp2tf

	@property
	def open_targets(self): 
		if self._open_targets is None: 
			self._open_targets = _load_opentargets(self.paths["open_targets"])
		return self._open_targets

	@property
	def PE_eqtl(self):
		if self._PE_eqtl is None: 
			self._PE_eqtl = _load_psychencode(self.paths["pe_eqtl"])
		return self._PE_eqtl

	@property
	def PE_enh(self):
		if self._PE_enh is None: 
			self._PE_enh = _load_psychencode_enhancers(self.paths["pe_enhancers"])
		return self._PE_enh

	@property
	def PM_gwas(self):
		if self._PM_gwas is None: 
			self._PM_gwas = _load_project_mine(self.paths["project_mine"])
		return self._PM_gwas

	@property
	def encode_tfs(self):
		if self._encode_tfs is None: 
			self._encode_tfs = _load_encode_tfs(self.paths["encode_tfs"])
		return self._encode_tfs


#----------------------------------------------------------------------------------------------------#
# Load data 
#----------------------------------------------------------------------------------------------------#
def _load_rsid_pickle(path): 
	"""See README.md for how this was generated. Section `Generate variant list`."""
	logger.write("Loading rsID file...")
	rsid = pd.read_pickle(path)
	rsid.loc[rsid.index[0]] # initialize loc - for some reason the first loc takes forever
	return rsid

def _load_rsid_parquet(path): 
	"""See README.md for how this was generated. Section `Generate variant list`."""
	logger.write("Loading rsID file...")
	rsid = pd.read_parquet(path)
	rsid.loc[rsid.index[0]] # initialize loc - for some reason the first loc takes forever
	return rsid

def _load_ensg(path):
	"""Loads series of gene symbols indexed by ENSG."""
	return pd.read_csv(path, sep="\t", names=["gene_id", "symbol"], skiprows=1, index_col=0)["symbol"]

def _load_snp2tf(path, collapse=True): 
	"""Load TF annotations for SNPs."""

	# df = pd.read_csv("../data/external/SNP2TF/snp2tfbs_JASPAR_CORE_2014_vert.bed.gz", sep='\t', names=['chrom', 'start', 'end', 'ref', 'alt', 'variant_id', 'matches', 'tfs', 'scores'])

	# df['variant_id'] = df['variant_id'].str.split(';')
	# df = explode(df, 'variant_id')

	# df['tfs'] = df['tfs'].str.split(',')
	# df['scores'] = df['scores'].str.split(',')
	# df = explode(df, ['tfs', 'scores'])
	# df['scores'] = df['scores'].astype(int)

	# df = df.merge(rsid, left_on='variant_id', right_index=True)

	# df = pd.concat([
	#     df.query('ref_x == ref_y & alt_x == alt_y')[['variant_id', 'tfs', 'scores']], 
	#     df.query('ref_x == alt_y & alt_x == ref_y')[['variant_id', 'tfs', 'scores']].assign(scores=lambda df: df['scores'] * -1)
	# ])
	
	import gzip

	logger.write("Loading snp2tf annotations...")

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
						results[var]["tfs"]	+= "," + tfs
						results[var]["scores"] += "," + scores

	results_df = pd.DataFrame.from_dict(results, orient="index")
	results_df.index.name = "variant_id"
	return results_df

def _load_opentargets(path):
	als = pd.read_csv(path)
	als.columns = als.columns.str.split(".").str[-1]
	als.set_index("symbol", inplace=True)
	return als

def _load_psychencode(path): 
	logger.write("Loading psychENCODE eQTLs...")
	PE_eqtls = pd.read_csv(path, sep='\t', usecols=["gene_id", "SNP_id", "nominal_pval", "regression_slope", "top_SNP"])
	PE_eqtls["SNP_id"] = "chr" + PE_eqtls["SNP_id"]
	PE_eqtls["gene_id"] = PE_eqtls["gene_id"].str.split(".").str[0]
	return PE_eqtls

def _load_psychencode_enhancers(path): 
	logger.write("Loading psychENCODE enhancers...")
	import pyranges as pr
	return pr.read_bed(str(path))

def _load_project_mine(path): 
	logger.write("Loading Project MinE GWAS...")
	return pd.read_pickle(path)

def _load_encode_tfs(path): 
	logger.write("Loading ENCODE merged TFs...")
	import pyranges as pr
	return pr.read_bed(str(path))



def _process_PM_gwas(path, rsid): 

	gwas = pd.read_parquet(path)
	gwas = gwas[gwas.index.isin(rsid.index)]

	_Ncases = 12577
	_Ncontrols = 23475
	_Pcase = _Ncases / (_Ncases + _Ncontrols)

	gwas["b"] *= -1  # slope now refers to a2. If positive, a2 occurs more frequently in ALS compared to CTR
	gwas.loc[gwas.a1 != gwas.index.map(rsid["ref"]), "b"] *= -1  # use alt allele from rsID. If mismatch, change slope
	gwas["LOD"] = np.log(((_Pcase + gwas.b) / (1 - _Pcase - gwas.b)) / (_Pcase / (1 - _Pcase)))

	gwas["a1"] = gwas["a1"].astype("category")
	gwas["a2"] = gwas["a2"].astype("category")

	return gwas.drop(columns=["chr", "bp"])

#----------------------------------------------------------------------------------------------------#
data = ExternalData(_external_data_paths)