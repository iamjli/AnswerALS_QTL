#!/usr/bin/env python3

import pandas as pd
import pyranges as pr

from src import base_dir, logger


_genomes_data_paths = {
	"chrom_sizes": base_dir / "data/external/hg38/hg38.chrom.sizes", 
	"gencode_gtf": base_dir / "data/external/gencode/gencode.v34.basic.annotation.gtf", 
}

class ReferenceGenome: 
	"""Data associated with hg38 genome build."""

	def __init__(self, paths): 

		self.paths = paths

		self._chroms = chrom_list()
		self._chrom_lex_map = lexicographic_map(self.chroms)

		self._gencode_gtf = None
		self._gencode_annos = None

	@property
	def chroms(self):
		return self._chroms

	@property
	def chrom_lex_map(self):
		return self._chrom_lex_map

	@property
	def chrom_sizes(self):
		return _load_chrom_sizes(self.paths["chrom_sizes"])

	@property
	def gencode_gtf(self):
		if self._gencode_gtf is None: 
			self._gencode_gtf = _load_gencode_gtf(self.paths["gencode_gtf"])
		return self._gencode_gtf

	@property
	def gencode_annos(self):
		if self._gencode_annos is None: 
			self._gencode_annos = _parse_gencode_annos(self.gencode_gtf, tss_slack=1000, tes_slack=1000)
		return self._gencode_annos


#----------------------------------------------------------------------------------------------------#
# Utility 
#----------------------------------------------------------------------------------------------------#
def chrom_list(): 
	return [f"chr{i}" for i in range(1,23)] + ["chrX", "chrY"]

def lexicographic_map(sorted_list): 
	return {key:i for i,key in enumerate(sorted_list)}

#----------------------------------------------------------------------------------------------------#
# Load data 
#----------------------------------------------------------------------------------------------------#
def _load_chrom_sizes(path): 
	chrom_sizes = pd.read_csv(path, sep="\t", names=["chrom", "end"])
	chrom_sizes["start"] = 0
	chrom_sizes = chrom_sizes[["chrom", "start", "end"]]
	chrom_sizes = chrom_sizes[chrom_sizes["chrom"].isin(chrom_list())]
	chrom_sizes.index = chrom_sizes["chrom"] + ":" + chrom_sizes["start"].astype(str) + "-" + chrom_sizes["end"].astype(str)
	return chrom_sizes

def _load_gencode_gtf(path): 
	logger.write("Loading Gencode annotations...")
	return pr.read_gtf(path)

def _parse_gencode_annos(gencode_gtf, tss_slack=1000, tes_slack=1000): 
	return pr.PyRanges(pd.concat([
		gencode_gtf.features.tss().slack(tss_slack).as_df().assign(Feature="tss"),
		gencode_gtf.features.tes().slack(tes_slack).as_df().assign(Feature="tes"), 
		gencode_gtf[gencode_gtf.Feature == "exon"].as_df(),
		gencode_gtf.features.introns(by="transcript").as_df()
	]))

#----------------------------------------------------------------------------------------------------#
hg38 = ReferenceGenome(_genomes_data_paths)