#!/usr/bin/env python3

import pandas as pd
import pyranges as pr

from src import base_dir, logger


_genomes_data_paths = {
	"chrom_sizes": base_dir / "data/external/hg38/hg38.chrom.sizes", 
	"gencode_gtf": base_dir / "data/external/gencode/gencode.v38.annotation.gtf", 
	# "gencode_gtf": base_dir / "data/external/gencode/gencode.v34.basic.annotation.gtf", 
}

class ReferenceGenome: 
	"""Data associated with hg38 genome build."""

	def __init__(self, paths): 

		self.paths = paths

		self._chroms = chrom_list()
		self._chrom_lex_map = lexicographic_map(self.chroms)

		self._gencode_gtf = None
		self._gencode_annos = None

		self._introns, self._exons, self._splice_junctions = None, None, None

	@property
	def chroms(self):
		return self._chroms

	@property
	def chrom_lex_map(self):
		return self._chrom_lex_map

	@property
	def chrom_sizes(self):
		return _load_chrom_sizes(self.paths["chrom_sizes"])

	def get_binned_chroms(self, size): 
		chrom_sizes = pr.PyRanges(hg38.chrom_sizes.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
		bins = chrom_sizes.window(size).as_df()
		bins.columns = ["chrom", "start", "end"]
		bins.index = bins["chrom"].astype(str) + ":" + bins["start"].astype(str) + "-" + bins["end"].astype(str)
		bins.index.name = "bin_id"
		return bins


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

	@property
	def exons(self):
		if self._exons is None: 
			self._exons = self.gencode_gtf[self.gencode_gtf.Feature == "exon"]
		return self._exons

	@property
	def introns(self):
		if self._introns is None:
			intron_by_transcript = self.gencode_gtf.features.introns(by='transcript').as_df().assign(Feature='intron-transcript')
			intron_by_gene = self.gencode_gtf.features.introns(by='gene').as_df().assign(Feature='intron-gene')
			# avoid duplicates from transcript-defined introns
			intron_by_gene = intron_by_gene.merge(intron_by_transcript[['Chromosome', 'Start', 'End', 'Strand', 'gene_id']], how='left', indicator='source').query('source == "left_only"').drop(columns=['source'])
			self._introns = pd.concat([intron_by_transcript, intron_by_gene])

			# self._introns = pr.PyRanges(pd.concat([
			# 	hg38.gencode_gtf.features.introns(by="transcript").as_df(), 
			# 	hg38.gencode_gtf.features.introns(by="gene").as_df()
			# ]))
		return self._introns

	@property
	def splice_junctions(self):
		if self._splice_junctions is None: 
			self._splice_junctions = pr.PyRanges(pd.concat([
				self.introns.five_end().as_df(), 
				self.introns.three_end().as_df(), 
			]))
		return self._splice_junctions




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
	gtf = pr.read_gtf(path)
	gtf = gtf.assign("gene_id", lambda df: df.gene_id.str.split(".").str[0])
	return gtf

def _parse_gencode_annos(gencode_gtf, tss_slack=1000, tes_slack=1000): 
	df = pd.concat([
		gencode_gtf.features.tss().slack(tss_slack).as_df().assign(Feature="tss"),
		gencode_gtf.features.tes().slack(tes_slack).as_df().assign(Feature="tes"), 
		gencode_gtf[gencode_gtf.Feature == "exon"].as_df(),
		gencode_gtf.features.introns(by="transcript").as_df()
	])
	df["gene_id"] = df["gene_id"].str.split(".").str[0]
	return pr.PyRanges(df)

#----------------------------------------------------------------------------------------------------#
hg38 = ReferenceGenome(_genomes_data_paths)