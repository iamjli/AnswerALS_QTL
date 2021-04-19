#!/usr/bin/env python3

import pandas as pd

from src import base_dir


def chrom_list(): 
	return [f"chr{i}" for i in range(1,23)] + ["chrX", "chrY"]

def lexicographic_map(sorted_list): 
	return {key:i for i,key in enumerate(sorted_list)}


_genomes_data_paths = {
	"chrom_sizes": base_dir / "data/external/hg38/hg38.chrom.sizes"
}

class ReferenceGenomeHG38: 
	"""Data associated with hg38 genome build."""

	def __init__(self, build): 

		assert build == "hg38"
		self.build = build

		self._chroms = chrom_list()
		self._chrom_lex_map = lexicographic_map(self.chroms)

		# TODO: chrom_sizes

	@property
	def chroms(self):
		return self._chroms

	@property
	def chrom_lex_map(self):
		return self._chrom_lex_map

	@property
	def chrom_sizes(self):
		return _load_chrom_sizes(_genomes_data_paths["chrom_sizes"])

	

hg38 = ReferenceGenomeHG38("hg38")




def _load_chrom_sizes(path): 
	chrom_sizes = pd.read_csv(path, sep="\t", names=["chrom", "end"])
	chrom_sizes["start"] = 0
	chrom_sizes = chrom_sizes[["chrom", "start", "end"]]
	chrom_sizes = chrom_sizes[chrom_sizes["chrom"].isin(chrom_list())]
	chrom_sizes.index = chrom_sizes["chrom"] + ":" + chrom_sizes["start"].astype(str) + "-" + chrom_sizes["end"].astype(str)
	return chrom_sizes