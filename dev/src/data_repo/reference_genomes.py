#!/usr/bin/env python3

import pandas as pd


def chrom_list(): 
	return [f"chr{i}" for i in range(1,23)] + ["chrX", "chrY"]

def lexicographic_map(sorted_list): 
	return {key:i for i,key in enumerate(sorted_list)}


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
	

hg38 = ReferenceGenomeHG38("hg38")