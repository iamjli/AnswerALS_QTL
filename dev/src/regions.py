#!/usr/bin/env python3
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)

import numpy as np
import pandas as pd
import pyranges as pr

from src import logger
from src import aals, hg38
from src.query.bam import get_pileups_in_interval, get_coverages_in_regions


def coords_to_tag(chrom, start, end): 
	return f"{chrom}:{start}-{end}"

def sort_regions(regions_df): 
	"""Lexicographical sort on bed-style regions."""
	return regions_df.sort_values(by=["start", "end"])\
		.sort_values(by=["chrom"], key=lambda col: col.map(hg38.chrom_lex_map))


#--------------------------------------------------------------------------------------------------#
class Interval: 
	"""Class for performing basic queries and manipulations on a single `Interval`."""

	def __init__(self, chrom, start, end, strand=None):

		self.chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
		self.start = int(start)
		self.end = int(end)
		self.strand = strand

		self._validate_inputs()

		self.is_stranded = self.strand == "+" or self.strand == "-"

		# position at genomic features
		self.mid = (self.start + self.end) // 2
		if self.is_stranded: 
			self.tss = self.start if self.strand == "+" else self.end
			self.tes = self.start if self.strand == "-" else self.end

	def _validate_inputs(self): 
		"""Check validity of constructor arguments."""
		assert self.end > self.start
		assert self.strand in ["+", "-", ".", None]
		# TODO: check bounds are within chromosome sizes

	#----------------------------------------------------------------------------------------------------#
	# Access genomic features as `Interval` objects
	#----------------------------------------------------------------------------------------------------#
	@property
	def as_start(self):
		return Interval(self.chrom, self.start, self.start+1, self.strand)

	@property
	def as_end(self):
		return Interval(self.chrom, self.end-1, self.end, self.strand)

	@property
	def as_mid(self):
		return Interval(self.chrom, self.mid, self.mid+1, self.strand)

	@property
	def as_tss(self):
		return Interval(self.chrom, self.tss, self.tss+1, self.strand)

	@property
	def as_tes(self):
		return Interval(self.chrom, self.tes-1, self.tes, self.strand)

	#----------------------------------------------------------------------------------------------------#
	# Operations to generate new `Interval` instances relative to `self`
	#----------------------------------------------------------------------------------------------------#
	def expand(self, window): 
		return Interval(self.chrom, self.start-window, self.end+window, self.strand)

	def shift(self, shift, wrt_strand=None): 

		if self.is_stranded and shift != 0 and wrt_strand is None: 
			raise ValueError("`wrt_strand` must be explicit if `Interval` is stranded.")

		if wrt_strand and self.strand == "-": shift = -shift 

		return Interval(self.chrom, self.start+shift, self.end+shift, self.strand)

	def transform(self, window=0, shift=0, wrt_strand=None): 
		"""Expand the region by `window`, shift the region downstream (3' direction) by `shift`. """
		return self.expand(window=window).shift(shift=shift, wrt_strand=wrt_strand)

	def get_pileups(self): 
		# TODO
		# get_pileups_in_interval
		pass

	#----------------------------------------------------------------------------------------------------#
	# Output formats
	#----------------------------------------------------------------------------------------------------#
	@property
	def tag(self):
		return coords_to_tag(self.chrom, self.start, self.end)

	@property
	def tuple3(self):
		return self.chrom, self.start, self.end

	@property
	def tuple(self):
		return self.chrom, self.start, self.end, self.strand

	def __repr__(self): 
		if self.is_stranded:
			return f"{self.tag}_{self.strand}"
		else:
			return self.tag




#----------------------------------------------------------------------------------------------------#
# Regions Accessor
#----------------------------------------------------------------------------------------------------#
@pd.api.extensions.register_dataframe_accessor("bed")
class RegionsAccessor:

	def __init__(self, regions): 

		self._validate(regions)
		self._regions = regions

		self._pos = None
		self._pr = None

	@staticmethod
	def _validate(regions):
		assert "chrom" in regions.columns
		assert "start" in regions.columns
		assert "end" in regions.columns

		if "strand" in regions: 
			assert regions["strand"].isin(["+", "-"]).all()

	#----------------------------------------------------------------------------------------------------#
	# Regions properties
	#----------------------------------------------------------------------------------------------------#
	@property
	def is_stranded(self):
		return "strand" in self._regions.columns

	@property
	def unstrand(self):
		return self._regions.drop(columns=["strand"], errors="ignore")

	@property
	def pos_columns(self):
		return ["chrom", "start", "end", "strand"] if self.is_stranded else ["chrom", "start", "end"]

	@property
	def pos(self):
		"""Returns a dataframe with additional position columns."""
		if self._pos is None: 
			positions = self._regions.copy()
			positions["mid"] = (positions["start"] + positions["end"]) // 2
			if self.is_stranded: 
				positions["tss"] = np.where(positions["strand"] == "+", positions["start"], positions["end"])
				positions["tes"] = np.where(positions["strand"] == "-", positions["start"], positions["end"]) 
			self._pos = positions
		return self._pos


	#----------------------------------------------------------------------------------------------------#
	# Operations to generate new regions
	#----------------------------------------------------------------------------------------------------#
	def recenter(self, feature):
		"""Returns bed dataframe centered around a genomic feature."""
		new_regions = self._regions.copy()
		if feature == "start": 
			new_regions["end"] = new_regions["start"] + 1
		elif feature == "end": 
			new_regions["start"] = new_regions["end"] - 1
		elif feature == "mid": 
			pos = self.pos["mid"]
			new_regions["start"], new_regions["end"] = pos, pos + 1
		elif feature == "tss": 
			pos = self.pos["tss"] + new_regions["strand"].replace({"+":0, "-":-1})
			new_regions["start"], new_regions["end"] = pos, pos + 1
		elif feature == "tes": 
			pos = self.pos["tes"] + new_regions["strand"].replace({"+":0, "-":1})
			new_regions["start"], new_regions["end"] = pos-1, pos
		else:
			raise ValueError
		return new_regions

	def expand(self, w): 

		new_regions = self._regions.copy()
		new_regions["start"] -= w
		new_regions["end"] += w
		return new_regions

	def shift(self, s): 

		new_regions = self._regions.copy()

		if self.is_stranded: 
			s = self._regions["strand"].replace({"+":s, "-":-s})
		new_regions["start"] += s
		new_regions["end"] += s
		return new_regions

	#----------------------------------------------------------------------------------------------------#
	# Make queries from other data
	#----------------------------------------------------------------------------------------------------#
	def get_coverages(self, omic, max_size=10): 
		assert self.regions.shape[1] <= max_size
		if omic == "atac":
			return get_coverages_in_regions(aals.atac_bams, self.regions)
		elif omic == "rna": 
			return get_coverages_in_regions(aals.atac_bams, self.regions)
		else: 
			raise ValueError

	#----------------------------------------------------------------------------------------------------#
	# Utility methods
	#----------------------------------------------------------------------------------------------------#
	@property
	def tags(self): 
		return self._regions["chrom"] + ":" + self._regions["start"].astype(str) + "-" + self._regions["end"].astype(str)

	def set_index_to_tags(self, name="peak_id"): 
		new_regions = self._regions.copy()
		new_regions.index = self.tags
		new_regions.index.name = name
		return new_regions

	def sort(self): 
		pass
		# return sort_regions(self)

	def as_pr(self): 
		if self._pr is None: 
			self._pr = df_to_pr(self._regions)
		return self._pr

	def as_bed(self, strand_fill="."): 
		"""Convert `Regions` object to BED format."""
		pass
		# regions_bed = pd.DataFrame(self.reset_index())
		# index_name = regions_bed.columns[0]
		# if "score" not in regions_bed.columns: 
		# 	regions_bed["score"] = "."
		# if "strand" not in regions_bed.columns: 
		# 	regions_bed["strand"] = strand_fill
		# return regions_bed[["chrom", "start", "end", index_name, "score", "strand"]]

	def to_bed(self, path):
		# TODO
		pass




def df_to_pr(df):
	"""Convert dataframe of regions to pyranges"""
	pos_columns = ["chrom", "start", "end", "strand"] if "strand" in df.columns else ["chrom", "start", "end"]
	# reorder columns to place positions first
	df = df.reset_index()
	df = df[pos_columns + df.columns.drop(pos_columns, errors="ignore").tolist()] 
	df = df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End", "strand": "Strand"})
	return pr.PyRanges(df)
