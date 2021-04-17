#!/usr/bin/env python3

import numpy as np
import pandas as pd

from src import logger, hg38


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

	#----------------------------------------------------------------------------------------------------#
	# Output formats
	#----------------------------------------------------------------------------------------------------#
	@property
	def tag(self):
		return coords_to_tag(self.chrom, self.start, self.end)

	@property
	def tuple(self):
		return self.chrom, self.start, self.end, self.strand

	def __repr__(self): 
		if self.is_stranded:
			return f"{self.tag}_{self.strand}"
		else:
			return self.tag


#--------------------------------------------------------------------------------------------------#
class Regions(pd.DataFrame): 
	"""
	Multiple genomic intervals stored as a dataframe.
	"""
	# TODO: need to be careful with copies. Implement properties that generate copies of each column.
	@property
	def _constructor(self):
		self._validate_inputs()
		return Regions
	
	def _validate_inputs(self): 
		"""Check validity of constructor arguments."""
		assert "chrom" in self.columns
		assert "start" in self.columns
		assert "end" in self.columns

		assert (self["end"] > self["start"]).all()

		if self.is_stranded: 
			assert self["strand"].isin(["+", "-"]).all()

		# TODO: check bounds are within chromosome sizes

	@property
	def is_stranded(self):
		return "strand" in self.columns
	
	def sort(self): 
		return sort_regions(self)

	def as_bed(self, strand_fill="."): 
		"""Convert `Regions` object to BED format."""
		regions_bed = pd.DataFrame(self.reset_index())
		index_name = regions_bed.columns[0]
		if "score" not in regions_bed.columns: 
			regions_bed["score"] = "."
		if "strand" not in regions_bed.columns: 
			regions_bed["strand"] = strand_fill
		return regions_bed[["chrom", "start", "end", index_name, "score", "strand"]]

	def to_bed(self, path):
		# TODO
		pass

	#----------------------------------------------------------------------------------------------------#
	# Access positions
	#----------------------------------------------------------------------------------------------------#
	@property
	def mid(self):
		return (self["start"] + self["end"]) // 2

	@property
	def tss(self):
		assert self.is_stranded
		return pd.Series(data=np.where(self["strand"] == "+", self["start"], self["end"]), index=self.index)

	@property
	def tss(self):
		assert self.is_stranded
		return pd.Series(data=np.where(self["strand"] == "-", self["start"], self["end"]), index=self.index)

	@property
	def as_start(self):
		df = self.copy()
		df["end"] = df["start"] + 1
		return Regions(df)

	@property
	def as_end(self):
		df = self.copy()
		df["start"] = df["end"] - 1
		return Regions(df)

	@property
	def as_mid(self):
		df, mid = self.copy(), self.mid
		df["start"], df["end"] = self.mid, self.mid + 1
		return Regions(df)

	@property
	def as_tss(self):
		assert self.is_stranded
		df, tss = self.copy(), self.tss
		df["start"], df["end"] = self.tss, self.tss + 1
		return Regions(df)
	
	@property
	def as_tes(self):
		assert self.is_stranded
		df, tes = self.copy(), self.tes
		df["start"], df["end"] = self.tes, self.tes + 1
		return Regions(df)

	#----------------------------------------------------------------------------------------------------#
	# Operations to generate new `Interval` instances relative to `self`
	#----------------------------------------------------------------------------------------------------#
	def expand(self, window): 
		df = self.copy()
		df["start"] -= window
		df["end"] += window
		return df

	def shift(self, shift, wrt_strand=None): 

		if self.is_stranded and shift != 0 and wrt_strand is None: 
			raise ValueError("`wrt_strand` must be explicit if `Regions` is stranded.")

		if wrt_strand:
			assert self.is_stranded
			shift = self["strand"].replace({"+":shift, "-":-shift})

		new_regions = self.copy()
		new_regions["start"] += shift
		new_regions["end"] += shift
		return Regions(new_regions)

	def transform(self, window=0, shift=0, wrt_strand=None): 
		"""Expand the region by `window`, shift the region downstream (3' direction) by `shift`. """
		return self.expand(window=window).shift(shift=shift, wrt_strand=wrt_strand)

	# def get_feature_regions(self, feature): 
	# 	"""
	# 	Gets 1bp length `Regions` with the start coordinate set to the feature of interest.
	# 	"""
	# 	feature_pos = self.get_feature_pos(feature, adjust_end=True)
	# 	if self.is_stranded: 
	# 		feature_regions = pd.DataFrame.from_dict({"chrom": self.chrom, "start": feature_pos, 
	# 			"end": feature_pos+1, "strand": self.strand})
	# 	else: 
	# 		feature_regions = pd.DataFrame.from_dict({"chrom": self.chrom, "start": feature_pos, 
	# 			"end": feature_pos+1})

	# 	return Regions(feature_regions)

