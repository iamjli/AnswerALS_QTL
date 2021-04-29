#!/usr/bin/env python3

import re
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)

import numpy as np
import pandas as pd
import pyranges as pr

from src import logger
from src.query import bam


__all__ = ["Interval", "Regions"]

tag_parser = re.compile(r"(?P<chrom>chr.{1,2}):(?P<start>\d*)-(?P<end>\d*)_?(?P<strand>[+-]?)")


#--------------------------------------------------------------------------------------------------#
class Interval: 
	"""Class for performing basic queries and manipulations on a single `Interval`."""

	def __init__(self, chrom, start, end, strand=None, name=None):

		self.chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
		self.start = int(start)
		self.end = int(end)
		self.strand = strand
		self.name = self.tag if name is None else name

		self._validate()

		self.is_stranded = self.strand == "+" or self.strand == "-"

		# position at genomic features
		self.mid = (self.start + self.end) // 2
		if self.is_stranded: 
			self.tss = self.start if self.strand == "+" else self.end
			self.tes = self.start if self.strand == "-" else self.end

	def _validate(self): 
		"""Check validity of constructor arguments."""
		assert self.end > self.start
		assert self.strand in ["+", "-", ".", None]
		# TODO: check bounds are within chromosome sizes

	@classmethod
	def load_tag(cls, tag): 
		parsed_tag = tag_parser.match(tag).groupdict()
		parsed_tag["start"], parsed_tag["end"] = int(parsed_tag["start"]), int(parsed_tag["end"])
		if parsed_tag["strand"] == "": parsed_tag["strand"] = None
		return cls(**parsed_tag)

	@classmethod
	def load_ensg(cls, gene): 
		from src.load import aals
		assert gene in aals.gene_coords.index
		chrom, start, end, strand = aals.gene_coords.loc[gene]
		return cls(chrom, start, end, strand, name=gene)

	@classmethod
	def load(cls, *args): 
		"""Lazy loading."""
		if len(args) == 1 and isinstance(args[0], Interval): 
			return args[0]
		elif len(args) == 1 and isinstance(args[0], str): 
			if args[0].startswith("chr"): 
				return cls.load_tag(args[0])
			elif args[0].startswith("ENSG"): 
				return cls.load_ensg(args[0])
			else:
				raise ValueError("Could not load Interval.")
		elif len(args) == 1 and isinstance(args[0], pd.Series): 
			return cls(**args[0], name=args[0].name)
		else: 
			return cls(*args[0])

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

	def get_pos(self, gf): 
		"""Access position by string or returns default genomic feature."""
		if gf == "ref": gf = "tss" if self.is_stranded else "mid"
		return getattr(self, f"{gf}")

	#----------------------------------------------------------------------------------------------------#
	# Operations to generate new `Interval` instances relative to `self`
	#----------------------------------------------------------------------------------------------------#
	def widen(self, w): 
		return Interval(self.chrom, self.start-w, self.end+w, self.strand)

	def slide(self, s, wrt_strand=None): 

		if self.is_stranded and s != 0 and wrt_strand is None: 
			raise ValueError("`wrt_strand` must be explicit if `Interval` is stranded.")

		if wrt_strand and self.strand == "-": s = -s 

		return Interval(self.chrom, self.start+s, self.end+s, self.strand)

	def transform(self, w=0, s=0, wrt_strand=None): 
		"""Expand the region by `window`, shift the region downstream (3' direction) by `shift`. """
		return self.widen(w=w).slide(s=s, wrt_strand=wrt_strand)

	#----------------------------------------------------------------------------------------------------#
	# Queries
	#----------------------------------------------------------------------------------------------------#
	def get_genotypes(self): 
		from src.analysis import vcf
		return vcf.query_interval(self.chrom, self.start, self.end)

	def get_rna_coverages(self): 
		coverages = self.as_Regions().get_rna_coverages()
		return coverages.iloc[0]

	def get_atac_coverages(self): 
		coverages = self.as_Regions().get_atac_coverages()
		return coverages.iloc[0]

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

	def unstrand(self): 
		if self.is_stranded:
			return Interval(self.chrom, self.start, self.end, name=self.name)
		else: 
			return Interval(self.chrom, self.start, self.end)

	def as_tuple3(self):
		return self.chrom, self.start, self.end

	def as_tuple(self):
		return self.chrom, self.start, self.end, self.strand

	def as_dict(self):
		return {"chrom": self.chrom, "start": self.start, "end": self.end, "strand": self.strand}

	def as_Regions(self):
		interval_s = pd.Series(self.dict, name=self.tag)
		return Regions(interval_s.to_frame().T)

	def __repr__(self): 
		if self.is_stranded:
			return f"{self.tag}_{self.strand}"
		else:
			return self.tag

	def length(self): 
		return self.end - self.start

# class Peak(Interval): 

# 	def __init__(self, peak_id): 

# 		parsed_tag = tag_parser.match(peak_id).groupdict()
# 		chrom, start, end = parsed_tag["chrom"], int(parsed_tag["start"]), int(parsed_tag["end"])
# 		super().__init__(chrom, start, end)

# class Gene(Interval): 

# 	def __init__(self, gene_id): 

# 		coords = aals.gene_coords.loc[gene_id]
# 		super().__init__(**coords)



#----------------------------------------------------------------------------------------------------#
# Regions subclass
#----------------------------------------------------------------------------------------------------#
class Regions(pd.DataFrame): 

	_df,_pr = None,None

	@property
	def _constructor(self): 
		# return Regions
		if "chrom" in self.columns and "start" in self.columns and "end" in self.columns: 
			return Regions
		else: 
			logger.update("Not formatted as `Regions`. Falling back to dataframe.")
			return pd.DataFrame

	@property
	def is_stranded(self):
		return "strand" in self.columns	

	@property
	def is_sorted(self):
		shifted = self.shift(fill_value=0)
		return ((self["start"] > shifted["start"]) | (self["chrom"] != shifted["chrom"])).all()

	#----------------------------------------------------------------------------------------------------#
	# Queries
	#----------------------------------------------------------------------------------------------------#
	def get_rna_coverages(self, max_size=10): 
		return bam.get_coverages_in_regions(aals.rna_bams, self)

	def get_atac_coverages(self, max_size=10): 
		return bam.get_coverages_in_regions(aals.atac_bams, self)

	#----------------------------------------------------------------------------------------------------#
	# Intersect with other regions
	#----------------------------------------------------------------------------------------------------#
	def in_interval(self, chrom, start, end): 
		return self[(self["chrom"] == "chrom") & (self["end"] > start) & (self["start"] < end)]

	def overlap_with(self, other): 
		"""Reports features that overlap with other."""
		other = _format_input_as_pyranges(other)
		overlap_idx = self.pr.overlap(other).__getattr__(self.index.name)
		return self.reindex(overlap_idx)

	def overlapping_idx(self, other, col_names=None, **kwargs): 
		"""Reports indices of overlapping intervals in self and other."""
		return _get_overlapping_regions(self, other, col_names=col_names)

	def adjacent_idx(self, hops): 
		"""Reports pairs indices of adjacent intervals. Distance is set by `hops`."""
		assert isinstance(hops, int) and hops != 0
		pos_hops = abs(hops)
		chrom_vals = self["chrom"].values
		chroms_1, chroms_2 = chrom_vals[:-pos_hops], chrom_vals[pos_hops:]
		same_chrom = chroms_1 == chroms_2
		names = self.index.name, f"{hops}_hop"
		if hops > 0: 
			return pd.DataFrame((row[:2] for row in zip(self.index[:-pos_hops], self.index[pos_hops:], same_chrom) if row[2]), columns=names)
		else: 
			return pd.DataFrame((row[:2] for row in zip(self.index[pos_hops:], self.index[:-pos_hops], same_chrom) if row[2]), columns=names)

	def k_adjacent(self, interval, k=5, gf=None, report_distance=True): 
		"""Gets the k nearest intervals in either direction."""
		interval = unpack_interval_arg(interval)
		contig_features = self[self["chrom"] == interval.chrom]
		nearest_feature = contig_features.k_nearest(interval, k=1, gf=gf, report_distance=False).index[0]
		nearest_idx = np.where(contig_features.index == nearest_feature)[0][0]

		lower_idx, upper_idx = max(nearest_idx-k, 0), min(nearest_idx+k, len(contig_features))
		regions = contig_features.iloc[lower_idx:upper_idx]

		if report_distance: regions = regions.assign(distance=self.distances_from_interval(interval, gf=gf))
		return regions

	def k_nearest(self, interval, k=10, gf=None, report_distance=True): 
		"""Gets k nearest features by absolute distance."""
		interval = unpack_interval_arg(interval)
		distances = self.distances_from_interval(interval, gf=gf)
		regions = self.reindex(distances.abs().sort_values()[:k].index).sort()
		if report_distance: regions = regions.assign(distance=self.distances_from_interval(interval, gf=gf))
		return regions

	def previous_feature(self, interval, n=1, gf=None, report_distance=True): 
		"""Gets k nearest features by absolute distance."""
		interval = unpack_interval_arg(interval)
		adjacent_intervals = self.k_adjacent(interval, k=1, gf=gf, report_distance=False)
		return Interval.load(adjacent_intervals.iloc[0])

	#----------------------------------------------------------------------------------------------------#
	# Converters and I/O
	#----------------------------------------------------------------------------------------------------#
	@property
	def annotate(self):
		# if self.omic == "rna": 
		# 	return self.assign(symbol=self["gene_id"].map(data.ensg))
		pass

	@property
	def pr(self): 
		if self._pr is None: 
			self._pr = df_to_pr(self.df)
		return self._pr

	def bed(self, strand_fill="."): 
		"""Convert `Regions` object to BED format."""
		pass

	@property
	def df(self): 
		if self._df is None: 
			self._df = pd.DataFrame(self)
		return self._df

	def write_bed(self, path):
		# TODO
		pass

	#----------------------------------------------------------------------------------------------------#
	# Utility methods
	#----------------------------------------------------------------------------------------------------#
	@property
	def tags(self): 
		return self["chrom"] + ":" + self["start"].astype(str) + "-" + self["end"].astype(str)

	def set_index_to_tags(self, name="peak_id"): 
		new_regions = self.copy()
		new_regions.index = self.tags
		new_regions.index.name = name
		return new_regions

	def sort(self): 
		return sort_regions(self)

	#----------------------------------------------------------------------------------------------------#
	# Constructors
	#----------------------------------------------------------------------------------------------------#
	def unstrand(self):
		return self.drop(columns=["strand"], errors="ignore")

	def widen(self, w): 
		"""Expand region by w."""
		new_regions = self.copy()
		new_regions["start"] -= w
		new_regions["end"] += w
		return new_regions

	def slide(self, s): 
		"""Slide region by s."""
		new_regions = self.copy()
		if self.is_stranded: 
			s = self["strand"].replace({"+":s, "-":-s})
		new_regions["start"] += s
		new_regions["end"] += s
		return new_regions

	def transform(self, w=0, s=0): 
		new_regions = self.copy()
		if self.is_stranded: 
			s = self["strand"].replace({"+":s, "-":-s})
		new_regions["start"] += s - w
		new_regions["end"] += s + w
		return new_regions
	
	#----------------------------------------------------------------------------------------------------#
	# Access positions
	#----------------------------------------------------------------------------------------------------# 
	@property
	def start(self):
		new_regions = self.copy()
		new_regions["end"] = new_regions["start"] + 1
		return new_regions

	@property
	def end(self):
		new_regions = self.copy()
		new_regions["start"] = new_regions["end"] - 1
		return new_regions

	@property
	def mid(self):
		new_regions = self.copy()
		new_regions["start"] = (new_regions["start"] + new_regions["end"]) // 2
		new_regions["end"] = new_regions["start"]
		return new_regions
	
	@property
	def tss(self):
		new_regions = self.copy()
		new_regions["start"] = np.where(new_regions["strand"] == "+", new_regions["start"], new_regions["end"]-1)
		new_regions["end"] = new_regions["start"] + 1
		return new_regions

	@property
	def tes(self):
		new_regions = self.copy()
		new_regions["start"] = np.where(new_regions["strand"] == "-", new_regions["start"], new_regions["end"]-1)
		new_regions["end"] = new_regions["start"] + 1 
		return new_regions

	@property
	def start_pos(self):
		return self["start"].copy()

	@property
	def end_pos(self):
		return self["end"].copy()

	@property
	def mid_pos(self):
		return ((self["start"] + self["end"]) // 2).rename("mid")

	@property
	def tss_pos(self):
		tss_pos = np.where(self["strand"] == "+", self["start"], self["end"])
		return pd.Series(data=tss_pos, index=self.index, name="tss")

	@property
	def tes_pos(self):
		tes_pos = np.where(self["strand"] == "-", self["start"], self["end"])
		return pd.Series(data=tes_pos, index=self.index, name="tes")

	def get_pos(self, gf=None): 
		"""Access position by string or returns default genomic feature."""
		if gf is None: gf = "tss" if self.is_stranded else "mid"
		return getattr(self, f"{gf}_pos")



	def distances_from_interval(self, interval, gf): 
		interval = unpack_interval_arg(interval)
		target_positions = self[self["chrom"] == interval.chrom].get_pos(gf=gf)
		distances = target_positions - interval.get_pos(gf=None)
		return distances*-1 if interval.strand == "-" else distances


def unpack_interval_arg(arg, regions=None): 
	if isinstance(arg, str) and arg in regions.index: 
		return Interval.load(regions.loc[arg])
	else: 
		return Interval.load(arg)

#----------------------------------------------------------------------------------------------------#
# BED operations
#----------------------------------------------------------------------------------------------------#
def _get_overlapping_regions(regions1, regions2, col_names=None, **kwargs): 

	regions1_id = _get_regions_id(regions1)
	regions2_id = _get_regions_id(regions2)

	regions1 = _format_input_as_pyranges(regions1)
	regions2 = _format_input_as_pyranges(regions2)

	joined = regions1.join(regions2, **kwargs).as_df()
	if regions1_id == regions2_id: regions2_id = regions2_id + "_b"
	pairs = joined[[regions1_id, regions2_id]]

	if col_names: pairs.columns = col_names
	return pairs


def _get_regions_id(obj): 
	"""Get id from regions."""
	if isinstance(obj, pr.PyRanges): 
		additional_cols = obj.columns.drop(["Chromosome", "Start", "End", "Strand"], errors="ignore")
		assert len(additional_cols) == 1, "cannot determine regions id"
		return additional_cols[0]
	else: 
		assert isinstance(obj, Regions), "input not formatted as pyranges or Regions object"
		return obj.index.name

def _format_input_as_pyranges(obj): 
	"""Formats pyranges or Regions object as pyranges"""
	if isinstance(obj, pr.PyRanges): 
		return obj
	else: 
		assert isinstance(obj, Regions), "input not formatted as pyranges or Regions object"
		return obj.pr
#----------------------------------------------------------------------------------------------------#
# Utility functions
#----------------------------------------------------------------------------------------------------#

def df_to_pr(df):
	"""Convert dataframe of regions to pyranges"""
	pos_columns = ["chrom", "start", "end", "strand"] if "strand" in df.columns else ["chrom", "start", "end"]
	# reorder columns to place positions first
	df = df.reset_index()
	df = df[pos_columns + df.columns.drop(pos_columns, errors="ignore").tolist()] 
	df = df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End", "strand": "Strand"})
	return pr.PyRanges(df)

def coords_to_tag(chrom, start, end): 
	return f"{chrom}:{start}-{end}"

def sort_regions(regions_df): 
	"""Lexicographical sort on bed-style regions."""
	from src.load import hg38
	
	tmp_regions = regions_df.copy()
	tmp_regions["chrom_tag"] = tmp_regions["chrom"].str[3:]
	tmp_regions["not_standard_chrom"] = ~tmp_regions["chrom"].isin(hg38.chroms)
	sorted_idx = tmp_regions.sort_values(["not_standard_chrom", "chrom_tag", "start", "end"]).index
	return regions_df.reindex(sorted_idx)

def parse_coords_index(coords): 
	regions = coords.str.extract(tag_parser)
	regions["start"] = pd.to_numeric(regions["start"], downcast="unsigned")
	regions["end"] = pd.to_numeric(regions["end"], downcast="unsigned")
	regions.index = coords 
	return regions



# #----------------------------------------------------------------------------------------------------#
# # Regions Accessor
# #----------------------------------------------------------------------------------------------------#
# @pd.api.extensions.register_dataframe_accessor("bed")
# class RegionsAccessor:

# 	def __init__(self, regions): 

# 		self._validate(regions)
# 		self._regions = regions

# 		self._pos = None
# 		self._pr = None

# 	@staticmethod
# 	def _validate(regions):
# 		assert "chrom" in regions.columns
# 		assert "start" in regions.columns
# 		assert "end" in regions.columns

# 		if "strand" in regions: 
# 			assert regions["strand"].isin(["+", "-"]).all()

# 	#----------------------------------------------------------------------------------------------------#
# 	# Regions properties
# 	#----------------------------------------------------------------------------------------------------#
# 	@property
# 	def is_stranded(self):
# 		return "strand" in self._regions.columns

# 	@property
# 	def unstrand(self):
# 		return self._regions.drop(columns=["strand"], errors="ignore")

# 	@property
# 	def pos_columns(self):
# 		return ["chrom", "start", "end", "strand"] if self.is_stranded else ["chrom", "start", "end"]

# 	@property
# 	def pos(self):
# 		"""Returns a dataframe with additional position columns."""
# 		if self._pos is None: 
# 			positions = self._regions.copy()
# 			positions["mid"] = (positions["start"] + positions["end"]) // 2
# 			if self.is_stranded: 
# 				positions["tss"] = np.where(positions["strand"] == "+", positions["start"], positions["end"])
# 				positions["tes"] = np.where(positions["strand"] == "-", positions["start"], positions["end"]) 
# 			self._pos = positions
# 		return self._pos


# 	#----------------------------------------------------------------------------------------------------#
# 	# Operations to generate new regions
# 	#----------------------------------------------------------------------------------------------------#
# 	def recenter(self, feature):
# 		"""Returns bed dataframe centered around a genomic feature."""
# 		new_regions = self._regions.copy()
# 		if feature == "start": 
# 			new_regions["end"] = new_regions["start"] + 1
# 		elif feature == "end": 
# 			new_regions["start"] = new_regions["end"] - 1
# 		elif feature == "mid": 
# 			pos = self.pos["mid"]
# 			new_regions["start"], new_regions["end"] = pos, pos + 1
# 		elif feature == "tss": 
# 			pos = self.pos["tss"] + new_regions["strand"].replace({"+":0, "-":-1})
# 			new_regions["start"], new_regions["end"] = pos, pos + 1
# 		elif feature == "tes": 
# 			pos = self.pos["tes"] + new_regions["strand"].replace({"+":0, "-":1})
# 			new_regions["start"], new_regions["end"] = pos-1, pos
# 		else:
# 			raise ValueError
# 		return new_regions

# 	def expand(self, w): 

# 		new_regions = self._regions.copy()
# 		new_regions["start"] -= w
# 		new_regions["end"] += w
# 		return new_regions

# 	def shift(self, s): 

# 		new_regions = self._regions.copy()
# 		if self.is_stranded: 
# 			s = self._regions["strand"].replace({"+":s, "-":-s})
# 		new_regions["start"] += s
# 		new_regions["end"] += s
# 		return new_regions

# 	def transform(self, w=0, s=0): 

# 		new_regions = self._regions.copy()
# 		if self.is_stranded: 
# 			s = self._regions["strand"].replace({"+":s, "-":-s})
# 		new_regions["start"] += s - w
# 		new_regions["end"] += s + w
# 		return new_regions

# 	#----------------------------------------------------------------------------------------------------#
# 	# Make queries from other data
# 	#----------------------------------------------------------------------------------------------------#
# 	def get_coverages(self, omic, max_size=10): 
# 		assert self.regions.shape[1] <= max_size
# 		if omic == "atac":
# 			return get_coverages_in_regions(aals.atac_bams, self.regions)
# 		elif omic == "rna": 
# 			return get_coverages_in_regions(aals.atac_bams, self.regions)
# 		else: 
# 			raise ValueError

# 	#----------------------------------------------------------------------------------------------------#
# 	# Utility methods
# 	#----------------------------------------------------------------------------------------------------#
# 	@property
# 	def tags(self): 
# 		return self._regions["chrom"] + ":" + self._regions["start"].astype(str) + "-" + self._regions["end"].astype(str)

# 	def set_index_to_tags(self, name="peak_id"): 
# 		new_regions = self._regions.copy()
# 		new_regions.index = self.tags
# 		new_regions.index.name = name
# 		return new_regions

# 	def sort(self): 
# 		pass
# 		# return sort_regions(self)

# 	def as_pr(self): 
# 		if self._pr is None: 
# 			self._pr = df_to_pr(self._regions)
# 		return self._pr

# 	def as_bed(self, strand_fill="."): 
# 		"""Convert `Regions` object to BED format."""
# 		pass
# 		# regions_bed = pd.DataFrame(self.reset_index())
# 		# index_name = regions_bed.columns[0]
# 		# if "score" not in regions_bed.columns: 
# 		# 	regions_bed["score"] = "."
# 		# if "strand" not in regions_bed.columns: 
# 		# 	regions_bed["strand"] = strand_fill
# 		# return regions_bed[["chrom", "start", "end", index_name, "score", "strand"]]

# 	def to_bed(self, path):
# 		# TODO
# 		pass





