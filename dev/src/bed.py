#!/usr/bin/env python3

import re

import pandas as pd
import numpy as np
import pyranges as pr

from . import DATA, logger


class Regions(pd.DataFrame): 

	# temporary properties
	_internal_names = pd.DataFrame._internal_names + ["_omic", "_phen_id", "_pos_cols"]
	_internal_names_set = set(_internal_names)

	_pr,_tss_pos,_tes_pos,_mid_pos,_ref = None, None, None, None, None

	@property
	def _constructor(self): 

		if self.index.name == "gene_id": 
			# Initialize as RNA-associated 
			self._omic = "rna"
			self._phen_id = "gene_id"
			self._pos_cols = ["chrom", "start", "end", "strand"]
		elif self.index.name == "peak_id": 
			# Initialize as ATAC-associated 
			self._omic = "atac"
			self._phen_id = "peak_id"
			self._pos_cols = ["chrom", "start", "end"]
		else: 
			pass
		return Regions

	@property
	def omic(self):
		if self.index.name == "gene_id": 
			return "rna"
		elif self.index.name == "peak_id": 
			return "atac"
		else: 
			pass	

	@property
	def annotate(self):
		if self.omic == "rna": 
			return self.assign(symbol=self["gene_id"].map(DATA.ENSG))

	@property
	def pr(self):
		"""Returns Regions dataframe as pyranges object."""
		if self._pr is None: 
			df = self.copy().reset_index()
			df = df[ self._pos_cols + df.columns.drop(self._pos_cols, errors="ignore").tolist() ]  # reorder columns
			df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End", "strand": "Strand"}, inplace=True)
			self._pr = pr.PyRanges(df)
		return self._pr

	def stdout(self):
		df = self.reset_index()[["chrom", "start", "end", self._phen_id]]
		return df.to_string(index=False, header=False).lstrip().replace("\n ", "\\n").replace("  ", "\\t")	

	@property
	def sort(self):
		return self.sort_values(["chrom", "start"])
	
	@property
	def mid_pos(self):
		"""Returns midpoints as series."""
		if self._mid_pos is None: 
			self._mid_pos = self[["start", "end"]].mean(axis=1).astype(int)
		return self._mid_pos

	@property
	def tss_pos(self):
		if self._tss_pos is None: 
			assert "strand" in self.columns
			tss_pos = self.start.copy()
			tss_pos.loc[self["strand"] == "-"] = self["end"]
			self._tss_pos = tss_pos.astype(int)
		return self._tss_pos

	@property
	def ref_pos(self): 
		if self.omic == "rna": 
			return self.tss_pos
		else: 
			return self.mid_pos

	@property
	def ref(self):
		if self._ref is None: 
			df = self.copy()
			df["start"] = self.ref_pos
			df["end"] = self.ref_pos+1
			self._ref = Regions(df)
		return self._ref
	
	def get_features_in_region(self, chrom, start, end): 
		"""Gets Regions within a specified region."""
		region = pr.from_dict({"Chromosome": [chrom], "Start": [start], "End": [end]})
		overlaps = self.pr.overlap(region)
		# return overlaps
		return pr_to_df(overlaps)
	
	def get_features_in_window(self, phen_id, w): 
		"""Get Regions within window of specified phenotype."""
		region = self.loc[phen_id]
		return self.get_features_in_window_by_coord(**region, w=w)

	def get_features_in_window_by_coord(self, chrom, start, end, strand=None, w=0):

		query_pr = get_coords_as_pr(chrom, start, end, strand)
		query_pos = get_reference_position(query_pr).as_df().iloc[0]["Start"]

		overlaps = pr_to_df(self.pr.overlap(query_pr.slack(w)))
		overlaps_pos = overlaps.ref_pos
		overlaps["distance"] = overlaps_pos - query_pos
		if strand == "+": 
			overlaps["rel_distance"] = overlaps["distance"]
		elif strand == "-": 
			overlaps["rel_distance"] = overlaps["distance"] * -1
		else: 
			pass
		return overlaps

	def get_k_nearest_features_by_coord(self, chrom, start, end, strand=None, k=10): 
		# get reference positions of query and regions
		query_pos_pr = get_reference_position(get_coords_as_pr(chrom, start, end, strand))
		regions_pos_pr = self.ref.pr

		# get nearest
		nearest = query_pos_pr.k_nearest(regions_pos_pr, k=k).as_df()
		nearest_regions = self.reindex(nearest.iloc[:,-2])
		nearest_regions["distance"] = nearest.iloc[:,-1].values
		return nearest_regions.sort_values("distance")

	def get_coords(self, phen_id): 
		return self.loc[phen_id]


def get_coords_as_pr(chrom, start, end, strand=None): 
	if strand == ".": 
		strand = None
	if strand: 
		return pr.from_dict({"Chromosome": [chrom], "Start": [start], "End": [end], "Strand": [strand]})
	else: 
		return pr.from_dict({"Chromosome": [chrom], "Start": [start], "End": [end]})

def get_reference_position(pr_region): 
	"""Get reference position for a set of coordinates"""
	if pr_region.stranded: 
		return pr_region.five_end()
	else: 
		region_s = pr_region.as_df().iloc[0]
		mid = (region_s["Start"]+region_s["End"]) // 2
		return pr.from_dict({"Chromosome": [region_s["Chromosome"]], "Start": [mid], "End": [mid+1]})


def pr_to_df(pr_obj, pos_cols=None): 
	"""Converts pyranges object to pandas dataframe."""
	if pos_cols is None: pos_cols = ["chrom", "start", "end", "strand"]

	df = pr_obj.as_df().rename(columns={"Chromosome": pos_cols[0], "Start": pos_cols[1], "End": pos_cols[2], "Strand": pos_cols[3]})

	if "gene_id" in df.columns: 
		df = df.set_index("gene_id")
		return Regions(df)
	elif "peak_id" in df.columns: 
		df = df.set_index("peak_id")
		return Regions(df)
	else: 
		logger.write("Format of `PyRanges` object not recognized")


	


## single operations

def get_point(regions, mode="tss"):  
	
	regions = regions.copy()

	if isinstance(regions, pd.DataFrame): 
		if mode == "tss": 
			assert (regions.strand != ".").all()
			# For rows on minus strand, set start to be one before end
			regions.loc[regions["strand"] == "-", "start"] = regions.loc[regions["strand"] == "-", "end"] - 1
			regions["end"] = regions["start"] + 1
		elif mode == "tes": 
			assert (regions.strand != ".").all()
			# For rows on minus strand, set end to be one after start
			regions.loc[regions["strand"] == "-", "end"] = regions.loc[regions["strand"] == "-", "start"] + 1
			regions["start"] = regions["end"] - 1
		elif mode == "midpoint": 
			regions["start"] = regions[["start", "end"]].mean(axis=1)
		else: 
			logger.write("No valid mode specified.")
			return 
		regions[["start", "end"]] = regions[["start", "end"]].astype(int)
		return regions

	else: 
		if mode == "tss": 
			assert (regions.Strand != ".").all()
			assert "Strand" in regions.columns
			regions = regions.five_end()
		elif mode == "tes": 
			assert (regions.Strand != ".").all()
			regions = regions.three_end()
		elif mode == "midpoint": 
			regions = regions.assign("midpoint", lambda df: ((df.Start + df.End)/2).astype(int))
		else: 
			logger.write("No valid mode specified.")
			return 
		return regions












## arithmetic
def join_regions_by_window(pr1, pr2, window=0): 
	"""Get pairs of regions within a specified window."""
	def get_distance(df): 
		"""Use 5' end if strand is specified (RNA). Otherwise use midpoint (ATAC)."""
		if (df.Strand == "+").all():	pos1 = df.Start
		elif (df.Strand == "-").all():	pos1 = df.End
		else:							pos1 = (df.Start+df.End)/2

		if (df.Strand_b == "+").all():	 pos2 = df.Start_b
		elif (df.Strand_b == "-").all(): pos2 = df.End_b
		else:					 		 pos2 = (df.Start_b+df.End_b)/2

		return (pos2-pos1).astype(int)

	out = pr1.join(pr2, slack=window, strandedness=False, how=None)
	out = out.assign("distance", get_distance)
	return out


def phenotype_snp_distance(pairs, phenotype_pos_df, variant_pos_df): 
	"""Get strand-specific distance between phenotype and variant position."""
	phenotype_pos = pairs["phenotype_id"].map(phenotype_pos_df["start"])
	variant_pos = pairs["variant_id"].map(variant_pos_df["pos"])

	distance = variant_pos - phenotype_pos

	if "strand" in phenotype_pos_df.columns: 
		distance.loc[pairs["phenotype_id"].map(phenotype_pos_df["strand"]) == "-"] *= -1

	return distance





