#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pyranges as pr




def df_to_pr(df, pos_cols=None, include_index=True):
	"""Converts pandas dataframe to pyranges object."""
	df = df.copy()

	if df.index.name == "variant_id":  # if this is a variant dataframe
		assert "end" not in df.columns
		df.rename(columns={"chrom": "Chromosome", "start": "Start", "pos": "Start"}, inplace=True)
		df = df.assign(End=(df.Start + 1).astype(int)).reset_index()
		df = df[ ["Chromosome", "Start", "End"] + df.columns.drop(["Chromosome", "Start", "End"]).tolist() ]
		return pr.PyRanges(df)

	else: 
		if pos_cols is None: pos_cols = ["chrom", "start", "end", "strand"]

		if include_index: 
			df.reset_index(inplace=True)
			df = df[ pos_cols + df.columns.drop(pos_cols).tolist() ]  # reorder columns

		df.rename(columns={pos_cols[0]: "Chromosome", pos_cols[1]: "Start", pos_cols[2]: "End", pos_cols[3]: "Strand"}, inplace=True)
		return pr.PyRanges(df)


def pr_to_df(pr_obj, pos_cols=None, index=None): 
	"""Converts pyranges object to pandas dataframe."""
	if pos_cols is None: pos_cols = ["chrom", "start", "end", "strand"]

	df = pr_obj.as_df().rename(columns={"Chromosome": pos_cols[0], "Start": pos_cols[1], "End": pos_cols[2], "Strand": pos_cols[3]})
	if index is not None: 
		df.set_index(index, inplace=True)
		df = df[ pos_cols + df.columns.drop(pos_cols).tolist() ]  # reorder columns

	return df



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
			print("No valid mode specified.")
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
			print("No valid mode specified.")
			return 
		return regions








	# is_df = isinstance(regions, pd.DataFrame)

	# if is_df: 
	# 	index_name = regions.index.name
	# 	regions = df_to_pr(regions)

	# if mode == "tss": 
	# 	assert "Strand" in regions.columns
	# 	regions = regions.five_end()
	# elif mode == "tes": 
	# 	assert "Strand" in regions.columns
	# 	regions = regions.three_end()
	# else: 
	# 	regions = regions.assign("midpoint", lambda df: ((df.Start + df.End)/2).astype(int))

	# if is_df: 
	# 	return pr_to_df(regions, index=index_name)
	# else: 
	# 	return regions

# 	if is_pyranges: return bed_to_pr(regions_df)
# 	else: return regions_df

# def get_midpoint(regions): 

# 	is_pyranges = isinstance(regions, pr.PyRanges)

# 	if is_pyranges: 
# 		regions_df = df_to_pr(regions)
# 	else: 
# 		regions_df = regions.copy()

# 	regions_df.loc["midpoint"] = regions_df.pos_df[["start", "end"]].mean(axis=1).astype(int)

# 	if is_pyranges: return bed_to_pr(regions_df)
# 	else: return regions_df

# def get_



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





