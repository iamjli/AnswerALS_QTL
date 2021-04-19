#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pyranges as pr

from src import data


#----------------------------------------------------------------------------------------------------#
# Annotate dataframes
#----------------------------------------------------------------------------------------------------#
@pd.api.extensions.register_dataframe_accessor("anno")
class Annotator:

	def __init__(self, df): 

		self._df = df

	def gene_symbol(self, col="ensg_id"): 
		pass











#--------------------------------------------------------------------------------------------------#
def annotate_peaks(regions_pr, annos_pr=None): 
	if annos_pr is None:
		annos_pr = data.gencode

	annos_merged_pr = annos_pr.merge(by=["Feature", "gene_id"], strand=False) # Merges overlapping intervals by feature_col and unstrand
	joined = regions_pr.join(annos_merged_pr, strandedness=None, report_overlap=True).as_df()

	# remove ENSG version number
	joined["gene_id"] = joined["gene_id"].str.split(".").str[0]

	return joined

def get_peak_annotations(regions, mode=None, annos_pr=None): 
	"""Returns collapsed feature annotations."""
	if annos_pr is None:
		annos_pr = data.gencode

	annotated_peaks = annotate_peaks(regions.bed.as_pr(), annos_pr)

	annotation_counts = annotated_peaks.groupby(["peak_id", "Feature"]).size().unstack(-1).reindex(regions.index).fillna(0).astype(int)
	annotation_counts["intergenic"] = (annotation_counts.sum(axis=1) == 0).astype(int)
	
	if mode == "counts_matrix": 
		return annotation_counts
	elif mode == "bool_matrix": 
		return annotation_counts > 0
	else mode == "unique":
		peak_features = annotation_counts > 0
		unique_features = pd.concat([
		    pd.Series("tss_only",   peak_features.index[peak_features["tss"]        & ~peak_features["tes"]]), 
		    pd.Series("tes_only",   peak_features.index[peak_features["tes"]        & ~peak_features["tss"]]), 
		    pd.Series("tss_tes",    peak_features.index[peak_features["tss"]        & peak_features["tes"]]), 
		    pd.Series("exon",       peak_features.index[peak_features["exon"]       & ~peak_features[["tss", "tes"]].any(axis=1)]), 
		    pd.Series("intron",     peak_features.index[peak_features["intron"]     & ~peak_features[["tss", "tes", "exon"]].any(axis=1)]), 
		    pd.Series("intergenic", peak_features.index[peak_features["intergenic"] & ~peak_features[["tss", "tes", "exon", "intron"]].any(axis=1)])
		])
		return unique_features.reindex(regions.index)

