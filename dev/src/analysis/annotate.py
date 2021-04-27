#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pyranges as pr

from pathlib import Path

from src import base_dir
from src import data



_anno_dir = base_dir / "data/peak_annos"

anno_paths = {
	"genomic_annos": _anno_dir / "genomic_annos.txt", 
	"chip_annos": _anno_dir / "chip_annos.txt", 
	"omic_counts_annos": _anno_dir / "omic_counts_annos.txt", 
	"corr_qtl_annos": _anno_dir / "corr_qtl_annos.txt", 
}

#----------------------------------------------------------------------------------------------------#
# Genomic annotations
#----------------------------------------------------------------------------------------------------#
def load_genomic_annos(path): 
	pass


#----------------------------------------------------------------------------------------------------#
def get_regions_annotations(regions, annos_pr, return_counts=False): 
	"""Returns collapsed feature annotations."""
	annotated_peaks = _join_regions_with_annos(regions, annos_pr)

	annotation_counts = annotated_peaks.groupby(["peak_id", "Feature"]).size().unstack(-1).reindex(regions.index).fillna(0).astype(int)
	annotation_counts["intergenic"] = (annotation_counts.sum(axis=1) == 0).astype(int)
	
	if return_counts: return annotation_counts

	anno_matrix = annotation_counts > 0
	anno_matrix["tss_only"] = anno_matrix["tss"] & ~anno_matrix["tes"]
	anno_matrix["tes_only"] = anno_matrix["tes"] & ~anno_matrix["tss"]
	anno_matrix["tss_tes"]  = anno_matrix["tss"] &  anno_matrix["tes"]

	unique_annos = pd.concat([
	    pd.Series("tss_only",   anno_matrix.index[anno_matrix["tss_only"]]),
	    pd.Series("tes_only",   anno_matrix.index[anno_matrix["tes_only"]]),
	    pd.Series("tss_tes",    anno_matrix.index[anno_matrix["tss_tes"]]), 
	    pd.Series("exon",       anno_matrix.index[anno_matrix["exon"]       & ~anno_matrix[["tss", "tes"]].any(axis=1)]), 
	    pd.Series("intron",     anno_matrix.index[anno_matrix["intron"]     & ~anno_matrix[["tss", "tes", "exon"]].any(axis=1)]), 
	    pd.Series("intergenic", anno_matrix.index[anno_matrix["intergenic"] & ~anno_matrix[["tss", "tes", "exon", "intron"]].any(axis=1)])
	]).rename("unique_annos")

	return anno_matrix, unique_annos


def _join_regions_with_annos(regions, annos_pr): 
	if annos_pr is None:
		annos_pr = data.gencode

	annos_merged_pr = annos_pr.merge(by=["Feature", "gene_id"], strand=False) # Merges overlapping intervals by feature_col and unstrand
	joined = regions.pr.join(annos_merged_pr, suffix="_gencode", strandedness=None, report_overlap=True).as_df()

	# remove ENSG version number
	joined["gene_id"] = joined["gene_id"].str.split(".").str[0]

	return joined



#----------------------------------------------------------------------------------------------------#
# ChIP annotations
#----------------------------------------------------------------------------------------------------#

def get_regions_overlap_with_tfs(regions, encode_pr): 
	"""Returns boolean matrix indicating whether regions overlaps with an ENCODE TF."""
	joined = regions.pr.join(encode_pr).as_df()
	is_overlap = joined.groupby(["peak_id", "Name"]).size().unstack(-1, fill_value=0) > 0
	is_overlap = is_overlap.reindex(regions.index, fill_value=False)
	return is_overlap


#----------------------------------------------------------------------------------------------------#
# Omic track features
#----------------------------------------------------------------------------------------------------#




#----------------------------------------------------------------------------------------------------#
# Omic correlation features
#----------------------------------------------------------------------------------------------------#

