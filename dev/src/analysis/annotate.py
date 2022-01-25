#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr

from src import base_dir
from src.load import data, hg38


_anno_paths = {
	"genomic_annos": base_dir / "data/peak_annos/genomic_annos.txt", 
	"unique_genomic_annos": base_dir / "data/peak_annos/unique_genomic_annos.txt", 
	"chip_annos": base_dir / "data/peak_annos/chip_annos.txt", 
	"omic_annos": base_dir / "data/peak_annos/omic_annos.txt", 
	"qtl_annos": base_dir / "data/peak_annos/qtl_annos.txt", 
	"tf_cluster_annos": base_dir / "data/peak_annos/peak_TF_clusters.txt",
}

class PeakAnnos: 

	def load_genomic_annos(): 
		g_ann = pd.read_csv(_anno_paths["genomic_annos"], sep="\t", index_col=0)
		g_ann["unique_annos"] = pd.Categorical(g_ann["unique_annos"], categories=["intergenic", "intron", "exon", "tes_only", "tss_tes", "tss_only"], ordered=True)
		return g_ann

	def load_unique_genomic_annos(): 
		g_ann_unique = pd.read_csv(_anno_paths["unique_genomic_annos"], sep="\t", index_col=0)
		return g_ann_unique.iloc[:, 0].astype("category")

	def load_chip_annos(): 
		return pd.read_csv(_anno_paths["chip_annos"], sep="\t", index_col=0)

	def load_omic_annos(): 
		return pd.read_csv(_anno_paths["omic_annos"], sep="\t", index_col=0)

	def load_tf_clusters_annos(): 
		c_ann_clusters = pd.read_csv(_anno_paths["tf_cluster_annos"], sep="\t", index_col=0)
		for col in c_ann_clusters.columns: 
			c_ann_clusters[col] = c_ann_clusters[col].astype("category")
		return c_ann_clusters

	def load_qtl_annos(): 
		return pd.read_csv(_anno_paths["qtl_annos"], sep="\t", index_col=0)



#----------------------------------------------------------------------------------------------------#
# Genomic annotations
#----------------------------------------------------------------------------------------------------#
def load_genomic_annos(path): 
	pass


#----------------------------------------------------------------------------------------------------#
def annotate_peaks(regions): 
	_tss_peaks = regions.pr.overlap(hg38.gencode_gtf.features.tss().slack(1000)).peak_id.unique()
	_tes_peaks = regions.pr.overlap(hg38.gencode_gtf.features.tes().slack(1000)).peak_id.unique()
	_tss_and_tes_peaks = set(_tss_peaks) & set(_tes_peaks)
	_splice_peaks = regions.pr.overlap(hg38.splice_junctions).peak_id.unique()
	_exon_peaks = regions.pr.overlap(hg38.exons).peak_id.unique()
	_intron_peaks = regions.pr.overlap(hg38.introns).peak_id.unique()

	peak_annos = pd.concat([
		pd.Series(data="tss_tes", index=_tss_and_tes_peaks), 
		pd.Series(data="tss", index=_tss_peaks), 
		pd.Series(data="tes", index=_tes_peaks), 
		pd.Series(data="splice", index=_splice_peaks), 
		pd.Series(data="exon", index=_exon_peaks), 
		pd.Series(data="intron", index=_intron_peaks), 
	])
	peak_annos = peak_annos[~peak_annos.index.duplicated()].reindex(regions.index, fill_value="intergenic")
	return peak_annos


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
# TF clusters
#----------------------------------------------------------------------------------------------------#
def generate_cluster_membership(processed_path): 

	from src.analysis.load_data import FeatureData
	from src import data, hg38

	from src.analysis.annotate import get_regions_overlap_with_tfs
	from kmodes.kmodes import KModes

	atac = FeatureData(feature="atac", counts_prefix="atac_pysam_200bp", qtl_dir="2_caQTL_within_10kb_atac_pysam_200bp_210422")
	peak_tf_annos = get_regions_overlap_with_tfs(atac.regions, data.encode_tfs)

	km = KModes(n_clusters=4, init='Huang', n_init=12, verbose=1, n_jobs=12)
	tf_clusters_4 = pd.Series(km.fit_predict(peak_tf_annos), index=peak_tf_annos.index, name="TF_clust_4")

	km = KModes(n_clusters=6, init='Huang', n_init=12, verbose=1, n_jobs=12)
	tf_clusters_4 = pd.Series(km.fit_predict(peak_tf_annos), index=peak_tf_annos.index, name="TF_clust_6")

	km = KModes(n_clusters=8, init='Huang', n_init=12, verbose=1, n_jobs=12)
	tf_clusters_4 = pd.Series(km.fit_predict(peak_tf_annos), index=peak_tf_annos.index, name="TF_clust_8")

	km = KModes(n_clusters=10, init='Huang', n_init=24, verbose=1, n_jobs=24)
	tf_clusters_10 = pd.Series(km.fit_predict(peak_tf_annos), index=peak_tf_annos.index, name="TF_clust_10")

	cluster_membership = pd.concat([tf_clusters_4, tf_clusters_6, tf_clusters_8, tf_clusters_10], axis=1)

	cluster_membership.to_csv(processed_path, sep="\t", index=True, header=True)
#----------------------------------------------------------------------------------------------------#
# Omic track features
#----------------------------------------------------------------------------------------------------#




#----------------------------------------------------------------------------------------------------#
# Omic correlation features
#----------------------------------------------------------------------------------------------------#

