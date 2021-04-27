#!/usr/bin/env python3

import numpy as np
import pandas as pd
from pathlib import Path

from src import base_dir, logger
from src import aals, hg38

__all__ = ["processed"]

#----------------------------------------------------------------------------------------------------#
# Generate processed data 
#----------------------------------------------------------------------------------------------------#
_paths = {
	"peak_TF_clusters": base_dir / "data/processed/peak_TF_clusters.txt"
}

def generate_cluster_membership(processed_path): 

	from src.analysis.load_data import FeatureData
	from src import data, hg38

	from src.analysis.annotate import get_regions_overlap_with_tfs
	from kmodes.kmodes import KModes

	atac = FeatureData(feature="atac", counts_prefix="atac_pysam_200bp", qtl_dir="2_caQTL_within_10kb_atac_pysam_200bp_210422")
	peak_tf_annos = get_regions_overlap_with_tfs(atac.regions, data.encode_tfs)

	km = KModes(n_clusters=4, init='Huang', n_init=12, verbose=1, n_jobs=12)
	tf_clusters_4 = pd.Series(km.fit_predict(peak_tf_annos), index=peak_tf_annos.index, name="TF_clust_4")

	km = KModes(n_clusters=10, init='Huang', n_init=24, verbose=1, n_jobs=24)
	tf_clusters_10 = pd.Series(km.fit_predict(peak_tf_annos), index=peak_tf_annos.index, name="TF_clust_10")

	cluster_membership = pd.concat([tf_clusters_4, tf_clusters_10], axis=1)

	cluster_membership.to_csv(processed_path, sep="\t", index=True, header=True)


def process_data(paths=_paths): 

	if paths["peak_TF_clusters"].is_file(): 
		logger.write("Existing peak TF cluster file exists. Skipping.")
	else: 
		generate_cluster_membership(paths["peak_TF_clusters"])



#----------------------------------------------------------------------------------------------------#
# Access processed data 
#----------------------------------------------------------------------------------------------------#
class ProcessedData: 
	"""Data associated with our samples such as sample paths and metadata."""

	def __init__(self, paths): 

		self.paths = paths

		self._peak_TF_clusters = None

	@property
	def peak_TF_clusters(self):
		if self._peak_TF_clusters is None: 
			self._peak_TF_clusters = _load_peak_TF_clusters(self.paths["peak_TF_clusters"])
		return self._peak_TF_clusters
	


#----------------------------------------------------------------------------------------------------#
# Load data 
#----------------------------------------------------------------------------------------------------#
def _load_peak_TF_clusters(path):
	"""Loads harmonized metadata."""
	return pd.read_csv(path, sep="\t", index_col=0)

#----------------------------------------------------------------------------------------------------#
processed = ProcessedData(_paths)