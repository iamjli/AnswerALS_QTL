#!/usr/bin/env python3

import numpy as np
import pandas as pd

from pathlib import Path

from src import logger



#----------------------------------------------------------------------------------------------------#
# Preprocess counts data
#----------------------------------------------------------------------------------------------------#

def write_phenotype_files(counts_data_obj, output_dir, prefix): 
	"""
	Output phenotypes in tensorqtl format. Produces 2 files: 
	 - tensorqtl phenotype input
	 - phenotype file that can be read by PEER
	"""
	regions, gtex = pd.DataFrame(counts_data_obj.regions).copy(), counts_data_obj.gtex.copy()
	regions = regions.reindex(gtex.index)

	# phenotype file for tensorqtl
	tensorqtl_phenotype_path = Path(output_dir) / f"{prefix}_gtex.bed.gz"
	if tensorqtl_phenotype_path.is_file(): 
		logger.write("Existing phenotype file found. Skipping...")
	else: 
		logger.write("Writing phenotype file for tensorqtl")
		tensorqtl_phenotype_df = get_counts_in_tensorqtl_fomat(regions, gtex)
		tensorqtl_phenotype_df.to_csv(tensorqtl_phenotype_path, sep="\t", index=False, compression="gzip")

	# PEER does not accept strand column, so write the same file without it
	PEER_phenotype_path = Path(output_dir) / f"{prefix}_gtex.for_PEER.bed.gz"
	if tensorqtl_phenotype_path.is_file(): 
		logger.write("Existing PEER file found. Skipping...")
	else: 
		logger.write("Writing phenotype file for PEER")
	PEER_phenotype_df = tensorqtl_phenotype_df.drop(columns=["strand"], errors="ignore")
	PEER_phenotype_df.to_csv(PEER_phenotype_path, sep="\t", index=False, compression="gzip")


def _get_counts_in_tensorqtl_fomat(phenotype_pos, phenotype_df): 
	"""
	Format dataframes
	"""
	assert phenotype_pos.index.equals(phenotype_df.index)

	# tensorqtl requires essentially a bed file prepended to the counts data
	tensorqtl_df = pd.concat([phenotype_pos, phenotype_df], axis=1)	
	tensorqtl_df["gene_id"] = phenotype_pos.index

	# reorder columns
	if "strand" in tensorqtl_df.columns: 
		metadata_cols = ["chrom", "start", "end", "gene_id", "strand"]
	else: 
		metadata_cols = ["chrom", "start", "end", "gene_id"]
	cols = metadata_cols + tensorqtl_df.columns.drop(metadata_cols).tolist()
	tensorqtl_df = tensorqtl_df[cols]

	# strip `chr`
	tensorqtl_df["chrom"] = tensorqtl_df["chrom"].str[3:]
	tensorqtl_df = tensorqtl_df.rename(columns={"chrom": "#chr"})

	return tensorqtl_df


def PEER_cmd(PEER_exec_path, phenotype_file, covariates_file, prefix, num_peer, output_dir):
	"""
	Command to execute PEER covariate correction. Be sure to use r-4.0.3
	""" 
	cmd = (
		"time Rscript {PEER_exec_path} {phenotype_file} {prefix} {num_peer} -o --covariates {covariates_file}"
	)
	return cmd
