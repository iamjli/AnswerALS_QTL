#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pysam

from pathos import multiprocessing

from src import logger


_gt_to_dosage_dict = {'0/0':0, '0/1':1, '1/1':2, './.':np.nan,
					 '0|0':0, '0|1':1, '1|0':1, '1|1':2, '.|.':np.nan}



def open_vcf_file(vcf_file, threads=2): 
	return pysam.TabixFile(str(vcf_file), threads=threads, parser=pysam.asVCF())

def _get_samples_names_from_file(file): 
	with open_vcf_file(file) as tbx: 
		sample_names = pd.Index(tbx.header[-1].split("\t")[9:])
	return sample_names

def _query_interval_from_file(file, chrom, start, end, sample_names, return_metadata=False): 
	"""
	Returns datafame of genotypes in an interval. Optionally returns metadata.
	"""
	metadata, genotypes = [], []
	with open_vcf_file(file) as tbx:
		for row in tbx.fetch(chrom, start-1, end): 
			metadata.append([row.id, row.contig, row.pos+1, row.ref, row.alt])
			genotypes.append([_gt_to_dosage_dict[gt.split(':')[0]] for gt in row[:len(sample_names)]])

	metadata_df = pd.DataFrame(metadata, columns=["variant_id", "chrom", "pos", "ref", "alt"]).set_index("variant_id")
	genotypes_df = pd.DataFrame(genotypes, index=metadata_df.index, columns=sample_names)
	assert metadata_df.index.is_unique

	if return_metadata: return genotypes_df, metadata_df
	return genotypes_df

def _query_regions_from_file(file, regions, sample_names, return_metadata=False): 
	"""
	Returns datafame of genotypes in an interval. Optionally returns metadata.
	"""
	metadata, genotypes = [], []
	with open_vcf_file(file) as tbx:
		for _,region in regions.iterrows(): 
			chrom, start, end = region[["chrom", "start", "end"]]
			for row in tbx.fetch(chrom, start-1, end): 
				metadata.append([row.id, row.contig, row.pos+1, row.ref, row.alt])
				genotypes.append([_gt_to_dosage_dict[gt.split(':')[0]] if gt.split(':')[0] in _gt_to_dosage_dict else 2-gt.split(':')[0].count('0') for gt in row[:len(sample_names) ]])

	metadata_df = pd.DataFrame(metadata, columns=["variant_id", "chrom", "pos", "ref", "alt"]).set_index("variant_id")
	genotypes_df = pd.DataFrame(genotypes, index=metadata_df.index, columns=sample_names)
	# assert metadata_df.index.is_unique

	if return_metadata: return genotypes_df, metadata_df
	return genotypes_df

def _query_many_regions_from_file(file, regions, sample_names, return_metadata=False): 

	split_regions = np.array_split(regions.copy(), 16)
	args = ((str(file), _regions, sample_names, return_metadata) for _regions in split_regions)
	with multiprocessing.ProcessingPool(16) as pool: 
		results = pool.map(lambda args: _query_regions_from_file(*args), args)
	if return_metadata: 
		snp_results, metadata_results = zip(*results)
		snp_results = pd.concat(snp_results)
		snp_results = snp_results[~snp_results.index.duplicated(keep=False)]
		metadata_results = pd.concat(metadata_results)
		metadata_results = metadata_results[~metadata_results.index.duplicated(keep=False)]
		assert (metadata_results.index.is_unique) & (snp_results.index.is_unique)
		return snp_results, metadata_results
	else: 
		snp_results = pd.concat(results)
		snp_results = snp_results[~snp_results.index.duplicated(keep=False)]
		assert snp_results.index.is_unique
		return snp_results


def _query_regions_from_project_mine_file(file, regions, sample_names): 
	"""
	Returns datafame of genotypes in an interval. Optionally returns metadata.
	"""
	metadata, genotypes = [], []
	with open_vcf_file(file) as tbx:
		for _,region in regions.iterrows(): 
			chrom, start, end = region[["chrom", "start", "end"]]
			for row in tbx.fetch(chrom, start-1, end): 
				gt = row[0].split(":")[0]
				gt = _gt_to_dosage_dict[gt] if gt in _gt_to_dosage_dict else np.nan
				metadata.append([row.id, row.contig, row.pos+1, row.ref, row.alt, gt])

	metadata_df = pd.DataFrame(metadata, columns=["variant_id", "chrom", "pos", "ref", "alt", "gt"])
	return metadata_df