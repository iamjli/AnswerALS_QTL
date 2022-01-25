#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pysam

from pathos import multiprocessing

from src import logger
from src.load import aals, data
from src.query import vcf_utils

_gt_to_dosage_dict = {'0/0':0, '0/1':1, '1/1':2, './.':np.nan,
					 '0|0':0, '0|1':1, '1|0':1, '1|1':2, '.|.':np.nan}

class QueryVCF: 

	def __init__(self, vcf_path=None, validate=True): 

		self.vcf_path = vcf_path
		self.sample_names = vcf_utils._get_samples_names_from_file(self.vcf_path)

		if validate: 
			self._validate()

	def _validate(self): 
		assert self.sample_names.equals(aals.sample_names)

	@classmethod
	def load(cls): 
		vcf_obj = cls(aals.paths["vcf"])
		return vcf_obj

	def query_interval(self, chrom, start, end, return_metadata=False): 
		"""
		Returns datafame of genotypes in an interval. Optionally returns metadata.
		"""
		return vcf_utils._query_interval_from_file(self.vcf_path, chrom, start, end, self.sample_names, return_metadata)

	def query_regions(self, regions, return_metadata=False): 
		if len(regions) < 200: 
			return vcf_utils._query_regions_from_file(self.vcf_path, regions, self.sample_names, return_metadata)
		else: 
			logger.write("Parallelizing...")
			return vcf_utils._query_many_regions_from_file(self.vcf_path, regions, self.sample_names, return_metadata)

	def query_pos(self, chrom, pos, return_metadata=False): 
		"""Query at a specific position."""
		results = self.query_interval(chrom, pos, pos, return_metadata)
		if return_metadata: 
			return results[0].iloc[0], results[1].iloc[0]
		else: 
			return results.iloc[0]

	def query_rsid(self, rsid, return_metadata=False, chr_prefix=True): 
		"""Query rsID genotypes."""
		chrom, pos, *_ = data.rsid.loc[rsid]
		if chr_prefix: 
			chrom = chrom if chrom.startswith('chr') else 'chr' + str(chrom)
		else: 
			chrom = chrom[3:] if chrom.startswith('chr') else chrom

		if return_metadata: 
			genotypes_s, metadata_s = self.query_pos(chrom, pos, return_metadata=True)
			genotypes_s.name = rsid
			return genotypes_s, metadata_s
		else: 
			genotypes_s = self.query_pos(chrom, pos, return_metadata=False)
			genotypes_s.name = rsid
			return genotypes_s

	def query_rsids(self, rsids, return_metadata=False, fast=False, chr_prefix=True): 

		rsid_regions = data.rsid.loc[rsids]
		rsid_regions["start"] = rsid_regions["pos"].astype(int)
		rsid_regions["end"] = rsid_regions["start"] + 1
		
		if chr_prefix and not rsid_regions.iloc[0]['chrom'].startswith('chr'):
			rsid_regions['chrom'] = 'chr' + rsid_regions['chrom'].astype(str)
		elif not chr_prefix and rsid_regions.iloc[0]['chrom'].startswith('chr'): 
			rsid_regions['chrom'] = rsid_regions['chrom'].str[3:] 

		return self.query_regions(rsid_regions, return_metadata)

		# if fast: 



	# def test_maf_dist(self, rsid=None, gt=None): 

	# 	if rsid is not None: 
	# 		gt = self.get_genotype(rsid)

	# 	# dataframe: index=["ALS", "CTR"], columns=[0,1,2]
	# 	maf_counts = gt.groupby(self.metadata["condition"]).value_counts().unstack(-1)
	# 	# expected maf counts: taken by multiplying the propoprtion of ALS/CTR and total maf counts 
	# 	proportion = maf_counts.sum(axis=1) / maf_counts.sum(axis=1).sum()  # proportion of ALS vs. CTR
	# 	exp_counts = np.outer(proportion, maf_counts.sum(axis=0))

	# 	results = stats.chisquare(maf_counts, f_exp=exp_counts, axis=0)
	# 	return pd.Series(data=results.pvalue, index=maf_counts.columns)



# def open_vcf_file(vcf_file, threads=2): 
# 	return pysam.TabixFile(str(vcf_file), threads=threads, parser=pysam.asVCF())

# def _get_samples_names_from_file(file): 
# 	with open_vcf_file(file) as tbx: 
# 		sample_names = pd.Index(tbx.header[-1].split("\t")[9:])
# 	return sample_names

# def _query_interval_from_file(file, chrom, start, end, sample_names, return_metadata=False): 
# 	"""
# 	Returns datafame of genotypes in an interval. Optionally returns metadata.
# 	"""
# 	metadata, genotypes = [], []
# 	with open_vcf_file(file) as tbx:
# 		for row in tbx.fetch(chrom, start-1, end): 
# 			metadata.append([row.id, row.contig, row.pos+1, row.ref, row.alt])
# 			genotypes.append([_gt_to_dosage_dict[gt] for gt in row[:len(sample_names)]])

# 	metadata_df = pd.DataFrame(metadata, columns=["variant_id", "chrom", "pos", "ref", "alt"]).set_index("variant_id")
# 	genotypes_df = pd.DataFrame(genotypes, index=metadata_df.index, columns=sample_names)
# 	assert metadata_df.index.is_unique

# 	if return_metadata: return genotypes_df, metadata_df
# 	return genotypes_df

# def _query_regions_from_file(file, regions, sample_names, return_metadata=False): 
# 	"""
# 	Returns datafame of genotypes in an interval. Optionally returns metadata.
# 	"""
# 	metadata, genotypes = [], []
# 	with open_vcf_file(file) as tbx:
# 		# for _,region in regions.iterrows(): 
# 		# 	chrom, start, end = region[["chrom", "start", "end"]]
# 		for region in regions: 
# 			chrom, start, end = region
# 			for row in tbx.fetch(chrom, start-1, end): 
# 				metadata.append([row.id, row.contig, row.pos+1, row.ref, row.alt])
# 				genotypes.append([_gt_to_dosage_dict[gt] for gt in row[:len(sample_names)]])

# 	metadata_df = pd.DataFrame(metadata, columns=["variant_id", "chrom", "pos", "ref", "alt"]).set_index("variant_id")
# 	genotypes_df = pd.DataFrame(genotypes, index=metadata_df.index, columns=sample_names)
# 	# assert metadata_df.index.is_unique

# 	if return_metadata: return genotypes_df, metadata_df
# 	return genotypes_df
