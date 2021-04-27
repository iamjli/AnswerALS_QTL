#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pysam

from src import data

gt_to_dosage_dict = {'0/0':0, '0/1':1, '1/1':2, './.':np.nan,
					 '0|0':0, '0|1':1, '1|0':1, '1|1':2, '.|.':np.nan}


class QueryVCF: 

	def __init__(self, vcf_path=None): 

		self.vcf_path = vcf_path
		# self.tbx = open_vcf_file(self.vcf_path)

		self.sample_names = _get_samples_names_from_file(self.vcf_path)
		self.n_samples = len(self.sample_names)

	@classmethod
	def load(cls): 
		from src import aals
		vcf_obj = cls(aals.paths["vcf"])
		assert vcf_obj.sample_names.equals(aals.sample_names)
		return vcf_obj

	def query_interval(self, chrom, start, end, return_metadata=False): 

		metadata, genotypes = [], []
		with open_vcf_file(self.vcf_path) as tbx:
			for row in tbx.fetch(chrom, start-1, end): 
				metadata.append([row.id, row.contig, row.pos+1, row.ref, row.alt])
				genotypes.append([gt_to_dosage_dict[gt] for gt in row[:self.n_samples]])

		metadata_df = pd.DataFrame(metadata, columns=["variant_id", "chrom", "pos", "ref", "alt"]).set_index("variant_id")
		genotypes_df = pd.DataFrame(genotypes, index=metadata_df.index, columns=self.sample_names)
		assert metadata_df.index.is_unique

		if return_metadata: return genotypes_df, metadata_df
		return genotypes_df

	def query_pos(self, chrom, pos, **kwargs): 
		return self.query_interval(chrom, pos, pos, **kwargs)

	def query_rsid(self, rsid, return_metadata=False): 
		chrom, pos, *_ = data.rsid.loc[rsid]
		if return_metadata: 
			genotypes_df, metadata_df = self.query_pos(chrom, pos, return_metadata=True)
			return genotypes_df.loc[rsid], metadata_df.loc[rsid]
		else: 
			genotypes_df = self.query_pos(chrom, pos, return_metadata=False)
			return genotypes_df.loc[rsid]

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



def open_vcf_file(vcf_file, threads=2): 
	return pysam.TabixFile(str(vcf_file), threads=threads, parser=pysam.asVCF())

def _get_samples_names_from_tbx(tbx): 
	return pd.Index(tbx.header[-1].split("\t")[9:])

def _get_samples_names_from_file(file): 
	with open_vcf_file(file) as tbx: 
		sample_names = pd.Index(tbx.header[-1].split("\t")[9:])
	return sample_names