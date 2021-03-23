#!/usr/bin/env python3
import re

import pandas as pd
import numpy as np
import pyranges as pr
import pysam

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from . import BASE_DIR, CHROMS, logger

# from .repository import _load_metadata, _load_bam_paths


gt_to_dosage_dict = {'0/0':0, '0/1':1, '1/1':2, './.':np.nan,
					 '0|0':0, '0|1':1, '1|0':1, '1|1':2, '.|.':np.nan}

class Genomic: 

	def __init__(self, vcf_path=None, metadata=None, bam_paths=None, rsid=None, R=None, initialize=False): 

		self.R = R

		if self.R is not None: 
			self.vcf_path  = self.R.vcf_path
			self.metadata  = self.R.metadata
			self.bam_paths = self.R.bam_paths
			# self.rsid      = self.R.rsid
		else: 
			from .repository import _load_metadata, _load_bam_paths, _load_rsid
			self.metadata  = _load_metadata()
			self.bam_paths = _load_bam_paths()
			# self.rsid      = _load_rsid()

		self.tbx = pysam.TabixFile(str(self.vcf_path))
		# self._chr_prefix = "chr" in self.tbx.references[0]

		# Ensure we're working with the same samples
		self.sample_names = pd.Index(self.tbx.header[-1].split("\t")[9:])
		assert self.sample_names.isin(self.metadata.index).all()
		assert self.sample_names.isin(self.bam_paths.index).all()
		self.metadata = self.metadata.reindex(self.sample_names)
		self.bam_paths = self.bam_paths.reindex(self.sample_names)

		self._rna_bams_pysam = None
		self._atac_bams_pysam = None

		self._rna_library_sizes = None
		self._atac_library_sizes = None

		if initialize: 
			self.rna_bams_pysam
			self.atac_bams_pysam


	## RNA and ATAC

	@property
	def rna_bams_pysam(self): 
		if self._rna_bams_pysam is None: 
			_rna_bams_pysam = {}
			for i,(guid,path) in enumerate(self.bam_paths["atac_bam"].items()): 
				logger.update("Loading {} of {} RNA bams as pysam objects...".format(i, len(self.bam_paths)))
				# _rna_bams_pysam[guid] = pysam.AlignmentFile(path, "rb")
			logger.flush()
			self._rna_bams_pysam = _rna_bams_pysam
			return self._rna_bams_pysam
		else: 
			return self._rna_bams_pysam

	@property
	def atac_bams_pysam(self): 
		if self._atac_bams_pysam is None: 
			_atac_bams_pysam = {}
			for i,(guid,path) in enumerate(self.bam_paths["atac_bam"].items()): 
				logger.update("Loading {} of {} RNA bams as pysam objects...".format(i, len(self.bam_paths)))
				# _atac_bams_pysam[guid] = pysam.AlignmentFile(path, "rb")
			logger.flush()
			self._atac_bams_pysam = _atac_bams_pysam
			return self._atac_bams_pysam
		else: 
			return self._atac_bams_pysam

	@property
	def rna_library_sizes(self):
		if self._rna_library_sizes is None: 
			sizes = {}
			for i,(guid,path) in enumerate(self.bam_paths["rna_bam"].iteritems()): 
				logger.update("Scanning {} of {} RNA files for library_sizes...".format(i, len(self.bam_paths)))
				results = pysam.idxstats(path)
				sizes[guid] = np.array([row.split("\t")[2:] for row in results.rstrip("\n").split("\n")], dtype=int).sum()
			logger.flush()
			self._rna_library_sizes = pd.Series(sizes)
		return self._rna_library_sizes

	@property
	def atac_library_sizes(self):
		if self._atac_library_sizes is None: 
			sizes = {}
			for i,(guid,path) in enumerate(self.bam_paths["atac_bam"].iteritems()): 
				logger.update("Scanning {} of {} ATAC files for library_sizes...".format(i, len(self.bam_paths)))
				results = pysam.idxstats(path)
				sizes[guid] = np.array([row.split("\t")[2:] for row in results.rstrip("\n").split("\n")], dtype=int).sum()
			logger.flush()
			self._atac_library_sizes = pd.Series(sizes)
		return self._atac_library_sizes
	

	def _match_chrom_prefixes(self, bam_pysam, chrom): 
		# Check format of chromosomes (i.e. if `chr1` or just `1`)
		if "chr" in bam_pysam.references[0]: 
			if "chr" not in chrom: 
				chrom = "chr" + str(chrom)
		else: 
			if "chr" in chrom: 
				chrom = chrom[3:]
		return chrom

	def _get_coverage(self, bam_pysam, chrom, start, end, **kwargs): 
		chrom = self._match_chrom_prefixes(bam_pysam, chrom)

		try:
			positions, pileup = zip(*[(result.pos, result.n) for result in bam_pysam.pileup(chrom, start, end)])
			coverage = pd.Series(data=pileup, index=positions)
			return coverage.reindex(np.arange(start, end), fill_value=0)
		except ValueError:
			return pd.Series(data=0, index=np.arange(start, end))

	def _get_coverages(self, bam_pysam_dic, chrom, start, end, **kwargs): 
		coverages = {}
		for i,(guid,bam) in enumerate(bam_pysam_dic.items()): 
			logger.update("Loading coverage for {} of {} bams...".format(i, len(bam_pysam_dic)))
			coverages[guid] = self._get_coverage(bam, chrom, start, end)
		logger.flush()
		results = pd.DataFrame.from_dict(coverages)
		return Coverage(results)

	def get_rna_coverage(self, chrom=None, start=None, end=None, region=None, w=0, norm_factors=1, **kwargs): 
		"""Returns dataframe of read pileup with rows indexed by position and columns indexed by GUID."""
		if region is not None: 
			chrom, start, end = re.search(r'(.*):(\d*)-(\d*)', region).groups()
			start, end = int(start), int(end)
		start, end = start-w, end+w
		results = pd.DataFrame.from_dict({ guid: self._get_coverage(bam, chrom, start, end) for guid,bam in self.rna_bams_pysam.items() })
		results = results.div(norm_factors, axis=1)
		return Coverage(results)

	def get_atac_coverage(self, chrom=None, start=None, end=None, region=None, w=0, norm_factors=1, **kwargs): 
		"""Returns dataframe of read pileup with rows indexed by position and columns indexed by GUID."""
		if region is not None: 
			chrom, start, end = re.search(r'(.*):(\d*)-(\d*)', region).groups()
			start, end = int(start), int(end)
		start, end = start-w, end+w

		coverages = self._get_coverages(self.atac_bams_pysam, chrom, start, end)
		coverages = coverages.div(norm_factors, axis=1)

		return Coverage(coverages)

	def plot_coverage_around_snp(self, omic, variant_id, w=2500, norm_factors=1, ax=None): 

		chrom,pos,_,_ = self.rsid.loc[variant_id]
		start,end = pos-w, pos+w

		if omic == "rna": 
			coverage = self.get_rna_coverage(chrom, start, end)
		else: 
			coverage = self.get_atac_coverage(chrom, start, end)

		coverage.index = (coverage.index - (start+end)/2).astype(int)
		coverage = coverage.div(norm_factors, axis=1)

		genotypes = self.get_genotype(variant_id)

		sns.lineplot(data=coverage.groupby(genotypes, axis=1).mean(), ax=ax)
		ax.axvline(pos - (start+end)/2)

	def plot_coverage_around_snp_by_condition(self, omic, variant_id, w=2500, norm_factors=1, ax=None): 

		chrom,pos,_,_ = self.rsid.loc[variant_id]
		start,end = pos-w, pos+w

		if omic == "rna": 
			coverage = self.get_rna_coverage(chrom, start, end)
		else: 
			coverage = self.get_atac_coverage(chrom, start, end)

		coverage.index = (coverage.index - (start+end)/2).astype(int)
		coverage = coverage.div(norm_factors, axis=1)

		genotypes = self.get_genotype(variant_id)

		sns.lineplot(data=coverage.loc[:,self.metadata["condition"] == "ALS"].groupby(genotypes, axis=1).mean(), ax=ax[0])
		ax[0].axvline(pos - (start+end)/2)

		sns.lineplot(data=coverage.loc[:,self.metadata["condition"] == "CTR"].groupby(genotypes, axis=1).mean(), ax=ax[1])
		ax[1].axvline(pos - (start+end)/2)


	## VCF

	def get_snp_metadata(self, chrom, start, end, **kwargs): 

		snps = [ row.split("\t")[:5] for row in self.tbx.fetch(chrom, start-1, end) ]
		snps = pd.DataFrame(snps, columns=["chrom", "start", "variant_id", "ref", "alt"])
		snps = snps[snps["variant_id"] != "."].set_index("variant_id")
		return snps

	def get_genotype(self, rsid, as_series=True): 
		
		gt = self.get_genotypes([rsid], as_df=False)[0]
		if as_series: 
			return pd.Series(data=gt, index=self.sample_names, name=rsid)
		else: 
			return gt

	def get_genotypes(self, rsid, as_df=True):
		"""Gets genotypes by rsid."""
		if as_df: 
			snps_012 = { rsid:self.get_genotypes_at_pos(row["chrom"], row["pos"], as_series=False) for rsid,row in self.rsid.loc[rsid].iterrows() }
			return pd.DataFrame.from_dict(snps_012, orient="index", columns=self.sample_names)
		else: 
			return np.array([ self.get_genotypes_at_pos(row["chrom"], row["pos"], as_series=False) for _,row in self.rsid.loc[rsid].iterrows() ])

	def get_genotypes_at_pos(self, chrom, pos, as_series=True, **kwargs): 

		result = next(self.tbx.fetch(chrom, pos-1, pos))
		values = result.split("\t")
		snp_id = values[2]
		snp_012 = np.array([gt_to_dosage_dict[gt] for gt in values[9:]])
		# snp_012 = [np.array(gt.replace(".", "0").split("/"), dtype=int).sum() for gt in values[9:]]

		if as_series: 
			return pd.Series(data=snp_012, index=self.sample_names, name=snp_id)
		else: 
			return snp_012

	def get_genotypes_in_region(self, chrom, start, end, w=0, as_df=True, **kwargs): 

		start = start - w
		end = end + w

		snps_012 = {} 
		for result in self.tbx.fetch(chrom, start-1, end): 
			values = result.split("\t")
			snp_id = values[2]
			if snp_id != ".": 
				snps_012[snp_id] = [np.array(gt.replace(".", "0").split("/"), dtype=int).sum() for gt in values[9:]]

		if as_df: 
			return pd.DataFrame.from_dict(snps_012, orient="index", columns=self.sample_names)
		else: 
			return snps_012

	def test_maf_dist(self, rsid=None, gt=None): 

		if rsid is not None: 
			gt = self.get_genotype(rsid)

		# dataframe: index=["ALS", "CTR"], columns=[0,1,2]
		maf_counts = gt.groupby(self.metadata["condition"]).value_counts().unstack(-1)
		# expected maf counts: taken by multiplying the propoprtion of ALS/CTR and total maf counts 
		proportion = maf_counts.sum(axis=1) / maf_counts.sum(axis=1).sum()  # proportion of ALS vs. CTR
		exp_counts = np.outer(proportion, maf_counts.sum(axis=0))

		results = stats.chisquare(maf_counts, f_exp=exp_counts, axis=0)
		return pd.Series(data=results.pvalue, index=maf_counts.columns)

		# return stats.chisquare(maf_counts, f_exp=exp_counts, axis=0)


from .omic import inverse_normal_transform

class Coverage(pd.DataFrame): 

	@property
	def _constructor(self):
		return Coverage

	def normalize(self, norm_factors, inverse_norm_transform=True): 
		df = self / norm_factors
		if inverse_norm_transform: 
			df = inverse_normal_transform(df)
		return df

	def residualize(self, residualizer): 
		self_t = torch.tensor(self.values, dtype=torch.float).to("cpu")
		return pd.DataFrame(residualizer.transform(self_t), index=self.index, columns=self.columns)

	def bin(self, bin_size): 
		df = self.groupby(self.index // bin_size).mean()
		df.index = df.index * bin_size
		return df







