#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pysam
from math import ceil
from collections import OrderedDict

from pathos.multiprocessing import ProcessingPool

from src import logger
from src import aals



class QueryBam: 
	"""Interface for querying bams and caching queries."""
	def __init__(self, bam_files, max_cache=20, n_cpus=24): 

		self.bam_files = bam_files
		self.max_cache = max_cache
		self.n_cpus = n_cpus
		self.reset_cache()

		self.sample_names = self.bam_files.index

	def reset_cache(self): 
		self._pileup_cache = list()
		self._coverage_cache = list()

	def _cache_pileup(self, interval, pileup, params): 
		"""Manage cacheing for pileup queries."""

		cache_entry = {"interval": interval.tag, "params": params, "pileup": pileup}
		self._pileup_cache.append(cache_entry)

		n_large = len([None for entry in self._pileup_cache if len(entry["pileup"]) > 1e4])
		if n_large > self.max_cache: 
			self._pileup_cache = self._pileup_cache[1:]

	def _check_pileup_cache(self, interval, params): 
		"""Check that a previous query has been cached."""
		for cache_entry in self._pileup_cache: 
			if interval.tag == cache_entry["interval"] and params == cache_entry["params"]: 
				return cache_entry["pileup"]
		return None

	def query_pileup(self, interval, 
		lib_factors=None, strand=None, 
		relative_to=None, as_df=False, 
		verbose=True, **params): 

		# Initialize query parameters
		chrom, start, end = interval.tuple3
		if "disjoint" not in params: params["disjoint"] = True

		# Check cache for previous results
		pileup = self._check_pileup_cache(interval, params)
		if pileup is not None:
			logger.write("loaded from cache", verbose=verbose)
		else: 
			# Otherwise generate results and cache
			pileup = get_pileups_in_interval(self.bam_files, chrom, start, end, **{**params, **{"n_cpus": self.n_cpus}})
			self._cache_pileup(interval, pileup, params)
			logger.write("wrote to cache", verbose=verbose)

		# Format results
		processed_pileup = process_pileup(pileup, lib_factors=lib_factors, strand=strand)
		if as_df: 
			processed_pileup = self.pileup_array_to_df(processed_pileup, interval, relative_to)
		
		return processed_pileup


	def pileup_array_to_df(self, pileup, interval=None, relative_to=None): 
		"""Converts numpy query result to dataframe."""
		# assert not (invert_index and relative_to is None), "Index should be centered before inverting"

		index = range(interval.start, interval.end)
		if pileup.ndim == 2: 
			pileup_df = pd.DataFrame(pileup, index=index, columns=self.sample_names)
		else: 
			assert pileup.ndim == 3 and pileup.shape[0] == 2
			pileup_df = pd.concat({
				"pos": pd.DataFrame(pileup[0], index=index, columns=self.sample_names), 
				"neg": pd.DataFrame(pileup[1], index=index, columns=self.sample_names), 
			}, axis=1)

		if relative_to: 
			center_pos = getattr(interval, relative_to)
			pileup_df.index = pileup_df.index - center_pos

		# if invert_index: 
		# 	pileup_df.index = pileup_df.index[::-1]
		# 	pileup_df.sort_index(inplace=True)

		return pileup_df



#----------------------------------------------------------------------------------------------------#
# QueryBam helpers
#----------------------------------------------------------------------------------------------------#
def process_pileup(pileup, lib_factors=None, strand=None): 
	
	processed_pileup = pileup.copy()

	# Get stranded pileups
	if strand == "pos": 
		processed_pileup = processed_pileup[0]
	elif strand == "neg": 
		processed_pileup = processed_pileup[1]
	elif strand == "both": 
		processed_pileup = processed_pileup.sum(axis=0)
	else: 
		pass

	# Normalize by library size
	if lib_factors is not None: 
		if isinstance(lib_factors, pd.Series): lib_factors = lib_factors.values
		processed_pileup = processed_pileup / lib_factors

	return processed_pileup

def invert_pileup(pileup): 
	pileup.index = pileup.index[::-1]
	pileup.sort_index(inplace=True)






def load_example_bam(omic="rna", i=0): 
	from src import aals
	bam_file = aals.bam_paths[f"{omic}_bam"][i]
	return open_bam_file(bam_file)

#--------------------------------------------------------------------------------------------------#
def get_depth_per_million(bam_files: pd.Series): 
	results = multiprocess_wrapper(_get_library_size, ((p,) for p in bam_files.values))
	return pd.Series(data=results, index=bam_files.index)


def _get_library_size(bam_file): 
	results = pysam.idxstats(bam_file)
	return np.array([row.split("\t")[2:] for row in results.rstrip("\n").split("\n")], dtype=int).sum() / 1e6





def get_pileups_in_interval(bam_files, chrom, start, end, disjoint=True, n_cpus=24): 
	"""
	Get pileups in interval by strand.
	"""
	args = ((path, chrom, start, end) for path in bam_files)

	query_func = _get_disjoint_pileups_in_interval if disjoint else _get_contiguous_pileups_in_interval
	results = multiprocess_wrapper(query_func, args, n_cpus=n_cpus)

	return np.moveaxis(np.stack(results, axis=1), 2, 0) # merge results (2, n_bases, n_samples)

def _get_disjoint_pileups_in_interval(bam_file, chrom, start, end): 
	"""
	Only mapped read segments are counted in pileups. 
	"""
	pysam_obj = open_bam_file(bam_file)
	chrom = chrom if pysam_chromosome_prefix(pysam_obj) == "chr" else chrom[3:]

	region_len = end - start
	slopes = [[0]*(region_len+1) for _ in range(2)]

	# Iterate through reads that overlap query region
	for read in pysam_obj.fetch(chrom, start, end): 

		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 

			# Count only mappable parts of the read which excludes introns
			for b_start,b_end in read.get_blocks(): 
				if not (b_start >= end or b_end <= start):  # skip reads that fall out of region
					rel_start = b_start - start if b_start >= start else 0
					rel_end = b_end - start if b_end <= end else region_len

					slopes[read.is_reverse][rel_start] += 1
					slopes[read.is_reverse][rel_end] -= 1

	pysam_obj.close()

	pileups = np.array(slopes).cumsum(axis=1)[:,:-1]
	return pileups.T  # (`region_len`, 2)

def _get_contiguous_pileups_in_interval(bam_file, chrom, start, end): 
	"""
	Fills in segments in between mapped portions of reads (i.e. intronic reads that have been spliced)
	"""
	pysam_obj = open_bam_file(bam_file)
	chrom = chrom if pysam_chromosome_prefix(pysam_obj) == "chr" else chrom[3:]

	region_len = end - start
	slopes = [[0]*(region_len+1) for _ in range(2)]

	# Iterate through reads that overlap query region
	for read in pysam_obj.fetch(chrom, start, end): 

		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 

			# Count full length of read including intronic regions
			b_start, b_end = read.reference_start, read.reference_end
			if not (b_start >= end or b_end <= start):  # skip reads that fall out of region
				rel_start = b_start - start if b_start >= start else 0
				rel_end = b_end - start if b_end <= end else region_len

				slopes[read.is_reverse][rel_start] += 1
				slopes[read.is_reverse][rel_end] -= 1

	pysam_obj.close()

	pileups = np.array(slopes).cumsum(axis=1)[:,:-1]
	return pileups.T  # (`region_len`, 2)

#----------------------------------------------------------------------------------------------------#
# Pileups
#----------------------------------------------------------------------------------------------------#
# def get_pileups_in_interval(bam_files, chrom, start, end, n_cpus=24): 
# 	args = ((path, chrom, start, end) for path in bam_files)
# 	results = multiprocess_wrapper(_get_pileup_in_interval, args, n_cpus=n_cpus)
# 	return reshape_bam_results(results, bam_files.index)


# def _get_pileup_in_interval(bam_file, chrom, start, end): 
# 	"""
# 	Pileups for a single bam file in an interval.
# 	 1. Count the number of times a read start or ends at each position. This is essentially the slope.
# 	 2. Count base coverage by taking the cumulative sum (integrate)

# 	# TODO: potentially faster solution using dicts commented at the end. Since most coordinates will have
# 	# no start/end reads, we work with much larger arrays than required. 
# 	"""
# 	pysam_obj = open_bam_file(bam_file)
# 	chrom = chrom if pysam_chromosome_prefix(pysam_obj) == "chr" else chrom[3:]

# 	region_len = end - start
# 	slopes_with_introns = [[0]*(region_len+1) for _ in range(2)]
# 	slopes_no_introns = [[0]*(region_len+1) for _ in range(2)]

# 	filtered_counts = 0

# 	# Iterate through reads that overlap query region
# 	for read in pysam_obj.fetch(chrom, start, end): 
# 		# if read.is_duplicate | read.is_qcfail | read.is_secondary | read.is_unmapped | (not read.is_proper_pair): 
# 		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 

# 			# Count full length of read including intronic regions
# 			b_start, b_end = read.reference_start, read.reference_end
# 			if not (b_start >= end or b_end <= start):  # skip reads that fall out of region
# 				rel_start = b_start - start if b_start >= start else 0
# 				rel_end = b_end - start if b_end <= end else region_len

# 				slopes_with_introns[read.is_reverse][rel_start] += 1
# 				slopes_with_introns[read.is_reverse][rel_end] -= 1

# 			# Count only mappable parts of the read which excludes introns
# 			for b_start,b_end in read.get_blocks(): 
# 				if not (b_start >= end or b_end <= start):  # skip reads that fall out of region
# 					rel_start = b_start - start if b_start >= start else 0
# 					rel_end = b_end - start if b_end <= end else region_len

# 					slopes_no_introns[read.is_reverse][rel_start] += 1
# 					slopes_no_introns[read.is_reverse][rel_end] -= 1

# 	pysam_obj.close()

# 	# Integrate
# 	slopes = slopes_with_introns + slopes_no_introns
# 	pileups = pd.DataFrame(np.array(slopes).cumsum(axis=1)[:,:-1].T, index=range(start,end))
	
# 	pos_pileups = pileups[[0,2]].rename(columns={0:"include_introns", 2:"exclude_introns"})
# 	neg_pileups = pileups[[1,3]].rename(columns={1:"include_introns", 3:"exclude_introns"})

# 	return pos_pileups, neg_pileups



#----------------------------------------------------------------------------------------------------#
# Coverage
#----------------------------------------------------------------------------------------------------#
def get_coverages_in_regions(bam_files, regions, n_cpus=24): 
	"""Coverage across all bams in a set of regions."""
	args = ((path, regions) for path in bam_files)
	results = multiprocess_wrapper(_get_coverage_in_regions, args, n_cpus=n_cpus)
	return reshape_bam_results(results, bam_files.index)


def _get_coverage_in_regions(bam_file, regions): 
	"""Coverage for a single bam file in a set of regions."""
	pysam_obj = open_bam_file(bam_file)

	n_regions = len(regions)
	pos_counts = {"complete_overlap": [0] * n_regions, "partial_overlap": [0] * n_regions}
	neg_counts = {"complete_overlap": [0] * n_regions, "partial_overlap": [0] * n_regions}

	# Get regions generator with correct chrom prefix
	regions = regions.copy()
	if pysam_chromosome_prefix(pysam_obj) == "": 
		regions["chrom"] = regions["chrom"].str[3:]
	regions_iterator = ((i,row["chrom"],row["start"],row["end"]) for i,(_,row) in enumerate(regions.iterrows()))

	for i, chrom, start, end in regions_iterator: 

		complete_overlap = [0,0]
		partial_overlap = [0,0]

		for read in pysam_obj.fetch(chrom, start, end): 
			if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 
				if read.reference_start >= start and read.reference_end <= end: 
					complete_overlap[read.is_reverse] += 1
				else: 
					partial_overlap[read.is_reverse] += 1

		pos_counts["complete_overlap"][i], neg_counts["complete_overlap"][i] = complete_overlap
		pos_counts["partial_overlap"][i], neg_counts["partial_overlap"][i] = partial_overlap

	pysam_obj.close()

	pos_counts_df = pd.DataFrame(pos_counts, index=regions.index)
	neg_counts_df = pd.DataFrame(neg_counts, index=regions.index)

	return pos_counts_df, neg_counts_df

#----------------------------------------------------------------------------------------------------#
# Binned coverage
#----------------------------------------------------------------------------------------------------#
def get_binned_coverages_in_interval(bam_files, chrom, start, end, bin_size, n_cpus=24): 
	args = ((path, chrom, start, end, bin_size) for path in bam_files)
	results = multiprocess_wrapper(_get_binned_coverage_in_interval, args, n_cpus=n_cpus)

	pos_results = pd.concat({guid:result[0] for guid,result in zip(bam_files.index, results)}, axis=1)
	neg_results = pd.concat({guid:result[1] for guid,result in zip(bam_files.index, results)}, axis=1)
	both_results = pos_results + neg_results

	combined_results = pd.concat({"pos": pos_results, "neg": neg_results, "both": both_results}, axis=1)
	return combined_results.reindex(columns=bam_files.index, level=-1)

def _get_binned_coverage_in_interval(bam_file, chrom, start, end, bin_size): 
	
	pysam_obj = open_bam_file(bam_file)
	chrom = chrom if pysam_chromosome_prefix(pysam_obj) == "chr" else chrom[3:]

	n_bins = ceil((end - start) / bin_size)
	# counts = np.zeros(n_bins, dtype=int)
	counts = [[0] * n_bins for _ in range(2)]

	for read in pysam_obj.fetch(chrom, start, end): 
		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 
			# assign read to the 5' end
			if read.is_reverse:
				if read.reference_end <= end: 
					counts[1][(read.reference_end - start) // bin_size] += 1
			else: 
				if read.reference_start >= start: 
					counts[0][(read.reference_start - start) // bin_size] += 1

	pysam_obj.close()

	# Get bin names
	boundaries = list(np.arange(start, end, step=bin_size)) + [end] 
	bin_ids = [f"{chrom}:{b_start}-{b_end}" for b_start,b_end in zip(boundaries[:-1], boundaries[1:])]

	pos_counts = pd.Series(data=counts[0], index=bin_ids)
	neg_counts = pd.Series(data=counts[1], index=bin_ids)

	return pos_counts, neg_counts



#----------------------------------------------------------------------------------------------------#
# Helper functions
#----------------------------------------------------------------------------------------------------#
def open_bam_file(bam_file, threads=2): 
	return pysam.AlignmentFile(bam_file, "rb", threads=threads)

def pysam_chromosome_prefix(pysam_obj): 
	return "chr" if "chr" in pysam_obj.references[0] else "" 

def multiprocess_wrapper(func, args, n_cpus=24): 
	"""Wrapper for pathos"""
	with ProcessingPool(n_cpus) as pool: 
		results = pool.map(lambda args: func(*args), args)
	return results

def reshape_bam_results(results, samples): 
	# Format for all individual runs should be a tuple of pos and neg counts dataframes with the form:
	# 	(count type) x (region feature)
	pos_results = pd.concat({guid:result[0] for guid,result in zip(samples, results)}, axis=1)
	neg_results = pd.concat({guid:result[1] for guid,result in zip(samples, results)}, axis=1)

	pos_results = pos_results.swaplevel(0,1,axis=1).sort_index(axis=1)
	neg_results = neg_results.swaplevel(0,1,axis=1).sort_index(axis=1)

	results_df = pd.concat({"pos":pos_results, "neg":neg_results, "both":pos_results+neg_results}, axis=1)
	results_df = results_df.swaplevel(0,1,axis=1).sort_index(axis=1)
	results_df.columns.set_names(["count_method", "strand", "GUID"], inplace=True)

	logger.write("Level 0: `count_method` - [ {} ]".format(", ".join(results_df.columns.levels[0])))
	logger.write("Level 1: `strand` - [ {} ]".format(", ".join(results_df.columns.levels[1])))
	logger.write("Level 2: `GUID`")
	return results_df.reindex(columns=samples, level=-1)





# def _get_pileup_in_interval_dic(bam_file, chrom, start, end, include_introns=False, **read_filters): 
# 	"""
# 	Pileups for a single bam file in an interval.
# 	 1. Count the number of times a read start or ends at each position. This is essentially the slope.
# 	 2. Count base coverage by taking the cumulative sum (integrate)
# 	"""
# 	pysam_obj = open_bam_file(bam_file)
# 	chrom = chrom if pysam_chromosome_prefix(pysam_obj) == "chr" else chrom[3:]

# 	region_len = end - start
# 	slopes = {}

# 	filtered_counts = 0

# 	# Iterate through reads that overlap query region
# 	for read in pysam_obj.fetch(chrom, start, end): 
# 		# if read.is_duplicate | read.is_qcfail | read.is_secondary | read.is_unmapped | (not read.is_proper_pair): 
# 		if read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary: 
# 			continue

# 		# Count full length of read including intronic regions
# 		b_start, b_end = read.reference_start, read.reference_end
# 		try:
# 			slopes[b_start][read.is_reverse] += 1
# 		except:
# 			slopes[b_start] = [0,0,0,0]
# 		try:
# 			slopes[b_end][read.is_reverse] -= 1
# 		except:
# 			slopes[b_end] = [0,0,0,0]

# 		# Count only mappable parts of the read which excludes introns
# 		for b_start,b_end in read.get_blocks(): 
# 			try:
# 				slopes[b_start][read.is_reverse+2] += 1
# 			except: 
# 				slopes[b_start] = [0,0,0,0]
# 			try:
# 				slopes[b_end][read.is_reverse+2] -= 1
# 			except:
# 				slopes[b_end] = [0,0,0,0]

# 	pos = np.array(list(slopes.keys()))
# 	vals = np.array(list(slopes.values()))

# 	pos_filt_idx = (pos >= start) & (pos <= end)
# 	pos = pos[pos_filt_idx]
# 	vals = vals[pos_filt_idx]

# 	pos_idx = pos.argsort()
# 	sorted_arr1 = pos[pos_idx[::-1]]
# 	counts = vals[pos_idx[::-1]]

# 	counts.cumsum(axis=0)

# 	return counts

# 	# Integrate
# 	pileups = np.array(slopes).cumsum(axis=1)
# 	pos_counts = pd.Series(data=pileups[0,:-1], index=range(start,end))
# 	neg_counts = pd.Series(data=pileups[1,:-1], index=range(start,end))

# 	return {"pos": pos_counts, "neg": neg_counts}










