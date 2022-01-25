#!/usr/bin/env python3

import math
from pathos import multiprocessing

import numpy as np
import pandas as pd
import pysam

from src import logger
from src.query import pysam_utils


class QueryBams: 
	"""Interface for querying bams and caching queries."""

	def __init__(self, omic, bam_paths, max_cache=20, n_cpus=24, n_threads=3): 

		self.omic = omic
		self._bam_paths = bam_paths
		self._max_cache = max_cache
		self._pysam_kwargs = dict(n_cpus=n_cpus, n_threads=n_threads)

		self.sample_names = self._bam_paths.index

		self._validate()
		self.reset_cache()

	def _validate(self): 

		from src.load import aals
		assert self.sample_names.equals(aals.sample_names)

	#----------------------------------------------------------------------------------------------------#
	# Manage cache
	#----------------------------------------------------------------------------------------------------#
	def reset_cache(self): 
		self._pileup_cache = list()
		self._coverage_cache = list()

	def _pileup_query_params(self, omic, interval, fill_deletions): 
		query_params = {"omic": omic, "fill_deletions": fill_deletions}
		return {**query_params, **interval.unstrand().as_dict()}

	def _access_pileup_cache(self, query_params): 
		"""
		Check that a previous query has been cached. Returns subsetted array, otherwise None.
		"""
		if self._max_cache == 0: return None
		# Compare query to every entry in the cache 
		for cached_pileupdata_obj in self._pileup_cache: 
			pileupdata_obj = cached_pileupdata_obj._compare_params_with(query_params)
			if pileupdata_obj: 
				logger.write(f"Loading cached {query_params['omic'].upper()} pileup data... ")
				return pileupdata_obj

		logger.write(f"No cached matches found. Generating {query_params['omic'].upper()} pileup data... ")
		return None

	def _update_pileup_cache(self, raw_pileups, query_params): 
		"""Maintains cache."""

		pileupdata_obj = PileupData(raw_pileups, query_params)
		if self._max_cache == 0: return pileupdata_obj

		self._pileup_cache.append(pileupdata_obj)

		# TODO: Prune by removing those in subset

		# Count large pileup queries and pop the first large one
		large_entry_idx = [i for i,entry in enumerate(self._pileup_cache) if entry.n_pos > 1e4]
		if len(large_entry_idx) > self._max_cache: 
			del self._pileup_cache[large_entry_idx[0]]  # remove the oldest large entry

		return pileupdata_obj

	def display_cache(self): 
		for entry in self._pileup_cache: 
			print(entry)

	#----------------------------------------------------------------------------------------------------#
	# Pileup queries
	#----------------------------------------------------------------------------------------------------#
	def query_raw_pileups(self, interval, fill_deletions): 
		"""Returns PileupData instance"""
		chrom, start, end = interval.as_tuple3()

		# Check for previously queried raw pileups
		query_params = self._pileup_query_params(self.omic, interval, fill_deletions)
		pileupdata_obj = self._access_pileup_cache(query_params)

		# If no matches exist, make query
		if not pileupdata_obj: 
			raw_pileups = get_pileups_in_interval(self._bam_paths, chrom, start, end, fill_deletions, **self._pysam_kwargs)
			pileupdata_obj = self._update_pileup_cache(raw_pileups, query_params)

		return pileupdata_obj

	def get_pileup_data(self, interval, fill_deletions, norm_kwargs=dict(), view_kwargs=dict()): 
		"""Ad hoc queries that return pileup data directly"""

		chrom, start, end = interval.as_tuple3()

		# Check for previously queried raw pileups
		query_params = self._pileup_query_params(self.omic, interval, fill_deletions)
		pileupdata_obj = self._access_pileup_cache(query_params)

		# If no matches exist, make query
		if not pileupdata_obj: 
			raw_pileups = get_pileups_in_interval(self._bam_paths, chrom, start, end, fill_deletions, **self._pysam_kwargs)
			pileupdata_obj = self._update_pileup_cache(raw_pileups, query_params)

		# Bin and nomalize data
		logger.write(f"Creating view...")
		view_kwargs["sample_names"] = self.sample_names
		data = pileupdata_obj.process(view_kwargs, norm_kwargs)
		return data



class PileupData: 

	def __init__(self, raw, params):

		self.raw = raw
		self.params = params

		self._validate()

		self.omic = self.params["omic"]
		self.chrom, self.start, self.end = self.params["chrom"], self.params["start"], self.params["end"]
		self.n_pos = self.end - self.start
		self.positions = np.arange(self.start, self.end)

		# self.normalized, self._lib_factors, self._residualizer = None, None, None

	def _validate(self): 

		for param in ["chrom", "start", "end", "omic", "fill_deletions"]: assert param in self.params
		assert self.raw.shape[1] == self.params["end"] - self.params["start"]

	def _compare_params_with(self, other_params): 
		"""Compares query params to see if there's a match. Coordinates may be a subset of previous query."""

		# These need to be exactly the same
		for param in ["omic", "fill_deletions", "chrom"]: 
			if other_params[param] != self.params[param]:
				return False

		# Queries that are subsets are allowable
		if other_params["start"] < self.params["start"] or other_params["end"] > self.params["end"]: 
			return False

		# Warn about possible other mismatches
		if set(other_params.keys()) != set(self.params.keys()): 
			logger.write("Warning: possible param mismatch.")
			print(other_params)
			print(self.params)

		# Return subsetted data
		start_idx, end_idx = other_params["start"] - self.start, other_params["end"] - self.start
		return PileupData(self.raw[:,start_idx:end_idx,:], other_params)

	def __repr__(self): 
		return f"{self.omic.upper()} pileup ({'filled' if self.params['fill_deletions'] else 'gapped'}) at {self.chrom}:{self.start}-{self.end} ({self.n_pos}bp)"

	#----------------------------------------------------------------------------------------------------#
	# Queries of raw pileup data
	#----------------------------------------------------------------------------------------------------#
	def get_raw_pileups_by_coords(self, start, end):
		"""
		Return self if query is the same, otherwise return a view of a subset. 
		"""
		start_idx, end_idx = start - self.start, end - self.start
		return self.raw[:,start_idx:end_idx,:]

	def process(self, view_kwargs, norm_kwargs): 

		# Get data by strand
		strand = view_kwargs.pop("strand", None)

		if strand is None: 
			data = self.raw.sum(axis=0)
		elif strand == "pos": 
			data = self.raw[0]
		elif strand == "neg": 
			data = self.raw[1]
		else:
			raise ValueError

		# Bin data
		if view_kwargs.pop("downsample", True):
			data, bin_positions = bin_pileup_data(data, self.positions, **view_kwargs)
		else: 
			bin_positions = self.positions.copy()

		# Normalize binned data
		data = normalize_pileup_data(data, **norm_kwargs)

		# Return specified format
		if view_kwargs.pop("as_df", True): 

			sample_names = view_kwargs.pop("sample_names", None)
			zero_pos = view_kwargs.pop("zero_pos", 0)

			if zero_pos: bin_positions -= zero_pos

			return pd.DataFrame(data, index=bin_positions, columns=sample_names)
		else: 
			return data


def bin_pileup_data(data, positions, n_bins=None, max_bins=1000, **kwargs): 

	n_pos, n_samples = data.shape[0], data.shape[1]

	if n_pos < max_bins: return data, positions

	if n_bins is None:
		n_bins = max_bins
		while n_pos / n_bins > max_bins: n_bins -= 1

	bin_edges = np.linspace(positions[0], positions[-1], num=n_bins+1)
	bins = np.digitize(positions, bin_edges) - 1

	# Get mean of values in each bin
	binned_data = np.zeros((n_bins, n_samples))
	bin_positions = np.zeros((n_bins))
	for i in range(n_bins): 
		binned_data[i] = data[bins == i].mean(axis=0)
		bin_positions[i] = positions[bins == i].mean()

	return binned_data, bin_positions.astype(int)

def normalize_pileup_data(data, lib_factors=None, residualizer=None): 

	if lib_factors is not None: 
		data = data / lib_factors.values

	if residualizer is not None: 
		data = residualizer.transform(data)

	return data




#--------------------------------------------------------------------------------------------------#
def get_depth_per_million(bam_files: pd.Series): 
	results = multiprocess_wrapper(_get_library_size, ((p,) for p in bam_files.values))
	return pd.Series(data=results, index=bam_files.index)

def _get_library_size(bam_file): 
	results = pysam.idxstats(bam_file)
	return np.array([row.split("\t")[2:] for row in results.rstrip("\n").split("\n")], dtype=int).sum() / 1e6


#----------------------------------------------------------------------------------------------------#
# Pileups
#----------------------------------------------------------------------------------------------------#
def get_pileups_in_interval(bam_files, chrom, start, end, fill_deletions, n_cpus=24, n_threads=3):
	"""
	Get pileups in interval by strand.
	"""
	assert fill_deletions is not None
	args = ((path, chrom, start, end, n_threads) for path in bam_files)

	func = (pysam_utils._get_pileups_from_interval_with_deletions if fill_deletions 
			else pysam_utils._get_pileups_from_interval_without_deletions)

	results = multiprocess_wrapper(func, args, n_cpus=n_cpus)
	results_np = np.zeros((*results[0].shape, len(results)), dtype=int)
	for i,result in enumerate(results): 
		results_np[:,:,i] = result
	return results_np

	# if return_params: 
	# 	query_params = {"chrom": chrom, "start": start, "end": end, "fill_detections": True}
	# 	return results, query_params
	# else: 
	# 	return results

# def get_pileups_in_interval(bam_files, chrom, start, end, disjoint=True, n_cpus=24): 
# 	"""
# 	Get pileups in interval by strand.
# 	"""
# 	args = ((path, chrom, start, end) for path in bam_files)

# 	query_func = _get_disjoint_pileups_in_interval if disjoint else _get_contiguous_pileups_in_interval
# 	results = multiprocess_wrapper(query_func, args, n_cpus=n_cpus)
# 	return results
# 	return np.moveaxis(np.stack(results, axis=1), 2, 0) # merge results (2, n_bases, n_samples)

# def _get_disjoint_pileups_in_interval(bam_file, chrom, start, end): 
# 	"""
# 	Only mapped read segments are counted in pileups. 
# 	"""
# 	reads = []
# 	count = 0

# 	pysam_obj = open_bam_file(bam_file)
# 	chrom = chrom if pysam_chromosome_prefix(pysam_obj) == "chr" else chrom[3:]

# 	region_len = end - start
# 	slopes = [[0]*(region_len+1) for _ in range(2)]

# 	# Iterate through reads that overlap query region
# 	for read in pysam_obj.fetch(chrom, start, end): 

# 		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 

# 			# Count only mappable parts of the read which excludes introns
# 			for b_start,b_end in read.get_blocks(): 
# 				if not (b_start >= end or b_end <= start):  # skip reads that fall out of region
# 					rel_start = b_start - start if b_start >= start else 0
# 					rel_end = b_end - start if b_end <= end else region_len

# 					slopes[read.is_reverse][rel_start] += 1
# 					slopes[read.is_reverse][rel_end] -= 1
# 					count += 1
# 				# reads.append([b_start - start, b_end - start, read.is_reverse])

# 	pysam_obj.close()

# 	# return reads
# 	# return count

# 	pileups = np.array(slopes).cumsum(axis=1)[:,:-1]
# 	return pileups.T  # (`region_len`, 2)

# def _get_contiguous_pileups_in_interval(bam_file, chrom, start, end): 
# 	"""
# 	Fills in segments in between mapped portions of reads (i.e. intronic reads that have been spliced)
# 	"""
# 	pysam_obj = open_bam_file(bam_file)
# 	chrom = chrom if pysam_chromosome_prefix(pysam_obj) == "chr" else chrom[3:]

# 	region_len = end - start
# 	slopes = [[0]*(region_len+1) for _ in range(2)]

# 	# Iterate through reads that overlap query region
# 	for read in pysam_obj.fetch(chrom, start, end): 

# 		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 

# 			# Count full length of read including intronic regions
# 			b_start, b_end = read.reference_start, read.reference_end
# 			if not (b_start >= end or b_end <= start):  # skip reads that fall out of region
# 				rel_start = b_start - start if b_start >= start else 0
# 				rel_end = b_end - start if b_end <= end else region_len

# 				slopes[read.is_reverse][rel_start] += 1
# 				slopes[read.is_reverse][rel_end] -= 1

# 	pysam_obj.close()

# 	pileups = np.array(slopes).cumsum(axis=1)[:,:-1]
# 	return pileups.T  # (`region_len`, 2)

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

	n_bins = math.ceil((end - start) / bin_size)
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
def open_bam_file(bam_file, n_threads=2): 
	return pysam.AlignmentFile(bam_file, "rb", threads=n_threads)

def load_example_bam(omic="rna", i=0): 
	from src import aals
	bam_file = aals.bam_paths[f"{omic}_bam"][i]
	return open_bam_file(bam_file)

def pysam_chromosome_prefix(pysam_obj): 
	return "chr" if "chr" in pysam_obj.references[0] else "" 

def multiprocess_wrapper(func, args, n_cpus=24): 
	"""Wrapper for pathos"""
	with multiprocessing.ProcessingPool(n_cpus) as pool: 
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










