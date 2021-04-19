#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pysam

from pathos.multiprocessing import ProcessingPool

from src import logger





#--------------------------------------------------------------------------------------------------#
def open_bam_file(bam_file, threads=2): 
	return pysam.AlignmentFile(bam_file, "rb", threads=threads)

def match_chromosome_prefix(pysam_obj, chrom): 
	"""Returns the correctly formatted chromosome ID corresponding to a bam file (i.e. `chr1` vs. `1`)."""
	if "chr" in pysam_obj.references[0]: 
		if "chr" not in chrom: 
			chrom = "chr" + str(chrom)
	else: 
		if "chr" in chrom: 
			chrom = chrom[3:]
	return chrom

def pysam_chromosome_prefix(pysam_obj): 
	return "chr" if "chr" in pysam_obj.references[0] else "" 

def get_regions_generator(pysam_obj, regions): 
	if pysam_chromosome_prefix(pysam_obj) == "":
		regions["chrom"] = regions["chrom"].str[3:]
	return ((i,row["chrom"],row["start"],row["end"]) for i,(index,row) in enumerate(regions.iterrows()))

def multiprocess(func, args, n_cpus=24): 
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
	return results_df

#--------------------------------------------------------------------------------------------------#
def _get_coverage_in_regions(bam_file, regions, read_overlap=None): 
	"""Coverage for a single bam file in a set of regions."""
	pysam_obj = open_bam_file(bam_file)

	n_regions = len(regions)
	pos_counts = {"complete_overlap": [0] * n_regions, "partial_overlap": [0] * n_regions}
	neg_counts = {"complete_overlap": [0] * n_regions, "partial_overlap": [0] * n_regions}

	for i, chrom, start, end in get_regions_generator(pysam_obj, regions): 

		complete_overlap = [0,0]
		partial_overlap = [0,0]

		for read in pysam_obj.fetch(chrom, start, end): 
			if (read.reference_start >= start) & (read.reference_end <= end): 
				complete_overlap[read.is_reverse] += 1
			else: 
				partial_overlap[read.is_reverse] += 1

		pos_counts["complete_overlap"][i], neg_counts["complete_overlap"][i] = complete_overlap
		pos_counts["partial_overlap"][i], neg_counts["partial_overlap"][i] = partial_overlap

	pysam_obj.close()

	pos_counts_df = pd.DataFrame(pos_counts, index=regions.index)
	neg_counts_df = pd.DataFrame(neg_counts, index=regions.index)

	return pos_counts_df, neg_counts_df


def get_coverages_in_regions(bam_files, regions, read_overlap=None): 
	"""Coverage across all bams in a set of regions."""
	args = ((path, regions) for path in bam_files)
	results = multiprocess(_get_coverage_in_regions, args)
	return reshape_bam_results(results, bam_files.index)


#--------------------------------------------------------------------------------------------------#
def _get_binned_coverage_in_interval(bam_file, interval, bin_size): 
	pass


def get_binned_coverages_in_interval(bam_files, interval, bin_size): 
	pass



#--------------------------------------------------------------------------------------------------#

def _get_pileup_in_interval(bam_file, chrom, start, end, **read_filters): 
	"""
	Pileups for a single bam file in an interval.
	 1. Count the number of times a read start or ends at each position. This is essentially the slope.
	 2. Count base coverage by taking the cumulative sum (integrate)
	"""
	pysam_obj = open_bam_file(bam_file)
	chrom = chrom if pysam_chromosome_prefix(pysam_obj) == "chr" else chrom[3:]

	region_len = end - start
	slopes_with_introns = [[0]*(region_len+1) for _ in range(2)]
	slopes_no_introns = [[0]*(region_len+1) for _ in range(2)]

	filtered_counts = 0

	# Iterate through reads that overlap query region
	for read in pysam_obj.fetch(chrom, start, end): 
		# if read.is_duplicate | read.is_qcfail | read.is_secondary | read.is_unmapped | (not read.is_proper_pair): 
		if read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary: 
			continue

		# Count full length of read including intronic regions
		b_start, b_end = read.reference_start, read.reference_end
		if (b_start >= end) or (b_end <= start):
			continue  # skip reads that fall out of region
		else: 
			rel_start = b_start - start if b_start >= start else 0
			rel_end = b_end - start if b_end <= end else region_len

			slopes_with_introns[read.is_reverse][rel_start] += 1
			slopes_with_introns[read.is_reverse][rel_end] -= 1

		# Count only mappable parts of the read which excludes introns
		for b_start,b_end in read.get_blocks(): 
			if (b_start >= end) or (b_end <= start):
				continue
			else:
				rel_start = b_start - start if b_start >= start else 0
				rel_end = b_end - start if b_end <= end else region_len

				slopes_no_introns[read.is_reverse][rel_start] += 1
				slopes_no_introns[read.is_reverse][rel_end] -= 1

	pysam_obj.close()

	# Integrate
	slopes = slopes_with_introns + slopes_no_introns
	pileups = pd.DataFrame(np.array(slopes).cumsum(axis=1)[:,:-1].T, index=range(start,end))
	
	pos_pileups = pileups[[0,2]].rename(columns={0:"include_introns", 2:"exclude_introns"})
	neg_pileups = pileups[[1,3]].rename(columns={1:"include_introns", 3:"exclude_introns"})

	return pos_pileups, neg_pileups


def get_pileups_in_interval(bam_files, chrom, start, end, **read_filters): 
	args = ((path, chrom, start, end) for path in bam_files)
	results = multiprocess(_get_pileup_in_interval, args)
	# return results
	return reshape_bam_results(results, bam_files.index)


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









