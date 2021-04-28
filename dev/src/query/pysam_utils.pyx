# cython: language_level=3

import numpy as np
cimport numpy as np

from pysam.libchtslib cimport *
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment

import pysam


def _get_pileups_from_interval_with_deletions(str bam_file, str chrom, int start, int end, int threads): 
	"""
	Fast implementation to get pileups from gapped reads (with deletions). 
	"""
	# open and prepare bam file
	cdef AlignmentFile pysam_obj = AlignmentFile(bam_file, "rb", threads=threads)
	cdef AlignedSegment read
	if "chr" not in pysam_obj.references[0]: chrom = chrom[3:]

	cdef int region_len = end - start
	cdef np.ndarray[np.int32_t, ndim=2] slopes = np.zeros((2, region_len+1), dtype=np.int32)

	cdef uint32_t k, pos, l
	cdef int op
	cdef uint32_t * cigar_p
	cdef int rel_start, rel_end, strand

	for read in pysam_obj.fetch(chrom, start, end): 
		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 
			# acess cigar attributes
			src = read._delegate
			cigar_p = pysam_bam_get_cigar(src)
			pos = src.core.pos  # start position of read
			strand = read.is_reverse
			l = 0  # position relative to start of read

			for k from 0 <= k < pysam_get_n_cigar(src):  # iterate through k read blocks

				op = cigar_p[k] & BAM_CIGAR_MASK
				l = cigar_p[k] >> BAM_CIGAR_SHIFT

				if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:

					# position of block relative to start query
					rel_start = pos - start
					rel_end = pos + l - start

					if not (rel_start > region_len or rel_end < 0):  # remove blocks that don't overlap query

						# clip block boundaries
						if rel_start < 0: rel_start = 0
						if rel_end > region_len: rel_end = region_len

						# use block boundaries to set slope based on strand
						slopes[strand, rel_start] += 1
						slopes[strand, rel_end] -= 1
				pos += l
	pysam_obj.close()

	# compute pileup as the cumulative sum of the slopes
	cdef int i
	for i in range(region_len): 
		slopes[0,i+1] += slopes[0,i]
		slopes[1,i+1] += slopes[1,i]

	return slopes[:,:region_len]


def _get_pileups_from_interval_without_deletions(str bam_file, str chrom, int start, int end, int threads): 
	"""
	Fast implementation to get pileups from filled reads (without deletions). 
	"""
	# open and prepare bam file
	cdef AlignmentFile pysam_obj = AlignmentFile(bam_file, "rb", threads=threads)
	cdef AlignedSegment read
	if "chr" not in pysam_obj.references[0]: chrom = chrom[3:]

	cdef int region_len = end - start
	cdef np.ndarray[np.int32_t, ndim=2] slopes = np.zeros((2, region_len+1), dtype=np.int32)

	cdef int rel_start, rel_end, strand

	for read in pysam_obj.fetch(chrom, start, end): 
		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 

			# position of bounds relative to start query
			rel_start = read.reference_start - start
			rel_end = read.reference_end - start
			strand = read.is_reverse

			if not (rel_start > region_len or rel_end < 0):  # remove reads that don't overlap query

				# clip boundaries
				if rel_start < 0: rel_start = 0
				if rel_end > region_len: rel_end = region_len

				# use boundaries to set slope based on strand
				slopes[strand, rel_start] += 1
				slopes[strand, rel_end] -= 1
	pysam_obj.close()

	# compute pileup as the cumulative sum of the slopes
	cdef int i
	for i in range(region_len): 
		slopes[0,i+1] += slopes[0,i]
		slopes[1,i+1] += slopes[1,i]

	return slopes[:,:region_len]

#----------------------------------------------------------------------------------------------------#
# Deprecated
#----------------------------------------------------------------------------------------------------#
def _get_pileups_in_interval(str bam_file, str chrom, int start, int end, int threads): 
	"""
	Get pileups from both filled and gapped reads, making one pass through the data. 
	Not performant. 
	"""
	cdef int region_len = end - start
	cdef np.ndarray[np.int32_t, ndim=2] slopes_fill = np.zeros((2, region_len+1), dtype=np.int32)
	cdef np.ndarray[np.int32_t, ndim=2] slopes_gap = np.zeros((2, region_len+1), dtype=np.int32)

	cdef AlignmentFile pysam_obj = AlignmentFile(bam_file, "rb", threads=threads)
	cdef AlignedSegment read

	if "chr" not in pysam_obj.references[0]: chrom = chrom[3:]

	cdef uint32_t k, pos, l
	cdef int op
	cdef uint32_t * cigar_p
	cdef bam1_t * src
	cdef int rel_start, rel_end, strand
	cdef int rel_min, rel_max

	for read in pysam_obj.fetch(chrom, start, end): 

		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 

			src = read._delegate
			pos = src.core.pos
			strand = read.is_reverse
			cigar_p = pysam_bam_get_cigar(src)
			l = 0

			for k from 0 <= k < pysam_get_n_cigar(src):
				op = cigar_p[k] & BAM_CIGAR_MASK
				l = cigar_p[k] >> BAM_CIGAR_SHIFT
				if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:

					rel_start = pos - start
					rel_end = pos + l - start
					# reads.append([rel_start, rel_end, read.is_reverse])
					if not (rel_start > region_len or rel_end < 0):

						if rel_start < 0: rel_start = 0
						if rel_end > region_len: rel_end = region_len
						slopes_gap[strand, rel_start] += 1
						slopes_gap[strand, rel_end] -= 1
				pos += l

			# position of bounds relative to start query
			rel_start = read.reference_start - start
			rel_end = read.reference_end - start

			if not (rel_start > region_len or rel_end < 0):  # remove reads that don't overlap query

				# clip boundaries
				if rel_start < 0: rel_start = 0
				if rel_end > region_len: rel_end = region_len

				# use boundaries to set slope based on strand
				slopes_fill[strand, rel_start] += 1
				slopes_fill[strand, rel_end] -= 1
	pysam_obj.close()

	cdef int i
	for i in range(region_len): 
		slopes_fill[0,i+1] += slopes_fill[0,i]
		slopes_fill[1,i+1] += slopes_fill[1,i]
		slopes_gap[0,i+1] += slopes_gap[0,i]
		slopes_gap[1,i+1] += slopes_gap[1,i]

	return slopes_fill[:,:region_len], slopes_gap[:,:region_len]


def _get_read_boundaries_in_interval(str bam_file, str chrom, int start, int end): 
	"""
	Tally positions of read boundaries using dictionaries. 
	Not performant. 
	"""
	cdef AlignmentFile pysam_obj = AlignmentFile(bam_file, "rb", threads=3)
	cdef AlignedSegment read

	if "chr" not in pysam_obj.references[0]: chrom = chrom[3:]

	cdef uint32_t k, pos, l
	cdef int op
	cdef uint32_t * cigar_p
	cdef bam1_t * src
	cdef int lower_bound, upper_bound

	pos_bound_counter = {}
	neg_bound_counter = {}

	for read in pysam_obj.fetch(chrom, start, end): 
		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 

			src = read._delegate
			pos = src.core.pos
			cigar_p = pysam_bam_get_cigar(src)
			l = 0

			for k from 0 <= k < pysam_get_n_cigar(src):
				op = cigar_p[k] & BAM_CIGAR_MASK
				l = cigar_p[k] >> BAM_CIGAR_SHIFT
				if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
					lower_bound = pos
					upper_bound = pos + l

					if read.is_reverse: 
						if lower_bound not in neg_bound_counter: 
							neg_bound_counter[lower_bound] = 0
						neg_bound_counter[lower_bound] += 1

						if upper_bound not in neg_bound_counter: 
							neg_bound_counter[upper_bound] = 0
						neg_bound_counter[upper_bound] -= 1

					else: 
						if lower_bound not in pos_bound_counter: 
							pos_bound_counter[lower_bound] = 0
						pos_bound_counter[lower_bound] += 1

						if upper_bound not in pos_bound_counter: 
							pos_bound_counter[upper_bound] = 0
						pos_bound_counter[upper_bound] -= 1

				pos += l
	pysam_obj.close()

	return pos_bound_counter, neg_bound_counter




# def _get_pileups_from_interval(str bam_file, str chrom, int start, int end, int threads): 

# 	cdef int region_len = end - start
# 	cdef np.ndarray[np.int32_t, ndim=2] slopes = np.zeros((2, region_len+1), dtype=np.int32)

# 	cdef AlignmentFile pysam_obj = AlignmentFile(bam_file, "rb", threads=threads)
# 	cdef AlignedSegment read

# 	if "chr" not in pysam_obj.references[0]: chrom = chrom[3:]

# 	cdef uint32_t k, pos, l
# 	cdef int op
# 	cdef uint32_t * cigar_p
# 	cdef bam1_t * src
# 	cdef int rel_start, rel_end, strand
# 	cdef int rel_min, rel_max

# 	cdef int counts = 0

# 	# reads = []

# 	for read in pysam_obj.fetch(chrom, start, end): 
# 		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 

# 			src = read._delegate
# 			pos = src.core.pos
# 			cigar_p = pysam_bam_get_cigar(src)
# 			l = 0

# 			for k from 0 <= k < pysam_get_n_cigar(src):
# 				op = cigar_p[k] & BAM_CIGAR_MASK
# 				l = cigar_p[k] >> BAM_CIGAR_SHIFT
# 				if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:

# 					rel_start = pos - start
# 					rel_end = pos + l - start
# 					# reads.append([rel_start, rel_end, read.is_reverse])
# 					if not (rel_start > region_len or rel_end < 0):

# 						if rel_start < 0: rel_start = 0
# 						if rel_end > region_len: rel_end = region_len

# 						# strand = read.is_reverse
# 						# strands.append(read.is_reverse)
# 						# # print(strand, pos, pos+l, rel_start, rel_end)
# 						# counts += 1
# 						strand = read.is_reverse
# 						# print(strand, rel_start, rel_end, region_len)
# 						slopes[strand, rel_start] += 1
# 						slopes[strand, rel_end] -= 1

# 				pos += l

# 	# return reads

# 	# return strands, counts
# 	cdef int i
# 	for i in range(region_len): 
# 		slopes[0,i+1] += slopes[0,i]
# 		slopes[1,i+1] += slopes[1,i]

# 	return slopes[:,:region_len]


# def _get_pileups_from_interval_filled(str bam_file, str chrom, int start, int end, int threads): 

# 	cdef int region_len = end - start
# 	cdef np.ndarray[np.int32_t, ndim=2] slopes = np.zeros((2, region_len+1), dtype=np.int32)

# 	cdef AlignmentFile pysam_obj = AlignmentFile(bam_file, "rb", threads=threads)
# 	cdef AlignedSegment read

# 	if "chr" not in pysam_obj.references[0]: chrom = chrom[3:]

# 	cdef uint32_t k, start_pos, end_pos
# 	cdef int op
# 	cdef uint32_t * cigar_p
# 	cdef bam1_t * src
# 	cdef int rel_start, rel_end, strand

# 	for read in pysam_obj.fetch(chrom, start, end): 

# 		if not (read.is_unmapped or read.is_qcfail or read.is_duplicate or read.is_secondary): 

# 			# rel_start = read.reference_start - start
# 			# rel_end = read.reference_end - start

# 			rel_start = read.reference_start - start
# 			rel_end = read.reference_end - start
			
# 			if not (rel_start > region_len or rel_end < 0):

# 				if rel_start < 0: rel_start = 0
# 				if rel_end > region_len: rel_end = region_len

# 				# strand = read.is_reverse
# 				# strands.append(read.is_reverse)
# 				# # print(strand, pos, pos+l, rel_start, rel_end)
# 				# counts += 1
# 				strand = read.is_reverse
# 				# print(strand, rel_start, rel_end, region_len)
# 				slopes[strand, rel_start] += 1
# 				slopes[strand, rel_end] -= 1

# 	cdef int i
# 	for i in range(region_len): 
# 		slopes[0,i+1] += slopes[0,i]
# 		slopes[1,i+1] += slopes[1,i]

# 	return slopes[:,:-1]




