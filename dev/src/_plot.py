#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pyranges as pr

import pysam

# Plotting
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(font_scale=1.5)
sns.set_style(style='white')





## Tracks

# def _get_base_coverage_from_bam(bam, chrom, start, end, w=0): 

# 	# Check format of chromosomes (i.e. if `chr1` or just `1`)
# 	if "chr" in bam.references[0]: 
# 		if "chr" not in chrom: 
# 			chrom = "chr" + str(chrom)
# 	else: 
# 		if "chr" in chrom: 
# 			chrom = chrom[3:]

# 	positions, pileup = [], []
# 	for pileupcolumn in bam.pileup(chrom, start-w, end+w, truncate=True):
# 		positions.append(pileupcolumn.pos)
# 		pileup.append(pileupcolumn.n)
# 	return pd.Series(data=pileup, index=positions)

# def get_base_coverages_from_bams(bams_dic, chrom, start, end, w=0, size_factors=1): 

# 	pileup = {}
# 	for guid,bam in bams_dic.items(): 
# 		pileup[guid] = _get_base_coverage_from_bam(bam, chrom, start, end, w)
# 	pileup = pd.DataFrame(pileup)

# 	return pileup.div(size_factors)




# from functools import partial
# from pathos.multiprocessing import ProcessingPool as Pool
# # from pathos.pools import ProcessPool

# def _get_base_coverage_from_bam_file(bam_path, chrom, start, end, w=0): 

# 	with pysam.AlignmentFile(bam_path, "rb") as bam_pysam:

# 		# Check format of chromosomes (i.e. if `chr1` or just `1`)
# 		if "chr" in bam_pysam.references[0]: 
# 			if "chr" not in chrom: 
# 				chrom = "chr" + str(chrom)
# 		else: 
# 			if "chr" in chrom: 
# 				chrom = chrom[3:]

# 		positions, pileup = [], []
# 		for pileupcolumn in bam_pysam.pileup(chrom, start-w, end+w, truncate=True):
# 			positions.append(pileupcolumn.pos)
# 			pileup.append(pileupcolumn.n)
# 		return pd.Series(data=pileup, index=positions)







	# from functools import partial
	# query = partial(test, start=start, end=end)
	# # query = partial(get_base_coverage_from_bam, chrom=chrom, start=start, end=end, w=w)

	# from pathos.multiprocessing import ProcessingPool as Pool
	# with Pool(2) as pool: 
	# 	results = pool.map(query, list(bams_dic.values()))
	# 	results = pool.map(query, chrom)
	# 	return results
	# 	# results = pool.map(query, list(bams_dic.values()))

	# return pd.DataFrame(dict(zip(bams_dic.keys(), results)))



