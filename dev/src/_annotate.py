#!/usr/bin/env python3

import pandas as pd
import numpy as np

from src.bed import df_to_pr
from src.utils import flatten


def annotate_regions_uniquely(regions_pr, index_name, annotations, feature_order, fill_na=None): 

	feature_col = annotations.columns.drop(["Chromosome", "Start", "End", "Strand"], errors="ignore")[0]
	annotations = annotations.merge(by=feature_col, strand=False) # Merges overlapping intervals by feature_col and unstrand

	joined = regions_pr.join(annotations, strandedness=None, report_overlap=True).as_df()
	joined[feature_col] = pd.Categorical(joined[feature_col], categories=feature_order, ordered=True)
	joined = joined.sort_values(["Chromosome", "Start", feature_col], ascending=[True, True, True]).drop_duplicates(subset=[index_name]).set_index(index_name)

	index = flatten([ df[index_name] for _,df in regions_pr.items() ]) 
	joined = joined.reindex(index)
	if fill_na is not None: 
		joined[feature_col].fillna(fill_na, inplace=True)

	return joined





def annotation_counts_per_region(regions, annotations, unique=False, feature_order=None, fraction=False): 
	"""Returns counts or fraction of regions that are assigned to a given annotation."""
	regions = regions.drop_duplicate_positions()
	feature_col = annotations.columns.drop(["Chromosome", "Start", "End", "Strand"], errors="ignore")[0]
	annotations = annotations.merge(by=feature_col, strand=False) # Merges overlapping intervals by feature_col and unstrand

	joined = regions.join(annotations, how="left", strandedness=None, report_overlap=True)

	if unique: 
		# TODO: need to break ties which happen when the region overlaps completely
		if feature_order is not None: 
			joined[feature_col] = pd.Categorical(joined[feature_col], categories=feature_order, ordered=True)

		joined = joined.as_df().sort_values(["Chromosome", "Start", feature_col], ascending=[True, True, True]).drop_duplicates(subset=["Chromosome", "Start", "End"])
		counts = joined[feature_col]
	else: 
		counts = joined.as_df().drop_duplicates(annotations.columns)[feature_col]

	if fraction: 
		return counts / len(regions)
	else: 
		return counts



# def annotation_counts_per_region(regions, annotations, unique=False, fraction=False): 
# 	"""Returns counts or fraction of regions that are assigned to a given annotation."""
# 	regions = regions.drop_duplicate_positions()
# 	feature_col = annotations.columns.drop(["Chromosome", "Start", "End", "Strand"], errors="ignore")[0]
# 	annotations = annotations.merge(by=feature_col, strand=False) # Merges overlapping intervals by feature_col and unstrand

# 	joined = regions.join(annotations, strandedness=False, report_overlap=True)

# 	if unique: 
# 		# TODO: need to break ties which happen when the region overlaps completely
# 		# joined[feature_col] = pd.Categorical(joined[feature_col], categories=["promoter", "intron", "exon"], ordered=True)
# 		joined = joined.as_df().sort_values(["Chromosome", "Start", "Overlap", feature_col], ascending=[True, True, False, True]).drop_duplicates(subset=["Chromosome", "Start", "End"])
# 		counts = joined[feature_col].value_counts()
# 	else: 
# 		counts = joined.as_df().drop_duplicates(annotations.columns)[feature_col].value_counts()

# 	if fraction: 
# 		return counts / len(regions)
# 	else: 
# 		return counts


def annotation_enrichments(bg_regions, fg_regions, annotations, unique): 
	"""Enrichment of annotations by region."""
	bg_counts = annotation_counts_per_region(bg_regions, annotations, unique=unique)
	fg_counts = annotation_counts_per_region(fg_regions, annotations, unique=unique)

	ct = pd.DataFrame.from_dict({
		"TP": fg_counts, 
		"FP": len(fg_regions) - fg_counts, 
		"TN": (len(bg_regions) - len(fg_regions)) - (bg_counts - fg_counts), 
		"FN": bg_counts - fg_counts
	})
	
	results = pr.stats.fisher_exact(ct.TP, ct.FP, ct.FN, ct.TN)
	results.index = ct.index
	
	return results







def get_variant_annos_as_df(variants):

	from biothings_client import get_client, alwayslist

	mv = get_client('variant')
	results = mv.getvariants(set(variants))

	cadd_dic = {}
	snpeff_dic = {}
	for result in results: 
		variant = result["query"]
		filtered_result = dict()


		# try: 
		# 	cadd_dic[variant]["annos"] = result["cadd"]["annotype"]
		# except: 
		# 	cadd_dic[variant]["annotype"] = dict()

		# try: 
		# 	snpeff_dic[variant]["ann"] = result["snpeff"]["ann"]
		# except: 
		# 	snpeff_dic[variant]["ann"] = dict()
		try: 
			snpeff_dic[variant] = result["snpeff"]
		except: 
			snpeff_dic = np.nan




		# if result["notfound"] == True: 
		# 	results_dic[variant] = np.nan
		# else: 
		# 	results_dic[variant] = {
		# 		"annotype": result[""]
		# 	}

	return snpeff_dic