#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pyranges as pr

from src import base_dir, logger
from src.load import aals, data, hg38
from src.qtl import preprocess, residualize, run_tensorqtl
from src.utils.regions import Regions


_omic_tag = {"atac": "atac", "rna": "rna", "erna": "rna", "splice": "rna"}  
_qtl_tag = {"atac": "caQTL", "rna": "eQTL", "erna": "erQTL", "splice": "sQTL"}  
_phen_tag = {"atac": "peak_id", "rna": "gene_id", "erna": "peak_id", "splice": "splice_id"}
_dist_tag = {"atac": "peak_dist", "rna": "tss_dist", "erna": "peak_dist", "splice": "splice_dist"}

class FeatureData: 
	"""Accessor for feature counts data and QTLs."""

	def __init__(self, feature, counts_prefix=None, counts_dir=None, qtl_dir=None, qtl_fdr=0.05): 

		self._validate(feature)

		self._feature = feature
		self._omic = _omic_tag[self._feature]
		self._qtl = _qtl_tag[self._feature]
		self._phen_id = _phen_tag[self._feature]
		self._dist_id = _dist_tag[self._feature]

		if counts_prefix: self.initialize_counts(counts_prefix, counts_dir)
		if qtl_dir: self.initialize_qtl(qtl_dir, qtl_fdr)

	def _validate(self, feature): 
		assert feature in _omic_tag

	def initialize_counts(self, counts_prefix, counts_dir):

		if counts_dir is None: counts_dir = 'phenotypes_210422' 

		# Set counts data paths
		self._counts_prefix = counts_prefix
		self._counts_data_paths = {
			"metadata": base_dir / f"tensorqtl_runs/{counts_dir}/{counts_prefix}.metadata.txt.gz", 
			"counts": base_dir / f"tensorqtl_runs/{counts_dir}/{counts_prefix}.counts.txt.gz", 
			"tpm": base_dir / f"tensorqtl_runs/{counts_dir}/{counts_prefix}.tpm.txt.gz", 
			"tmm": base_dir / f"tensorqtl_runs/{counts_dir}/{counts_prefix}.tmm.txt.gz", 
			"gtex": base_dir / f"tensorqtl_runs/{counts_dir}/{counts_prefix}.gtex.txt.gz", 
			"tmm_factors": base_dir / f"tensorqtl_runs/{counts_dir}/{counts_prefix}.tmm_factors.txt.gz", 
		}
		# for key,path in self._counts_data_paths.items(): 
		# 	assert path.is_file(), f"{path} does not exist"

		# Initialize counts property attributes
		self._regions, self._lengths = [None] * 2
		self._counts, self._tpm, self._tmm, self._gtex = [None] * 4
		self._mask = None
		self._tmm_factors = None

	def initialize_qtl(self, qtl_dir, qtl_fdr): 

		# Set QTL data paths
		self._qtl_dir, self._qtl_fdr = qtl_dir, qtl_fdr
		self._qtl_config_path = base_dir / f"tensorqtl_runs/{qtl_dir}/config.txt"
		assert self._qtl_config_path.is_file(), f"{self._qtl_config_path} does not exist"

		# Initialize QTL property attributes
		self._qtl_accessor = run_tensorqtl.TensorQTLManager.load_config(self._qtl_config_path)
		self._qtl_sig, self._qtl_leads, self._qtl_all = [None] * 3
		self._qtl_sig_phen, self._qtl_sig_rsid, self._qtl_all_phen = [None] * 3

		self._covariates_path = self._qtl_accessor.covariates_path
		assert self._covariates_path.is_file(), f"{self._covariates_path} does not exist"
		self._residualizer = None


	#----------------------------------------------------------------------------------------------------#
	# Access counts data
	#----------------------------------------------------------------------------------------------------#
	@property
	def regions(self):
		if self._regions is None: 
			self._regions, self._lengths = preprocess.load_counts_metadata(self._counts_data_paths["metadata"])
			self._regions.index.name = self._phen_id
		return self._regions

	@property
	def lengths(self):
		if self._lengths is None: self._regions, self._lengths = preprocess.load_counts_metadata(self._counts_data_paths["metadata"])
		return self._lengths
	
	@property
	def counts(self):
		if self._counts is None: self._counts = pd.read_csv(self._counts_data_paths["counts"], sep="\t", index_col=0)
		return self._counts

	@property
	def tpm(self):
		if self._tpm is None: self._tpm = pd.read_csv(self._counts_data_paths["tpm"], sep="\t", index_col=0)
		return self._tpm

	@property
	def tmm(self):
		if self._tmm is None: self._tmm = pd.read_csv(self._counts_data_paths["tmm"], sep="\t", index_col=0)
		return self._tmm

	@property
	def gtex(self):
		if self._gtex is None: self._gtex = pd.read_csv(self._counts_data_paths["gtex"], sep="\t", index_col=0)
		return self._gtex

	@property
	def mask(self):
		if self._mask is None: 
			self._mask = self.counts.index.isin(self.gtex.index)
		return self._mask
	
	@property
	def tmm_factors(self):
		if self._tmm_factors is None: self._tmm_factors = pd.read_csv(self._counts_data_paths["tmm_factors"], sep="\t", index_col=0)["tmm"]
		return self._tmm_factors
	
	@property
	def residualizer(self):
		if self._residualizer is None: self._residualizer = residualize.Residualizer.load_covariates(path=self._covariates_path)
		return self._residualizer

	def residualize(self, df): 
		return self.residualizer.transform(df)

	#----------------------------------------------------------------------------------------------------#
	# Access QTL results
	#----------------------------------------------------------------------------------------------------#
	@property
	def qtl_sig(self):
		if self._qtl_sig is None: 
			self._qtl_sig = self._qtl_accessor.load_cis_nominal_results(self._qtl_fdr)
			self._reformat_qtls(self._qtl_sig)
		return self._qtl_sig

	@property
	def qtl_sig_phen(self):
		if self._qtl_sig_phen is None: 
			self._qtl_sig_phen = set(self.qtl_sig[self._phen_id])
		return self._qtl_sig_phen

	@property
	def qtl_sig_rsid(self):
		if self._qtl_sig_rsid is None: 
			self._qtl_sig_rsid = set(self.qtl_sig["variant_id"].unique())
		return self._qtl_sig_rsid

	@property
	def qtl_leads(self):
		if self._qtl_leads is None: 
			self._qtl_leads = self._qtl_accessor.load_cis_results(self._qtl_fdr)
			self._reformat_qtls(self._qtl_leads)
		return self._qtl_leads

	@property
	def qtl_all(self):
		if self._qtl_all is None: 
			self._qtl_all = self._qtl_accessor.load_cis_nominal_results_unfiltered()
			self._reformat_qtls(self._qtl_all)
		return self._qtl_all

	@property
	def qtl_all_phen(self):
		if self._qtl_all_phen is None: 
			self._qtl_all_phen = set(self.qtl_leads.index)
		return self._qtl_all_phen

	
	def _reformat_qtls(self, qtl_results): 
		"""Reformat column names to reflect feature."""

		qtl_results.rename(columns={"phenotype_id": self._phen_id, "tss_distance": self._dist_id}, inplace=True)
		if "gene_id" in qtl_results.columns: 
			# qtl_results["symbol"] = get_gene_symbols(qtl_results.gene_id)
			qtl_results["symbol"] = fops.gene_to_symbol(qtl_results.gene_id)
		if "splice_id" in qtl_results.columns: 
			# qtl_results["cluster_id"] = get_cluster_id(qtl_results.splice_id)
			qtl_results["cluster_id"] = fops.splice_to_cluster(qtl_results.splice_id)

		if qtl_results.index.name == "phenotype_id": 
			qtl_results.index.name = self._phen_id
			if qtl_results.index.name == "gene_id": 
				qtl_results["symbol"] = fops.gene_to_symbol(qtl_results.index)
			if qtl_results.index.name == "splice_id": 
				qtl_results["cluster_id"] = fops.splice_to_cluster(qtl_results.index)

		# if self._feature == "splice": 
		# 	qtl_results["splice_dist"] = correct_splice_dist(qtl_results, self.regions).values

	def print_qtl_summary(self): 

		logger.write(f"---- {self._qtl} at FDR={self._qtl_fdr} ----")
		n_qtl_sig = self.qtl_sig.shape[0]
		n_qtl_total = self.qtl_leads["num_var"].sum()
		frac_qtl = n_qtl_sig / n_qtl_total * 100

		n_phen_sig = self.qtl_sig[self._phen_id].nunique()
		n_phen_total = self.qtl_leads.index.nunique()
		frac_phen = n_phen_sig / n_phen_total * 100

		n_rsid_sig = self.qtl_sig['variant_id'].nunique()

		logger.write(f"{'# QTLs:'.ljust(20)}  {n_qtl_sig:>8,d} / {n_qtl_total:>12,d}  ({frac_qtl:5.2f}%)")
		logger.write(f"{'# unique phenotypes:'.ljust(20)}  {n_phen_sig:>8,d} / {n_phen_total:>12,d}  ({frac_phen:5.2f}%)")
		logger.write(f"{'# unique variants:'.ljust(20)}  {n_rsid_sig:>8,d}")


#----------------------------------------------------------------------------------------------------#
# Feature operations
#----------------------------------------------------------------------------------------------------#
def get_feature_type(s): 
	if s[0].startswith("ENSG"): 
		return "gene_id"
	elif "clu" in s[0]: 
		return "splice_id"
	elif s[0].startswith("rs"): 
		return "variant_id"
	elif s[0].startswith("chr"):
		return "peak_id"
	elif s[0] in data.ensg: 
		return "symbol"
	else: 
		assert False

class FeatureOps: 

	def __init__(self):
		self._splice_id_to_genes, self._splice_id_to_obs_genes = None, None
		pass
		# self.gene_regions, self.peak_regions = None, None

	@property
	def genes(self):
		return self.gene_regions.index

	@property
	def peaks(self):
		return self.peak_regions.index

	@property
	def splice_sites(self):
		return self.splice_regions.index
	
	@property
	def splice_id_to_genes(self):
		if self._splice_id_to_genes is None: 
			self._splice_id_to_genes = _map_splice_id_to_genes(self.splice_regions)
		return self._splice_id_to_genes

	@property
	def splice_id_to_obs_genes(self):
		if self._splice_id_to_obs_genes is None: 
			self._splice_id_to_obs_genes = self.splice_id_to_genes[self.splice_id_to_genes.gene_id.isin(self.genes)]
		return self._splice_id_to_obs_genes

	#----------------------------------------------------------------------------------------------------#
	# Direct namespace mappings
	#----------------------------------------------------------------------------------------------------#
	def gene_to_symbol(self, s): 
		return _gene_id_to_symbol(s, data.ensg)

	def splice_to_cluster(self, s): 
		return _splice_id_to_cluster_id(s)

	def splice_to_gene(self, s): 
		return s.map(self.splice_id_to_obs_genes.drop_duplicates("splice_id")["gene_id"])

	def splice_to_symbol(self, s): 
		return self.gene2symbol(self.splice2gene(s))


	def get_pos(self, s, ref=None):
		feat = get_feature_type(s) 
		if feat == "gene_id": 
			return s.map(self.gene_regions.get_pos(ref))
		elif feat == "peak_id": 
			return s.map(self.peak_regions.get_pos(ref))
		elif feat == "splice_id": 
			return s.map(self.splice_regions.get_pos(ref))
		elif feat == "variant_id": 
			return s.map(data.rsid.pos)

def _gene_id_to_symbol(s, ensg_to_symbol): 
	return s.map(ensg_to_symbol)

def _splice_id_to_cluster_id(s): 
	return s.str.split(":").str[3].str.split("_").str[1].astype(int)


def correct_splice_dist(splice_df, splice_regions): 
	splice_df = splice_df.reset_index()
	variant_pos = splice_df.variant_id.map(data.rsid.pos)
	splice_strand = splice_df.splice_id.map(splice_regions.strand.replace({"+":1, "-":-1}))
	splice_distances = pd.concat({
		"snp_to_five": (variant_pos - splice_df.splice_id.map(splice_regions.tss_pos)) * splice_strand, 
		"snp_to_three": (variant_pos - splice_df.splice_id.map(splice_regions.tes_pos)) * splice_strand, 
	}, axis=1)
	splice_distances["closest"] = splice_distances["snp_to_five"].copy()
	splice_distances.loc[splice_distances.snp_to_three.abs() < splice_distances.snp_to_five.abs(), "closest"] = splice_distances.snp_to_three
	return splice_distances["closest"]



def _map_splice_id_to_genes(splice_regions): 
	exons = hg38.gencode_gtf[hg38.gencode_gtf.Feature == "exon"]
	introns = pr.PyRanges(pd.concat([
		hg38.gencode_gtf.features.introns(by="transcript").as_df(), 
		hg38.gencode_gtf.features.introns(by="gene").as_df()
	]))

	splice_regions = splice_regions.assign(cluster_id=_splice_id_to_cluster_id(splice_regions.index))

	# Splice IDs that exactly match an annotated exon
	_splice_intron_overlap = splice_regions.pr.join(introns[["gene_id"]]).as_df()
	_perfect_overlap_intron = _splice_intron_overlap[(_splice_intron_overlap.Start == _splice_intron_overlap.Start_b) & (_splice_intron_overlap.End == _splice_intron_overlap.End_b + 1)]

	# Splice IDs whose boundaries align with exons from the same gene
	correct_dtype_mismatch = lambda df: df.assign(Start=pd.to_numeric(df.Start, downcast="integer"), End=pd.to_numeric(df.End, downcast="integer"))
	# correct_dtype_mismatch = lambda df: df.assign(Start=df.Start.astype(np.int64), End=df.End.astype(np.int64))
	_splice_five_end_exon_overlap = splice_regions.tss.pr.join(exons.three_end().apply(correct_dtype_mismatch)).as_df()
	_splice_three_end_exon_overlap = splice_regions.tes.pr.join(exons.five_end().apply(correct_dtype_mismatch)).as_df()
	_perfect_overlap_exon = _splice_three_end_exon_overlap[["splice_id", "gene_id", "cluster_id"]].merge(_splice_five_end_exon_overlap[["splice_id", "gene_id"]], on=["splice_id", "gene_id"])

	perfect_overlap = pd.concat([_perfect_overlap_intron[["splice_id", "cluster_id", "gene_id"]], _perfect_overlap_exon]).drop_duplicates()

	unmapped_splice = splice_regions[~splice_regions.index.isin(perfect_overlap.splice_id)]

	splice_gene = pd.concat([
		perfect_overlap.assign(match="exact"), 
		_splice_five_end_exon_overlap[_splice_five_end_exon_overlap.splice_id.isin(unmapped_splice.index)][["splice_id", "cluster_id", "gene_id"]].drop_duplicates().assign(match="5_end"), 
		_splice_three_end_exon_overlap[_splice_three_end_exon_overlap.splice_id.isin(unmapped_splice.index)][["splice_id", "cluster_id", "gene_id"]].drop_duplicates().assign(match="3_end")
	])

	return splice_gene




def align_splice_id_to_gencode(splice_regions): 	

	introns = hg38.introns

	# shift end position to align with Gencode 
	splice_regions = splice_regions.assign(
		cluster_id = lambda df: df.index.str.split('_').str[1].astype(int), 
		end = lambda df: df['end'] - 1
	).df.reset_index()

	# Perform a series of matches between the junction sites and known boundaries
	exact_match = splice_regions.merge(introns, left_on=['chrom', 'start', 'end', 'strand'], right_on=['Chromosome', 'Start', 'End', 'Strand']).drop_duplicates(['splice_id', 'transcript_id'])
	# Sometimes the strand is flipped but coordinates match
	complement_match = splice_regions.merge(introns, left_on=['chrom', 'start', 'end'], right_on=['Chromosome', 'Start', 'End']).query('strand != Strand')
	# Find instances where only one junction site matches
	start_only_match = splice_regions.merge(introns, left_on=['chrom', 'start', 'strand'], right_on=['Chromosome', 'Start', 'Strand']).query('end != End')
	end_only_match = splice_regions.merge(introns, left_on=['chrom', 'end', 'strand'], right_on=['Chromosome', 'End', 'Strand']).query('start != Start')
	# In some cases, both ends of a excision event map to exon boundaries from the same gene but have not been annotated yet. These may be exon skipping events
	both_match = start_only_match.merge(end_only_match[['splice_id', 'gene_id']].drop_duplicates(), on=['splice_id', 'gene_id']).pipe(lambda df: df[~df['splice_id'].isin(exact_match.splice_id)]).drop_duplicates(['splice_id', 'gene_id'])

	# Merge all matches starting with most confident assignments
	matches = pd.concat({
		'exact': exact_match, 
		'exon_skip': both_match.pipe(lambda df: df[~df.splice_id.isin(exact_match.splice_id)]), 
	}, names=['match']).reset_index(level=0)
	matches = pd.concat([
		matches, 
		complement_match.pipe(lambda df: df[~df.splice_id.isin(matches.splice_id)]).assign(match='complement')
	])
	matches = pd.concat([
		matches, 
		start_only_match.pipe(lambda df: df[~df.splice_id.isin(matches.splice_id)]).assign(match='start_only'), 
		end_only_match.pipe(lambda df: df[~df.splice_id.isin(matches.splice_id)]).assign(match='end_only'), 
	])
	matches = matches.drop(columns=['Score', 'Frame', 'transcript_support_level', 'exon_number', 'exon_id', 'ont', 'havana_gene', 'ccdsid', 'protein_id', 'hgnc_id', 'level', 'havana_transcript'], errors='ignore')

	return matches

def map_splice_id_to_genes(splice_regions, matches, filt_genes, filt_splices): 
	"""
	We will sequentially subset the matches starting with most confident ones to assign genes to clusters to splice_ids. 
	"""
	# This checks for splice_ids with the same cluster_id that all point to the same gene and returns a series of cluster_id -> gene_id
	extract_unambiguous_clusters = lambda df: df.drop_duplicates(['cluster_id', 'gene_id']).pipe(lambda df: df[df.cluster_id.map(df.cluster_id.value_counts()) == 1]).set_index('cluster_id')['gene_id']

	# Unambiguous mappings using the strongest pieces of evidence
	unambiguous_clusters = matches.pipe(lambda df: df[df['match'].isin(['exact', 'complement', 'exon_skip'])]).pipe(extract_unambiguous_clusters)

	# Next get the remaining "ambiguous" clusters, and prioritize those that were measured in RNA or passed prior filtering. Of these, check for unambiguous pairs
	unambiguous_clusters = pd.concat([unambiguous_clusters, (
	    matches.pipe(lambda df: df[~df.cluster_id.isin(unambiguous_clusters.index)])
	    .pipe(lambda df: df[df.gene_id.isin(filt_genes)])
	    .pipe(lambda df: df[df.splice_id.isin(filt_splices)])
	    .pipe(lambda df: df[df['match'].isin(['exact', 'complement', 'exon_skip'])])
	    .pipe(extract_unambiguous_clusters)
	)])

	# Finally, if there are any clusters whose junctions align with the end of an exon, add those as well
	unambiguous_clusters = pd.concat([unambiguous_clusters, (
	    matches.pipe(lambda df: df[~df.cluster_id.isin(unambiguous_clusters.index)])
	    .pipe(lambda df: df[df.gene_id.isin(filt_genes)])
	    .pipe(lambda df: df[df.splice_id.isin(filt_splices)])
	    .pipe(extract_unambiguous_clusters)
	)])

	# We are left with matches with no clear concordance across the cluster. Genes will be assigned to splice_ids individually
	unmapped_matches = (
	    matches.pipe(lambda df: df[~df.cluster_id.isin(unambiguous_clusters.index)])
	    .pipe(lambda df: df[df.gene_id.isin(filt_genes)])
	    .pipe(lambda df: df[df.splice_id.isin(filt_splices)])
	)

	cluster_ids = splice_regions.assign(cluster_id = lambda df: df.index.str.split('_').str[1].astype(int))['cluster_id']
	splice_to_gene = pd.concat([
		cluster_ids.map(unambiguous_clusters).dropna().rename('gene_id').to_frame().reset_index(), 
		unmapped_matches[['splice_id', 'gene_id']]
	])
	return splice_to_gene




fops = FeatureOps()