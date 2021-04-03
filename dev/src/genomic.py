#!/usr/bin/env python3
import re
from collections import OrderedDict

import pandas as pd
import numpy as np
import pyranges as pr
import pysam

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from . import DATA, BASE_DIR, RESULTS_PATHS, CHROMS, logger
from .bed import Interval
from .qtl import QTL_pairwise, Residualizer
from .utils import flatten


gt_to_dosage_dict = {'0/0':0, '0/1':1, '1/1':2, './.':np.nan,
					 '0|0':0, '0|1':1, '1|0':1, '1|1':2, '.|.':np.nan}


def parse_interval_args(func): 
	"""Decorator to allow interval arguments as opposed to specifying coordinates."""
	def _interval_parser(self, *args, **kwargs): 
		if type(args[0]).__name__ == "Interval": 
			return func(self, **{**args[0].as_dict, **kwargs})
		else:
			return func(self, *args, **kwargs)
	return _interval_parser


class BamQuery: 

	def __init__(self, omic, bam_paths, residualizer=None, initialize=True, prefiltered=False, max_cache=3): 

		self.omic = omic
		self.bam_paths = bam_paths
		self.residualizer = residualizer
		self.prefiltered = prefiltered

		self.sample_names = self.bam_paths.index

		self._sam_files_dic = None
		self._library_sizes = None

		self.norm_factors = 1

		self._max_cache = max_cache
		self._pileup_cache = OrderedDict()
		self._coverage_cache = OrderedDict()

		if initialize: 
			self.sam_files_dic

	@classmethod
	def load_rna(cls, bam_paths=DATA.bam_paths, **kwargs): 
		logger.write("Loading RNA bams specified in: {}".format(RESULTS_PATHS["bams"].relative_to(BASE_DIR)))
		residualizer = Residualizer.load_rna()
		return cls(omic="rna", bam_paths=bam_paths["rna_bam"], residualizer=residualizer, **kwargs)

	@classmethod
	def load_atac(cls, bam_paths=DATA.bam_paths, **kwargs): 
		logger.write("Loading ATAC bams specified in: {}".format(RESULTS_PATHS["bams"].relative_to(BASE_DIR)))
		residualizer = Residualizer.load_rna()
		return cls(omic="atac", bam_paths=bam_paths["atac_bam"], residualizer=residualizer, **kwargs)

	def _update_coverage_cache(self, region, coverages):
		self._coverage_cache[region] = coverages

	def _update_pileup_cache(self, region, pileups): 
		"""Updates cache with pileups formated as a dict keyed by `pos` and `neg`."""
		if len(self._pileup_cache) >= self._max_cache: 
			del self._pileup_cache[next(iter(self._pileup_cache))]
		self._pileup_cache[region] = pileups 

	def _reset_caches(self): 
		self._pileup_cache = OrderedDict()
		self._coverage_cache = OrderedDict()

	def set_norm_factors(self, norm_factors=1): 
		self.norm_factors = norm_factors

	@property
	def sam_files_dic(self): 
		"""Dictionary of pysam objects."""
		if self._sam_files_dic is None: 
			_sam_files_dic = {}
			for i,(guid,path) in enumerate(self.bam_paths.items()): 
				logger.update("Loading {} bams as pysam objects ({} of {})...".format(self.omic.upper(), i, len(self.bam_paths)))
				_sam_files_dic[guid] = pysam.AlignmentFile(path, "rb", threads=4)
			logger.flush()
			self._sam_files_dic = _sam_files_dic
			return self._sam_files_dic
		else: 
			return self._sam_files_dic

	@property
	def library_sizes(self):
		"""Total library size, including reads that don't pass QC."""
		if self._library_sizes is None: 
			sizes = {}
			for i,(guid,path) in enumerate(self.bam_paths.iteritems()): 
				logger.update("Scanning {} files for library_sizes ({} of {})...".format(self.omic.upper(), i, len(self.bam_paths)))
				results = pysam.idxstats(path)
				sizes[guid] = np.array([row.split("\t")[2:] for row in results.rstrip("\n").split("\n")], dtype=int).sum()
			logger.flush()
			self._library_sizes = pd.Series(sizes)
		return self._library_sizes

	## Genomic coverage operations
	@staticmethod
	def _match_chrom_prefix(pysam_obj, chrom): 
		"""Returns the correctly formatted chromosome ID corresponding to a bam file (i.e. `chr1` vs. `1`)."""
		if "chr" in pysam_obj.references[0]: 
			if "chr" not in chrom: 
				chrom = "chr" + str(chrom)
		else: 
			if "chr" in chrom: 
				chrom = chrom[3:]
		return chrom

	@staticmethod
	def _get_coverage(pysam_obj, chrom, start, end, prefiltered=False): 

		chrom = BamQuery._match_chrom_prefix(pysam_obj, chrom)

		complete_overlap = {"pos": 0, "neg": 0}
		partial_overlap = {"pos": 0, "neg": 0}

		for read in pysam_obj.fetch(chrom, start, end):
			if not prefiltered: 
				if read.is_duplicate | read.is_qcfail | read.is_secondary | read.is_unmapped | (not read.is_proper_pair): 
					continue # move on to next read

			if (read.reference_start >= start) & (read.reference_end <= end): 
				if read.is_reverse: 
					complete_overlap["neg"] += 1
				else: 
					complete_overlap["pos"] += 1
			else: 
				if read.is_reverse: 
					partial_overlap["neg"] += 1
				else: 
					partial_overlap["pos"] += 1

		any_overlap = {
			"pos": complete_overlap["pos"]+partial_overlap["pos"], 
			"neg": complete_overlap["neg"]+partial_overlap["neg"], 
		}

		# return counts relative to strand, if specified
		return {"complete_overlap": complete_overlap, "any_overlap": any_overlap}

	@staticmethod
	def _get_pileup(pysam_obj, chrom, start, end, include_introns=False, prefiltered=False, **kwargs): 
		"""
		Gets base-by-base read counts within a region from a pysam bam object.
		"""
		chrom = BamQuery._match_chrom_prefix(pysam_obj, chrom)

		region_len = end - start
		pos_markers = np.zeros(region_len+1, dtype=int)
		neg_markers = np.zeros(region_len+1, dtype=int)
		filtered_counts = 0

		# iterate through reads that overlap query region
		for read in pysam_obj.fetch(chrom, start, end): 
			if not prefiltered: 
				if read.is_duplicate | read.is_qcfail | read.is_secondary | read.is_unmapped | (not read.is_proper_pair): 
					# filter out reads that don't pass QC. This has already been done for ATAC bams, but not for RNA. 
					filtered_counts += 1
					continue # move on to next read

			# Some reads are spliced and are composed of non-contiguous regions. This specifies whether we should includ those regions.
			if include_introns:
				blocks = [(read.reference_start, read.reference_end)] 
			else: 
				blocks = read.get_blocks()

			for b_start,b_end in blocks: 
				if (b_start <= end) & (b_end >= start):
					if b_start < start: 
						rel_start = 0
					else: 
						rel_start = b_start - start
					if b_end > end: 
						rel_end = region_len
					else: 
						rel_end = b_end - start

					if read.is_reverse: 
						neg_markers[rel_start] += 1
						neg_markers[rel_end] -= 1
					else: 
						pos_markers[rel_start] += 1
						pos_markers[rel_end] -= 1

		pos_counts = pd.Series(data=pos_markers.cumsum()[:-1], index=range(start,end))
		neg_counts = pd.Series(data=neg_markers.cumsum()[:-1], index=range(start,end))

		return {"pos": pos_counts, "neg": neg_counts}

	@parse_interval_args
	def get_coverages(self, chrom, start, end, strand=None, report_strand=None, overlap="any", norm_factors=None, verbose=True):

		# load data if cached, otherwise continue
		region = "{}:{}-{}".format(chrom, start, end)
		if region in self._coverage_cache: 
			logger.write("Loading {} coverages from cache".format(self.omic))
			coverages = self._coverage_cache[region]
		else: 
			coverages = {}
			for i,(guid,bam) in enumerate(self.sam_files_dic.items()): 
				logger.update("Loading {} bam coverages ({} of {})...".format(self.omic, i, len(self.sam_files_dic)), verbose=verbose)
				coverages[guid] = self._get_coverage(bam, chrom, start, end)
			logger.flush(verbose=verbose)
			self._update_coverage_cache(region, coverages)

		# report partial or complete overlaps
		if overlap == "any": 
			results = { guid:val["any_overlap"] for guid,val in coverages.items() }
		elif overlap == "complete": 
			results = { guid:val["complete_overlap"] for guid,val in coverages.items() }
		else: 
			logger.write("overlap argument must be `any` or `complete`")

		# get results relative to strand
		if strand == "-": 
			pos_results = pd.Series({ guid:val["neg"] for guid,val in results.items() })
			neg_results = pd.Series({ guid:val["pos"] for guid,val in results.items() })
		else: 
			pos_results = pd.Series({ guid:val["pos"] for guid,val in results.items() })
			neg_results = pd.Series({ guid:val["neg"] for guid,val in results.items() })

		if isinstance(norm_factors, pd.Series): 
			logger.write("Using custom normalization")
			pass
		elif norm_factors is None: 
			logger.write("Using normalization set during initialization")
			norm_factors = self.norm_factors
		elif norm_factors == 1: 
			norm_factors = 1
		else: 
			assert norm_factors == "lib_size"
			logger.write("Using library normalization (consider using edgeR normalization)")
			norm_factors = self.library_sizes / 1e6

		pos_results /= norm_factors
		neg_results /= norm_factors

		if report_strand is None: 
			return pos_results + neg_results
		elif report_strand == "forward": 
			return pos_results
		elif report_strand == "reverse": 
			return neg_results
		else: 
			return pos_results, neg_results

	@parse_interval_args
	def get_pileups(self, chrom, start, end, strand=None, report_strand=None, include_introns=False, norm_factors=None, verbose=True): 
		# load data if cached, otherwise continue
		region = "{}:{}-{}".format(chrom, start, end)
		if include_introns: region += "_introns"
		if region in self._pileup_cache: 
			logger.write("Loading {} pileups from cache".format(self.omic))
			pileups = self._pileup_cache[region]
		else: 
			pileups = {}
			for i,(guid,bam) in enumerate(self.sam_files_dic.items()): 
				logger.update("Loading {} bam pileups ({} of {})...".format(self.omic, i, len(self.sam_files_dic)), verbose=verbose)
				pileups[guid] = self._get_pileup(bam, chrom, start, end, include_introns=include_introns)
			logger.flush(verbose=verbose)
			self._update_pileup_cache(region, pileups)

		# get results relative to strand
		if strand == "-": 
			pos_results = pd.DataFrame.from_dict({guid:x["neg"] for guid,x in pileups.items()})
			neg_results = pd.DataFrame.from_dict({guid:x["pos"] for guid,x in pileups.items()})
		else: 
			pos_results = pd.DataFrame.from_dict({guid:x["pos"] for guid,x in pileups.items()})
			neg_results = pd.DataFrame.from_dict({guid:x["neg"] for guid,x in pileups.items()})

		if isinstance(norm_factors, pd.Series): 
			logger.write("Using custom normalization")
			pass
		elif norm_factors is None: 
			logger.write("Using normalization set during initialization")
			norm_factors = self.norm_factors
		elif norm_factors == 1: 
			norm_factors = 1
		else: 
			assert norm_factors == "lib_size"
			logger.write("Using library normalization (consider using edgeR normalization)")
			norm_factors = self.library_sizes / 1e6

		pos_results /= norm_factors
		neg_results /= norm_factors

		if report_strand is None: 
			return Coverage(pos_results + neg_results)
		elif report_strand == "forward": 
			return Coverage(pos_results)
		elif report_strand == "reverse": 
			return Coverage(neg_results)
		else: 
			return Coverage(pos_results), Coverage(neg_results)

	@parse_interval_args
	def plotter(self, chrom, start, end, norm_factors=None, residualizer="self", **kwargs): 
		"""
		Returns object for plotting pileups.
		"""
		kwargs["strand"] = "+"
		pos_pileups, neg_pileups = self.get_pileups(chrom, start, end, report_strand="both", norm_factors=norm_factors, **kwargs)

		if residualizer == "self": 
			residualizer = self.residualizer

		return TrackPlotter(pos_pileups, neg_pileups, residualizer)


class TrackPlotter: 

	def __init__(self, pos_tracks, neg_tracks, residualizer=None): 

		self.pos_tracks = Coverage(pos_tracks)
		self.neg_tracks = Coverage(neg_tracks)
		self.tracks = self.pos_tracks + self.neg_tracks
		self.residualizer = residualizer

		self.positions = self.tracks.index

	@classmethod
	def load(cls, genomic_data, residualizer, interval_obj, **kwargs): 
		chrom, start, end = interval_obj.coords
		pos_tracks, neg_tracks = genomic_data.get_pileups(chrom, start, end, strand="+", report_strand="both", **kwargs)
		return cls(pos_tracks, neg_tracks, residualizer)

	def estimate_coverage(self, bounds, by_strand=False, residualize=True, n_bins=1): 

		if by_strand: 
			pos_coverage = TrackPlotter.coursen(self.pos_tracks, bounds, n_bins=n_bins)
			neg_coverage = TrackPlotter.coursen(self.neg_tracks, bounds, n_bins=n_bins)
			if residualize: 
				pos_coverage = self.residualizer.transform(pos_coverage)
				neg_coverage = self.residualizer.transform(neg_coverage)
			return pos_coverage, neg_coverage

		else: 
			# get coverage by non-overlapping intervals. When n_bins is 1, it's just the mean over the entire interval
			coverage = TrackPlotter.coursen(self.tracks, bounds, n_bins=n_bins)
			if residualize: 
				coverage = self.residualizer.transform(coverage)
			return coverage

	@staticmethod
	def get_bounds(data, bounds=None, w=0, shift=0): 
		if bounds is None: 
			center_pos = data.index[len(data) // 2]
			l_bound = center_pos
			r_bound = center_pos + 1
		else: 
			l_bound, r_bound = bounds
		lower = l_bound - w + shift
		upper = r_bound + w + shift
		return lower, upper

	@staticmethod
	def coursen(data, bounds=None, bin_size=None, n_bins=1000): 
		if bounds is not None: 
			data = data.loc[bounds[0]: bounds[1]]
		if bin_size is None: 
			if len(data) < n_bins: 
				return data
			bin_size = len(data) / n_bins
		# create bins of size `bin_size` starting at the first entry
		bins = (data.index - data.index.min()) // bin_size
		# scale back to original positions
		bins = (bins * bin_size).astype(int) + data.index.min()
		return data.groupby(bins).mean()

	@staticmethod
	def subselect(data, bounds=None, w=0, shift=0, n_bins=1000, **kwargs): 

		if "interval" in kwargs: 
			bounds = kwargs["interval"].bounds

		# get plotting window
		bounds = TrackPlotter.get_bounds(data, bounds, w, shift)

		# get subsetted data
		data = TrackPlotter.coursen(data, bounds=bounds, n_bins=n_bins)

		return data

	@staticmethod
	def rotate(x, y):
		origin = x[0], 0

		# shift to origin
		x1 = x - origin[0]
		y1 = y - origin[1]

		#rotate
		x2 = -y1
		y2 = x1

		# shift back
		x3 = x2 + origin[1]
		y3 = y2 + origin[0]

		return x3, y3 


	def _plotter(self, track_data, residualizer=None, genotypes=None, bounds=None, w=None, shift=0, subselect=True, ax=None, orient="h", **kwargs): 
		"""
		track_data - dataframe indexed by position and with sample headers
		"""
		# get plotting window
		# bounds = TrackPlotter.get_bounds(track_data, bounds, w, shift)

		# # get subsetted data
		# track_data = TrackPlotter.coursen(track_data, bounds=bounds, n_bins=1000)
		if "interval" in kwargs: 
			bounds = kwargs["interval"].bounds

		# reduce data dimension
		if subselect: 
			track_data = TrackPlotter.subselect(track_data, bounds, w, shift)

		bounds = (track_data.index.min(), track_data.index.max())

		# residualize
		if residualizer is not None: 
			track_data = residualizer.transform(track_data)

		# collapse dataframe, otherwise assume data is a series
		if isinstance(track_data, pd.DataFrame): 
			# get average (across genotypes if provided)
			if genotypes is not None: 
				track_data = track_data.groupby(genotypes, axis=1).mean()
			else: 
				track_data = track_data.mean(axis=1)

		if orient=="h": 
			sns.lineplot(data=track_data, ax=ax, legend=False)
			ax.set(xlim=bounds)
		else: 
			x,y = TrackPlotter.rotate(track_data.index, track_data.values)
			ax.plot(x,y)
			ax.set(ylim=bounds)

	def plot_track(self, strand=None, by_strand=False, **kwargs): 
		if strand == "pos": 
			self._plotter(self.pos_tracks, residualizer=self.residualizer, **kwargs)
		elif strand == "neg": 
			self._plotter(-1 * self.neg_tracks, residualizer=self.residualizer, **kwargs)
		elif strand == "both": 
			self._plotter(self.pos_tracks, residualizer=self.residualizer, **kwargs)
			self._plotter(-1 * self.neg_tracks, residualizer=self.residualizer, **kwargs)
		else: 
			self._plotter(self.tracks, residualizer=self.residualizer, **kwargs)

	# def plot_by_strand(self, **kwargs): 

	# 	self._plotter(self.pos_tracks, residualizer=self.residualizer, **kwargs)
	# 	self._plotter(-1 * self.neg_tracks, residualizer=self.residualizer, **kwargs)

	def plot_normalized_variance(self, **kwargs): 
		pass

	def plot_explained_variance(self, genotypes=None, by_strand=False, **kwargs):

		if by_strand: 
			pos_tracks = TrackPlotter.subselect(self.pos_tracks, **kwargs)
			pos_results = QTL_pairwise(genotypes, pos_tracks, self.residualizer.C.T)
			pos_results = pos_results.sort_values("phenotype_id").set_index("phenotype_id")["r2"].fillna(0)
			self._plotter(pos_results, **kwargs)

			neg_tracks = TrackPlotter.subselect(self.neg_tracks, **kwargs)
			neg_results = QTL_pairwise(genotypes, neg_tracks, self.residualizer.C.T)
			neg_results = neg_results.sort_values("phenotype_id").set_index("phenotype_id")["r2"].fillna(0)
			self._plotter(neg_results, **kwargs)
		else: 
			tracks = TrackPlotter.subselect(self.tracks, **kwargs)
			results = QTL_pairwise(genotypes, tracks, self.residualizer.C.T)
			results = results.sort_values("phenotype_id").set_index("phenotype_id")["r2"].fillna(0)
			self._plotter(results, **kwargs)
			return results

	def plot_correlation_with_feature(self, feature, residualize=True, **kwargs): 
		"""Feature-centric plot. How does the nearby area correlate with this feature?"""
		if residualize: 
			results = QTL_pairwise(feature, self.tracks, self.residualizer.C.T)
		else: 
			tracks = self.residualizer.transform(self.tracks)
			results = QTL_pairwise(feature, tracks)
		results = results.rename(columns={"variant_id": "feature_bin", "phenotype_id": "position"})
		results = results.sort_values("position").set_index("position")
		self._plotter(results["r2"].fillna(0) * np.sign(results["slope"]), **kwargs)

	def plot_correlation_with_feature_heatmap(self, feature_bins, n_bins=1000, r2_max=1, ax=None, **kwargs): 
		feature_bins.index = (feature_bins.index - np.median(feature_bins.index)).astype(int)
		tracks = TrackPlotter.subselect(self.tracks, n_bins=n_bins, **kwargs)
		results = QTL_pairwise(
			feature_bins, 
			tracks, 
			self.residualizer.C.T
		).rename(columns={"variant_id": "feature_bin", "phenotype_id": "position"})
		results["r2"] *= np.sign(results["slope"])
		results_2d = results.pivot(index="feature_bin", columns="position", values="r2")
		results_2d = results_2d.reindex(index=feature_bins.index, columns=tracks.index)

		sns.heatmap(data=results_2d, center=0, vmin=-r2_max, vmax=r2_max, cmap="bwr", cbar=False, ax=ax)

	def plot_qtl_pvals(self, qtl_results=None, ax=None, **kwargs): 
		# get unique pvals for each variant
		qtl_results = qtl_results.sort_values("pval_nominal").drop_duplicates("variant_id")
		# add position information
		qtl_results["pos"] = qtl_results["variant_id"].map(DATA.rsid["pos"])
		# select rows where variant is within plotting range
		bounded_df = qtl_results[qtl_results["pos"].between(*ax.get_xlim())]

		plot_s = -np.log10(qtl_results.set_index("pos")["pval_nominal"])
		sns.scatterplot(x=plot_s.index, y=plot_s.values, ax=ax)
		ax.set_ylim([0, plot_s.values.max() * 1.1])


	def draw_variant_line(self, variant_id=None, ax=None, **kwargs): 

		variants = np.array(variant_id).flatten()
		positions = DATA.rsid.loc[variants, "pos"]
		for p in positions: 
			ax.axvline(x=p)





from .omic import inverse_normal_transform

class Coverage(pd.DataFrame): 
	"""
	Methods for coverage dataframe. The difference between this and `omic` is that this is indexed by genome position. 
	"""
	@property
	def _constructor(self):
		return Coverage

	def normalize(self, norm_factors=1, inverse_norm_transform=True): 
		df = self / norm_factors
		if inverse_norm_transform: 
			df = inverse_normal_transform(df)
		return df

	def residualize(self, residualizer): 
		self_t = torch.tensor(self.values, dtype=torch.float).to("cpu")
		return pd.DataFrame(residualizer.transform(self_t), index=self.index, columns=self.columns)

	def bin(self, bin_size=None, n=None): 
		if n: 
			bin_size = int(len(self) / n)
		df = self.groupby(self.index // bin_size).mean()
		df.index = df.index * bin_size
		return df



class QueryVariants: 

	def __init__(self, vcf_path=None, rsid=None, initialize=False): 

		self.vcf_path = vcf_path
		self.rsid = rsid

		self.tbx = pysam.TabixFile(str(self.vcf_path))




class Genomic: 

	def __init__(self, vcf_path=None, metadata=None, bam_paths=None, rsid=None, R=DATA, initialize=False): 

		self.DATA = DATA

		if self.DATA is not None: 
			self.vcf_path  = RESULTS_PATHS["vcf"]
			self.metadata  = self.DATA.metadata
			self.bam_paths = self.DATA.bam_paths
			self.rsid	  = self.DATA.rsid
		else: 
			from .repository import _load_metadata, _load_bam_paths, _load_rsid
			self.metadata  = _load_metadata()
			self.bam_paths = _load_bam_paths()
			# self.rsid	  = _load_rsid()

		self.tbx = pysam.TabixFile(str(self.vcf_path))
		# self._chr_prefix = "chr" in self.tbx.references[0]

		# Ensure we're working with the same samples
		self.sample_names = pd.Index(self.tbx.header[-1].split("\t")[9:])
		assert self.sample_names.isin(self.metadata.index).all()
		assert self.sample_names.isin(self.bam_paths.index).all()
		self.metadata = self.metadata.reindex(self.sample_names)

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

		# self.bam_paths = self.bam_paths.reindex(self.sample_names)

		# self._rna_bams_pysam = None
		# self._atac_bams_pysam = None

		# self._rna_library_sizes = None
		# self._atac_library_sizes = None

		# if initialize: 
		# 	self.rna_bams_pysam
		# 	self.atac_bams_pysam





	# @parse_interval_args
	# def get_coverages_old(self, chrom, start, end, strand=None, **kwargs): 

	# 	strand_specific = kwargs.pop("strand_specific", True)
	# 	normalize = kwargs.pop("normalize", "lib_size")
	# 	verbose = kwargs.pop("verbose", True)
	# 	# logger.write("Unused kwargs:", kwargs)

	# 	coverages = {}
	# 	for i,(guid,bam) in enumerate(self.sam_files_dic.items()): 
	# 		logger.update("Loading bam coverage ({} of {})...".format(i, len(self.sam_files_dic)), verbose=verbose)
	# 		coverages[guid] = self._get_coverage(bam, chrom, start, end, strand)
	# 	logger.flush(verbose=verbose)

	# 	pos_results = pd.Series({ guid:val["pos"] for guid,val in coverages.items() })
	# 	neg_results = pd.Series({ guid:val["neg"] for guid,val in coverages.items() })

	# 	if isinstance(normalize, pd.Series): 
	# 		norm_factors = normalize
	# 	elif normalize == "lib_size": 
	# 		norm_factors = self.library_sizes / self.library_sizes.median()
	# 	else: 
	# 		norm_factors = 1
	# 	pos_results = pos_results / norm_factor
	# 	neg_results = neg_results / norm_factor

	# 	if strand_specific: 
	# 		return pos_results, neg_results
	# 	else: 
	# 		return pos_results + neg_results

	# @parse_interval_args
	# def get_pileups_old(self, chrom, start, end, strand=None, **kwargs): 

	# 	strand_specific = kwargs.pop("strand_specific", False)
	# 	include_introns = kwargs.pop("include_introns", False)
	# 	normalize = kwargs.pop("normalize", "lib_size")
	# 	verbose = kwargs.pop("verbose", True)
	# 	# logger.write("Unused kwargs:", kwargs)

	# 	pileups = {}
	# 	for i,(guid,bam) in enumerate(self.sam_files_dic.items()): 
	# 		logger.update("Loading bam pileups ({} of {})...".format(i, len(self.sam_files_dic)), verbose=verbose)
	# 		pileups[guid] = self._get_pileup(bam, chrom, start, end, strand, include_introns=include_introns)
	# 	logger.flush(verbose=verbose)

	# 	pos_results = pd.DataFrame.from_dict({guid:x["pos"] for guid,x in pileups.items()})
	# 	neg_results = pd.DataFrame.from_dict({guid:x["neg"] for guid,x in pileups.items()})

	# 	if isinstance(normalize, pd.Series): 
	# 		norm_factors = normalize
	# 	elif normalize == "lib_size": 
	# 		norm_factors = self.library_sizes / self.library_sizes.median()
	# 	else: 
	# 		norm_factors = 1
	# 	pos_results = pos_results / norm_factor
	# 	neg_results = neg_results / norm_factor

	# 	if strand_specific: 
	# 		return pos_results, neg_results
	# 	else: 
	# 		return pos_results + neg_results
			
	# 	if strand_specific: 
	# 		return Coverage(pos_results), Coverage(neg_results)
	# 	else: 
	# 		return Coverage(pos_results + neg_results)