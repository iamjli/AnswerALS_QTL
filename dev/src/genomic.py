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

	def __init__(self, omic, bam_paths, residualizer=None, initialize=True, prefiltered=False, max_cache=12): 

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

	def get_coverages_pathos(self, chrom, start, end): 
		coverages = {}
		for i,(guid,bam) in enumerate(self.sam_files_dic.items()): 
			coverages[guid] = self._get_coverage(bam, chrom, start, end, prefiltered=True)["complete_overlap"]
		results = { guid:val["pos"]+val["neg"] for guid,val in coverages.items() }
		return pd.Series(results)

	@parse_interval_args
	def get_coverages(self, chrom, start, end, strand=None, report_strand=None, overlap="any", norm_factors=None, verbose=True):

		# load data if cached, otherwise continue
		region = "{}:{}-{}".format(chrom, start, end)
		if region in self._coverage_cache: 
			logger.write("Loading {} coverages from cache".format(self.omic), verbose=verbose)
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
			logger.write("overlap argument must be `any` or `complete`", verbose=verbose)

		# get results relative to strand
		if strand == "-": 
			pos_results = pd.Series({ guid:val["neg"] for guid,val in results.items() })
			neg_results = pd.Series({ guid:val["pos"] for guid,val in results.items() })
		else: 
			pos_results = pd.Series({ guid:val["pos"] for guid,val in results.items() })
			neg_results = pd.Series({ guid:val["neg"] for guid,val in results.items() })

		if isinstance(norm_factors, pd.Series): 
			logger.write("Using custom normalization", verbose=verbose)
			pass
		elif norm_factors is None: 
			logger.write("Using normalization set during initialization", verbose=verbose)
			norm_factors = self.norm_factors
		elif norm_factors == 1: 
			norm_factors = 1
		else: 
			assert norm_factors == "lib_size"
			logger.write("Using library normalization (consider using edgeR normalization)", verbose=verbose)
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
			logger.write("Loading {} pileups from cache".format(self.omic), verbose=verbose)
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
			logger.write("Using custom normalization", verbose=verbose)
			pass
		elif norm_factors is None: 
			logger.write("Using normalization set during initialization", verbose=verbose)
			norm_factors = self.norm_factors
		elif norm_factors == 1: 
			norm_factors = 1
		else: 
			assert norm_factors == "lib_size"
			logger.write("Using library normalization (consider using edgeR normalization)", verbose=verbose)
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

		# Tracks should be library-normalized but not residualized yet
		self.pos_tracks = Coverage(pos_tracks) 
		self.neg_tracks = Coverage(neg_tracks)
		self.tracks = self.pos_tracks + self.neg_tracks

		self.residualizer = residualizer
		self.covariates = self.residualizer.C.T

		# self.positions = self.tracks.index

	@classmethod
	def load(cls, genomic_data, residualizer, interval_obj, **kwargs): 
		chrom, start, end = interval_obj.coords
		pos_tracks, neg_tracks = genomic_data.get_pileups(chrom, start, end, strand="+", report_strand="both", **kwargs)
		return cls(pos_tracks, neg_tracks, residualizer)

	def scale(self, scale): 
		return TrackPlotter(self.pos_tracks * scale, self.neg_tracks * scale, self.residualizer)

	# @staticmethod
	# def center_positions()

	@staticmethod
	def coursen(track_data, bounds=None, center=False, bin_size=None, n_bins=999, **kwargs): 
		"""
		Subsamples the track evenly for more manageable plotting. Can handle series or dataframe inputs.
		"""
		if bounds is None: 
			bounds = (track_data.index.min(), track_data.index.max())

		# subset the track
		track_data = track_data.loc[bounds[0]: bounds[1]+1]

		if bin_size is not None: 
			n_bins = len(track_data) // bin_size

		if len(track_data) > n_bins: 
			# series of original positions to bin positions
			bins_s = pd.Series(index=track_data.index, data=pd.cut(track_data.index, bins=n_bins, retbins=False))
			# bins = np.arange(track_data.index[0], track_data.index[-1], step=)
			# bins_s = pd.Series(index=track_data.index, )
			# set bin labels to midpoints
			bins_s = bins_s.apply(lambda pos: pos.mid.round()).astype(int)
			track_data = track_data.groupby(bins_s).mean()

		if center: 
			midpoint = (bounds[0] + bounds[1]) // 2
			track_data.index = track_data.index - midpoint

		return track_data

	def get_tracks(self, strand=None, residualize=True, coursen=True, **kwargs):
		if strand == "pos": 
			tracks = self.pos_tracks
		elif strand == "neg": 
			tracks = self.neg_tracks
		elif strand is None: 
			tracks = self.tracks

		if coursen: 
			tracks = TrackPlotter.coursen(tracks, **kwargs)

		if residualize: 
			tracks = self.residualizer.transform(tracks)

		return tracks

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

	@staticmethod
	def _plotter(track_data_series, bounds=None, coursen=True, scale=1, ax=None, invert=False, orient="h", **kwargs):
		"""
		Plots a single track with given bounds. Subsamples the data for faster plotting. 

		`track_data_series` - series indexed by position
		"""
		if invert: track_data_series *= -1
		track_data_series *= scale

		if coursen: 
			track_data_series = TrackPlotter.coursen(track_data_series, bounds=bounds, **kwargs)
		
		plot_bounds = (track_data_series.index.min(), track_data_series.index.max())

		if orient=="h": 
			# plot track horizontally
			sns.lineplot(data=track_data_series, ax=ax, legend=False)
			ax.set(xlim=plot_bounds)
		else: 
			# plot track vertically
			x,y = TrackPlotter.rotate(track_data_series.index, track_data_series.values)
			ax.plot(x,y)
			ax.set(ylim=plot_bounds)

	def plot_track_mean(self, strand=None, residualize=True, plot_feature="mean", **kwargs): 

		if strand == "both": 
			self.plot_track_mean(strand="pos", residualize=True, plot_feature=plot_feature, **kwargs)
			self.plot_track_mean(strand="neg", residualize=True, plot_feature=plot_feature, invert=True, **kwargs)
		else: 
			tracks = self.get_tracks(strand=strand, residualize=True, **kwargs)
			if plot_feature == "mean": 
				track = tracks.mean(axis=1)
			elif plot_feature == "variance_norm": 
				track = tracks.std(axis=1) / tracks.mean(axis=1)
			else: 
				track = tracks.iloc[:,plot_feature]
			TrackPlotter._plotter(track, coursen=False, **kwargs)

	def plot_tracks_by_genotype(self, genotypes, strand=None, **kwargs): 

		if strand == "both":
			self.plot_tracks_by_genotype(genotypes=genotypes, strand="pos", **kwargs)
			self.plot_tracks_by_genotype(genotypes=genotypes, strand="neg", invert=True, **kwargs)
		else: 
			tracks = self.get_tracks(strand=strand, residualize=True, **kwargs)
			# group tracks by genotypes
			grouped_tracks = tracks.groupby(genotypes, axis=1).mean()
			for col in grouped_tracks.columns.sort_values(): 
				self._plotter(grouped_tracks[col], coursen=False, **kwargs)

	def plot_explained_variance(self, genotypes, strand=None, **kwargs):

		if strand == "both": 
			self.plot_explained_variance(genotypes=genotypes, strand="pos", **kwargs)
			self.plot_explained_variance(genotypes=genotypes, strand="neg", **kwargs)
		else: 
			tracks = self.get_tracks(strand=strand, residualize=False, **kwargs)
			results = QTL_pairwise(genotypes, tracks, self.covariates)
			results = results.sort_values("phenotype_id").set_index("phenotype_id")["r2"].fillna(0)
			self._plotter(results, coursen=False, **kwargs)


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


	def draw_variant_line(self, variant_id=None, ax=None, **kwargs): 

		variants = np.array(variant_id).flatten()
		positions = DATA.rsid.loc[variants, "pos"]
		if kwargs.get("center", False):
			bounds = kwargs["bounds"]
			mid = (bounds[0] + bounds[1]) // 2
			positions -= mid
		for p in positions: 
			ax.axvline(x=p)

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

	@staticmethod
	def plot_correlation(x, y, x_plot_params, y_plot_params, r_max=1, ax=None): 

		# residualize tracks here, then don't provide residual later
		x_tracks = x.get_tracks(strand=None, residualize=True, coursen=True, **x_plot_params)
		y_tracks = y.get_tracks(strand=None, residualize=True, coursen=True, **y_plot_params)

		r_nominal_df = QTL_pairwise(x_tracks, y_tracks, covariates_df=None, return_r_matrix=True)
		r_nominal_df = r_nominal_df.T[::-1]

		sns.heatmap(r_nominal_df, 
					center=0, vmin=-r_max, vmax=r_max, cmap="bwr", cbar=False, 
					xticklabels=False, yticklabels=False, ax=ax)

		return r_nominal_df

	@staticmethod
	def figure_track_correlation(atac_tracks, rna_tracks, x_plot_params, y_plot_params, r_max=1, scale_rna=1, figsize=(15,15), x_ratios=[1,6,1], y_ratios=[1,6,1]): 

		with plt.rc_context(dict(**sns.plotting_context("notebook", font_scale=1.5))):
			
			fig, axes = plt.subplots(3, 3, figsize=figsize, sharex=False, sharey=False, 
														 gridspec_kw={'width_ratios': x_ratios, 'height_ratios': y_ratios, 
																	 'wspace':0, 'hspace':0})

			axes[0][0].axis("off")
			axes[0][2].axis("off")
			axes[2][0].axis("off")
			axes[2][2].axis("off")

			x = rna_tracks if x_plot_params.get("omic", "rna") == "rna" else atac_tracks
			y = atac_tracks if y_plot_params.get("omic", "atac") == "atac" else rna_tracks

			x_plot_params["center"] = x_plot_params.get("center", True)
			y_plot_params["center"] = y_plot_params.get("center", True)

			x_plot_params["n_bins"] = x_plot_params.get("n_bins", 500)
			y_plot_params["n_bins"] = y_plot_params.get("n_bins", 250)

			# plot pos track on top
			y.plot_track_mean(**x_plot_params, strand="pos", ax=axes[0][1])
			x.plot_track_mean(**x_plot_params, strand="pos", ax=axes[0][1])
			# neg track on bottom
			y.plot_track_mean(**x_plot_params, strand="neg", ax=axes[2][1], invert=True)
			x.plot_track_mean(**x_plot_params, strand="neg", ax=axes[2][1], invert=True)

			y.plot_track_mean(**y_plot_params, strand="pos", orient="v", ax=axes[1][0])
			x.plot_track_mean(**y_plot_params, strand="pos", orient="v", ax=axes[1][0])

			y.plot_track_mean(**y_plot_params, strand="neg", orient="v", ax=axes[1][2], invert=True)
			x.plot_track_mean(**y_plot_params, strand="neg", orient="v", ax=axes[1][2], invert=True)
			axes[1][2].set_yticklabels([])

			TrackPlotter.plot_correlation(x, y, x_plot_params, y_plot_params, r_max=r_max, ax=axes[1][1])

			# print(axes[1][1].get_ylim())
			# print(axes[1][0].get_ylim())


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


