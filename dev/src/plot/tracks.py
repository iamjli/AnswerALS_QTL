#!/usr/bin/env python3

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from src import base_dir, logger
from src.utils.regions import Interval


class TrackPlotter: 

	def __init__(self, bam_query, norm_factors, residualizer): 

		self.bam_query = bam_query
		self._norm_factors = norm_factors
		self._residualizer = residualizer

		self.sample_names = self.bam_query.sample_names

	@classmethod
	def load_atac(cls): 
		from src.analysis import atac_bam, atac_tmm_factors, atac_residualizer
		return cls(atac_bam, atac_tmm_factors, atac_residualizer)

	@classmethod
	def load_rna(cls): 
		from src.analysis import rna_bam, rna_tmm_factors, rna_residualizer
		return cls(rna_bam, rna_tmm_factors, rna_residualizer)

	def set_normalization(self, norm_factors, residualizer): 
		self._norm_factors = norm_factors
		self._residualizer = residualizer

		self._raw_pileups, self._norm_pileups, self._res_pileups = None, None, None

	def set_interval(self, interval, disjoint=True, downsample=True, max_positions=1e5, **kwargs): 

		logger.update("Fetching pileups...")
		interval = Interval.load(interval)
		self._raw_pileups = self.bam_query.query_pileup(interval, disjoint=disjoint) # numpy (2, n_bases, 181)
		self._norm_pileups = self._raw_pileups / self._norm_factors.values
		self._res_pileups = self._residualizer.transform(self._norm_pileups)

		self._pileup_index = range(interval.start, interval.end)

	def get_pileup(self, mode="res", strand=None, reformat=True, **kwargs): 
		"""
		Returns pileup dataframe.
		"""
		# Get pileups under various normalization approaches
		if mode == "raw": 
			pileup_data = self._raw_pileups
		elif mode == "norm": 
			pileup_data = self._norm_pileups
		elif mode == "res": 
			pileup_data = self._res_pileups
		else: 
			raise ValueError

		# Get specific strand of pileup
		if strand == "pos": 
			pileup_data = pileup_data[0]
		elif strand == "neg": 
			pileup_data = pileup_data[1]
		else: 
			pileup_data = pileup_data.sum(axis=0)

		# Cast as dataframe
		pileup_data = pd.DataFrame(pileup_data, index=self._pileup_index, columns=self.sample_names)

		if reformat: 
			pileup_data = self.reformat_pileups(pileup_data, **kwargs)

		return pileup_data

	def reformat_pileups(self, pileup_data, bounds=None, zero_pos=None, bin_size=None, n_bins=999, **kwargs): 
		"""
		Subsamples the track evenly for more manageable plotting. Can handle series or dataframe inputs.
		"""
		if bounds is None: bounds = (pileup_data.index.min(), pileup_data.index.max())

		# subset the pileups to region of interest
		pileup_data = pileup_data.loc[bounds[0]:bounds[1]+1]

		if bin_size is not None: n_bins = pileup_data.shape[1] // bin_size

		if pileup_data.shape[1] > n_bins: 
			# map original positions to bin labels
			bins_s = pd.Series(index=pileup_data.index, data=pd.cut(pileup_data.index, bins=n_bins, retbins=False))
			bins_s = bins_s.apply(lambda pos: pos.mid.round()).astype(int)

			# subsample pileups by taking mean of each bin
			pileup_data = pileup_data.groupby(bins_s).mean()

		# set index to center of bounds
		if zero_pos: 
			# midpoint = (bounds[0] + bounds[1]) // 2
			pileup_data.index = pileup_data.index - zero_pos

		return pileup_data


	def _basic_plotter(self, pileup_data_s, invert=False, reformat=True, orient="h", ax=None, **kwargs): 
		"""
		Plots a single track with given bounds. Must be provided as a series.
		"""
		# TODO: reformatting params?

		if invert: pileup_data_s *= -1

		plot_bounds = (pileup_data_s.index.min(), pileup_data_s.index.max())

		if orient == "h": 
			# plot track horizontally
			sns.lineplot(data=pileup_data_s, ax=ax, legend=False)
			ax.set(xlim=plot_bounds)
			if invert: 
				ax.axhline(0, color="black", alpha=0.5, linewidth=1)
		else: 
			# plot track vertically
			x,y = _rotate_axes(pileup_data_s.index, pileup_data_s.values)
			ax.plot(x,y)
			ax.set(ylim=plot_bounds)

	def plot_track(self, strand=None, feature="mean", **kwargs): 

		if strand == "both": 
			self.plot_track(strand="pos", feature=feature, **kwargs)
			self.plot_track(strand="neg", feature=feature, invert=True, **kwargs)

		else: 
			pileup = self.get_pileup(strand=strand, **kwargs)
			if feature == "mean": 
				pileup = pileup.mean(axis=1)
			elif feature == "norm_variance": 
				pileup = pileup.std(axis=1) / pileup.mean(axis=1)
			else: 
				# use a specific sample as the feature
				assert feature in pileup.columns
				pileup = pileup[feature]
			self._basic_plotter(pileup, reformat=False, **kwargs)

	def plot_by_genotype(self, genotypes, strand=None, **kwargs): 

		if strand == "both": 
			pass
		else: 
			pileup = self.get_pileup(strand=strand, **kwargs)
			grouped_pileups = pileup.groupby(genotypes, axis=1)



#----------------------------------------------------------------------------------------------------#
# Helpers
#----------------------------------------------------------------------------------------------------#

def _rotate_axes(x, y):
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



#----------------------------------------------------------------------------------------------------#
# Others
#----------------------------------------------------------------------------------------------------#
def plot_pileups_as_heatmap(pileups): 
	pass







#----------------------------------------------------------------------------------------------------#
# Old implementation
#----------------------------------------------------------------------------------------------------#
class _TrackPlotter: 

	def __init__(self, pos_tracks, neg_tracks, norm_factors=None, residualizer=None): 

		# Tracks should be library-normalized but not residualized yet
		self.pos_tracks = pos_tracks
		self.neg_tracks = neg_tracks
		self.tracks = self.pos_tracks + self.neg_tracks

		self.norm_factors = norm_factors
		self.residualizer = residualizer
		self.covariates = self.residualizer.C.T

		# self.positions = self.tracks.index

	# @classmethod
	# def load(cls, genomic_data, residualizer, interval_obj, **kwargs): 
	# 	chrom, start, end = interval_obj.coords
	# 	pos_tracks, neg_tracks = genomic_data.get_pileups(chrom, start, end, strand="+", report_strand="both", **kwargs)
	# 	return cls(pos_tracks, neg_tracks, residualizer)

	def scale(self, scale): 
		return TrackPlotter(self.pos_tracks * scale, self.neg_tracks * scale, self.norm_factors, self.residualizer)

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

	# @staticmethod
	# def normalize(track_data, depth=True, residualize=True): 
	# 	if depth: 
	# 		track_data /= self.norm_factors
	# 	if residualize: 
	# 		track_data = self.residualize.transform(track_data)
	# 	return track_data

	def get_tracks(self, strand=None, depth=True, residualize=True, coursen=True, **kwargs):
		if strand == "pos": 
			tracks = self.pos_tracks
		elif strand == "neg": 
			tracks = self.neg_tracks
		elif strand is None: 
			tracks = self.tracks

		if coursen: 
			tracks = TrackPlotter.coursen(tracks, **kwargs)

		if depth: 
			tracks /= self.norm_factors

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

	def plot_track(self, strand=None, feature="mean", **kwargs): 

		if strand == "both": 
			self.plot_track(strand="pos", feature=feature, **kwargs)
			self.plot_track(strand="neg", feature=feature, invert=True, **kwargs)
		else: 
			tracks = self.get_tracks(strand=strand, **kwargs)
			if feature == "mean": 
				track = tracks.mean(axis=1)
			elif feature == "variance_norm": 
				track = tracks.std(axis=1) / tracks.mean(axis=1)
			else: 
				track = tracks.iloc[:,feature]
			TrackPlotter._plotter(track, coursen=False, **kwargs)

	def plot_tracks_by_genotype(self, genotypes, strand=None, **kwargs): 

		if strand == "both":
			self.plot_tracks_by_genotype(genotypes=genotypes, strand="pos", **kwargs)
			self.plot_tracks_by_genotype(genotypes=genotypes, strand="neg", invert=True, **kwargs)
		else: 
			tracks = self.get_tracks(strand=strand, **kwargs)
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
			results = results.sort_values("phenotype_id").set_index("phenotype_id")["r"].fillna(0)
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
		positions = data.rsid.loc[variants, "pos"]
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
		qtl_results["pos"] = qtl_results["variant_id"].map(data.rsid["pos"])
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
			y.plot_track(**x_plot_params, strand="pos", ax=axes[0][1])
			x.plot_track(**x_plot_params, strand="pos", ax=axes[0][1])
			# neg track on bottom
			y.plot_track(**x_plot_params, strand="neg", ax=axes[2][1], invert=True)
			x.plot_track(**x_plot_params, strand="neg", ax=axes[2][1], invert=True)

			y.plot_track(**y_plot_params, strand="pos", orient="v", ax=axes[1][0])
			x.plot_track(**y_plot_params, strand="pos", orient="v", ax=axes[1][0])

			y.plot_track(**y_plot_params, strand="neg", orient="v", ax=axes[1][2], invert=True)
			x.plot_track(**y_plot_params, strand="neg", orient="v", ax=axes[1][2], invert=True)
			axes[1][2].set_yticklabels([])

			TrackPlotter.plot_correlation(x, y, x_plot_params, y_plot_params, r_max=r_max, ax=axes[1][1])