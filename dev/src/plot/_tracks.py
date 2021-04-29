
class TrackPlotter: 

	def __init__(self, pos_tracks, neg_tracks, norm_factors=None, residualizer=None): 

		# Tracks should be library-normalized but not residualized yet
		self.pos_tracks = pos_tracks
		self.neg_tracks = neg_tracks
		self.tracks = self.pos_tracks + self.neg_tracks

		self.norm_factors = norm_factors
		self.residualizer = residualizer
		self.covariates = self.residualizer.C.T

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

