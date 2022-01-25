#!/usr/bin/env python3

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from src import logger


#----------------------------------------------------------------------------------------------------#
# DataFrames
#----------------------------------------------------------------------------------------------------#
def clip_extreme_values(df, lower=None, upper=None): 
	if upper is None: 
		upper = df.replace(np.inf, np.nan).max().max()
		df = df.replace(np.inf, upper)
	if lower is None: 
		if (df.values.min() < 0).any(): 
			lower = df.replace(-np.inf, np.nan).min().min()
			df = df.replace(-np.inf, lower)
		else: 
			lower = df.replace(0, np.nan).min().min()
			df = df.replace(0, lower)
	return df


def clip_extreme_values_by_axis(df, axis): 

	if axis == 0: 
		logger.write("Replacing with column min and max")
		# replace max
		max_values = df.replace()
	pass


#----------------------------------------------------------------------------------------------------#
# Heatmaps
#----------------------------------------------------------------------------------------------------#
def add_column_spacing(ax, axis, color="white", linewidth=1, offset=0): 

	if axis == "col":
		ticks = ax.get_xticks()
		for i in range(len(ticks) - 1): 
			pos = (ticks[i] + ticks[i+1]) / 2 + offset
			ax.axvline(pos, color=color, lw=linewidth)
	else: 
		ticks = ax.get_yticks()
		for i in range(len(ticks) - 1): 
			pos = (ticks[i] + ticks[i+1]) / 2 + offset
			ax.axhline(pos, color=color, lw=linewidth)


def _sized_heatmap_from_longform(data, x, y, hue, size, **kwargs): 

	pass


def map_continuous_colors(s, vmin=None, vmax=None, center=None, cmap='bwr'): 
    if center == 0: 
        abs_max = s.abs().max()
        vmin, vmax = -abs_max, abs_max
    else: 
        vmin, vmax = s.min(), s.max()
        
    n_colors = 256
    palette = sns.color_palette(cmap, n_colors=n_colors)
    bin_colors = {i:palette[i] for i in range(n_colors)}
    bin_colors[-1] = (0.5,0.5,0.5)
    # bin_colors = {i:tuple(x*256 for x in val) for i,val in bin_colors.items()}
    
    bin_edges = np.linspace(vmin, vmax, num=n_colors+1)
    # return bin_edges
    
    binned_vals = pd.cut(s, bins=bin_edges, labels=np.arange(n_colors), include_lowest=True).values.add_categories(-1).fillna(-1)
    binned_vals = pd.Series(binned_vals, index=s.index).astype(int)
    # return bin_colors
    # return binned_vals, bin_colors
    return binned_vals.map(bin_colors)



def heatmap_with_sizes(vals_df, size_df, s_scale, ax):


	# Mapping from column names to integer coordinates
	xlabels, ylabels = vals_df.columns, vals_df.index
	# xlabel_to_num = {label:i for i,label in enumerate(xlabels)}
	xlabel_to_idx = pd.Series(index=xlabels, data=range(len(xlabels)))
	ylabel_to_idx = pd.Series(index=ylabels, data=range(len(ylabels)))

	vals_df = vals_df.reset_index()
	vals_melted = vals_df.melt(id_vars=vals_df.columns[0])
	vals_melted.columns = ["y", "x", "vals"]

	colors = map_continuous_colors(vals_melted["vals"], center=0)

	size_df = size_df.reset_index()
	size_melted = size_df.melt(id_vars=size_df.columns[0])
	size_melted.columns = ["y", "x", "sizes"]

	ax.scatter(
		x=vals_melted["x"].map(xlabel_to_idx),
		y=vals_melted["y"].map(ylabel_to_idx),
		s=size_melted["sizes"] * s_scale,
		color=colors
	)
	sns.despine()
	# return sns.scatterplot(
	# 	x=vals_melted["x"].map(xlabel_to_idx),
	# 	y=vals_melted["y"].map(ylabel_to_idx),
	# 	color=colors, ax=ax
	# )


	
	# x_labels = [v for v in sorted(x.unique())]
	# y_labels = [v for v in sorted(y.unique())]
	# x_to_num = {p[1]:p[0] for p in enumerate(x_labels)} 
	# y_to_num = {p[1]:p[0] for p in enumerate(y_labels)} 
	
	# size_scale = 500
	# ax.scatter(
	# 	x=x.map(x_to_num), # Use mapping for x
	# 	y=y.map(y_to_num), # Use mapping for y
	# 	s=size * size_scale, # Vector of square sizes, proportional to size parameter
	# 	marker='s' # Use square as scatterplot marker
	# )
	
	# # Show column labels on the axes
	# ax.set_xticks([x_to_num[v] for v in x_labels])
	# ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='right')
	# ax.set_yticks([y_to_num[v] for v in y_labels])
	# ax.set_yticklabels(y_labels)



#----------------------------------------------------------------------------------------------------#
# Color mapping
#----------------------------------------------------------------------------------------------------#
def get_color_map(vals, palette=None): 
	"""
	Returns a dictionary mapping category to color. Expects a unique list of categories
	"""
	assert len(vals) == len(set(vals)), print("Values must be unique")
	
	if palette is None: palette = sns.color_palette()

	return {val:palette[i] for i,val in enumerate(vals)}

def map_series_to_colors(s, palette=None, sort=True): 
	"""
	Returns a series of colors.
	"""
	categories = s.unique()
	if sort: categories = np.sort(categories)

	cmap = get_color_map(categories, palette=palette)

	return s.map(cmap)

def display_cmap(color_dict): 
	
	n_colors = len(color_dict)
	cmap = mpl.colors.LinearSegmentedColormap.from_list("custom", sns.color_palette()[:n_colors], N=4)

def batch_cmaps(df): 
	cmaps = [ColorMap.load_data(df[col]) for col in df.columns]
	row_colors = pd.concat([cmap.mapped for cmap in cmaps], axis=1)
	return cmaps, row_colors

class ColorMap: 

	def __init__(self, labels, palette=None, name=None): 

		self.labels = labels
		self.palette = palette if palette is not None else sns.color_palette()
		self.name = name if name is not None else "custom"
		assert len(self.labels) <= len(self.palette)

		self.cmap = get_color_map(self.labels, self.palette)

	@classmethod
	def load_data(cls, data, order=None, **kwargs): 

		# if data.dtype == "category": 
		# 	try: 
		# 		data = pd.to_numeric(data)
		# 	except: 
		data = data.astype(object)

		if order is None: 
			labels = data.unique()
			labels = np.sort(labels)
		else: 
			labels = order

		cls_obj = cls(labels, name=data.name, **kwargs)
		cls_obj.data = data
		cls_obj.mapped = cls_obj.map_series(cls_obj.data)
		return cls_obj

	@property
	def colors(self):
		return [self.cmap[label] for label in self.labels]

	def map_series(self, s): 
		return s.map(self.cmap)

	def display_cmap(self): 

		N = len(self.labels)
		_cmap = mpl.colors.LinearSegmentedColormap.from_list(self.name, self.colors, N=N)
		_norm = mpl.colors.BoundaryNorm(np.arange(-0.5,N), _cmap.N)

		with plt.rc_context(dict(**sns.plotting_context("notebook", font_scale=1))):

			fig, ax = plt.subplots(figsize=(6, 1))
			fig.subplots_adjust(bottom=0.5)
			x = fig.colorbar(mpl.cm.ScalarMappable(cmap=_cmap, norm=_norm), ticks=np.linspace(0,N-1,N), cax=ax, orientation="horizontal")
			x.set_ticklabels(self.labels) 