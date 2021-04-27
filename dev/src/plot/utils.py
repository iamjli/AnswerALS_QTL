#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
def add_column_spacing(ax, color="white", linewidth=1, offset=0): 

	ticks = ax.get_xticks()
	for i in range(len(ticks) - 1): 
		pos = (ticks[i] + ticks[i+1]) / 2 + offset
		ax.axvline(pos, color=color, lw=linewidth)

def sized_heatmap():
	pass


def _sized_heatmap_from_longform(data, x, y, hue, size, **kwargs): 

	pass



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

		labels = data.unique()
		if order is None: labels = np.sort(labels)

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