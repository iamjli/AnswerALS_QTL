#!/usr/bin/env python3

import numpy as np
import pandas as pd

from pathlib import Path
import pickle
import warnings

from src import logger, aals, hg38
from src.regions import RegionsAccessor



class CountsData: 
	"""
	Methods for uniform processing of omics data.
	"""
	chroms = [f"chr{i}" for i in range(1,23)] + ["chrX", "chrY"]

	def __init__(self, counts, regions, lengths, sample_frac_threshold=0.2, count_threshold=6, tpm_threshold=0.1): 
		
		self.counts = counts.copy()
		self.regions = regions.copy()
		self.lengths = lengths.copy()

		self.sample_names = self.counts.columns

		# Library norm factors 
		self._cpm_factors = None
		self._tmm_factors = None

		# Counts matrices
		self._tpm = None
		self._tmm = None
		self._gtex = None

		# Thresholding parameters
		self.sample_frac_threshold = sample_frac_threshold
		self.count_threshold = count_threshold
		self.tpm_threshold = tpm_threshold
		self._mask = None

	def validate(self): 
		assert self.counts.index.isin(self.regions.index).all()
		assert self.sample_names.equals(aals.sample_names)

	def initialize(self): 
		_ = self.tpm, self.tmm, self.gtex, self.mask

	@property
	def mask(self):
		"""Filter low read counts."""
		if self._mask is None: 
			n_threshold = self.counts.shape[1] * self.sample_frac_threshold
			self._mask = (
				((self.counts >= self.count_threshold).sum(axis=1) >= n_threshold) & 
				((self.tpm >= self.tpm_threshold).sum(axis=1) >= n_threshold) & 
				(self.regions["chrom"].isin(self.chroms))
			)
		return self._mask

	#----------------------------------------------------------------------------------------------------#
	# Library size normalization factors
	#----------------------------------------------------------------------------------------------------#
	@property
	def cpm_factors(self):
		if self._cpm_factors is None:
			cpm = self.counts.sum(axis=0) / 1e6
			self._cpm_factors = cpm.rename("cpm")
		return self._cpm_factors

	@property
	def tmm_factors(self):
		if self._tmm_factors is None: 
			logger.write("Computing TMM factors...")
			# Abusing nomenclature. Technically, effective library size is TMM * CPM
			tmm_factors = edgeR_calcNormFactors(self.counts) * self.cpm_factors
			self._tmm_factors = pd.Series(data=tmm_factors, index=self.sample_names, name="tmm")
		return self._tmm_factors

	#----------------------------------------------------------------------------------------------------#
	# Normalized counts matrices
	#----------------------------------------------------------------------------------------------------#
	@property
	def tpm(self):
		if self._tpm is None: 
			self._tpm = counts_to_tpm(self.counts, self.lengths)
			rpk = self.counts.div(self.lengths / 1000, axis=0)
			self._tpm = rpk.div(rpk.sum() / 1e6, axis=1)
		return self._tpm
	
	@property
	def tmm(self):
		"""edgeR normalization: normalize by library size (counts per million) and TMM factors."""
		if self._tmm is None: 
			self._tmm = self.counts / self.tmm_factors
		return self._tmm

	@property
	def gtex(self):
		if self._gtex is None: 
			logger.write("Computing GTEX matrix...")
			self._gtex = inverse_normal_transform(self.tmm[self.mask])
		return self._gtex

	#----------------------------------------------------------------------------------------------------#
	# Save/load
	#----------------------------------------------------------------------------------------------------#
	def to_pickle(self, output_dir, prefix): 

		self.initialize()
		self.prefix = prefix
		pickle_path = Path(output_dir) / f"{self.prefix}.pkl"
		with open(pickle_path, "wb") as f:
			pickle.dump(self, f)

	@classmethod
	def load_pickle(cls, output_dir, prefix):

		pickle_path = Path(output_dir) / f"{prefix}.pkl"
		with open(pickle_path, "rb") as f: 
			return pickle.load(f)


#----------------------------------------------------------------------------------------------------#
# DataFrame wrapper for normalization
#----------------------------------------------------------------------------------------------------#
@pd.api.extensions.register_dataframe_accessor("norm")
class NormAccessor:

	def __init__(self, counts): 

		self._counts = counts

		self._inv_norm = None
		self._q_norm = None
		self._res = None

	@property
	def inverse_norm(self):
		if self._inv_norm is None: 
			self._inv_norm = inverse_normal_transform(self._counts)
		return self._inv_norm

	@property
	def quantile_norm(self):
		if self._q_norm is None: 
			self._q_norm = quantile_normalize(self._counts)
		return self._q_norm

	def residualize(self, residualizer):
		return residualizer.transform(self._counts)
	


#----------------------------------------------------------------------------------------------------#
# Normalization functions
#----------------------------------------------------------------------------------------------------#
def counts_to_tpm(counts, lengths): 
	rpk = counts.div(lengths / 1000, axis=0)
	return rpk.div(rpk.sum() / 1e6, axis=1)


def inverse_normal_transform(M):
	"""
	Transform rows to a standard normal distribution
	"""
	from scipy import stats

	if isinstance(M, pd.Series): 
		M_ = M.to_frame().T
	else: 
		M_ = M

	R = stats.mstats.rankdata(M_, axis=1)  # ties are averaged
	Q = stats.norm.ppf(R/(M_.shape[1]+1))

	if isinstance(M, pd.DataFrame):
		return pd.DataFrame(Q, index=M.index, columns=M.columns)
	elif isinstance(M, pd.Series): 
		return pd.Series(Q.squeeze(), index=M.index)
	else:
		Q = stats.norm.ppf(R/(M_.shape[1]+1))
	return Q

def normalize_quantiles(df):
	"""
	Quantile normalization to the average empirical distribution
	Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")
	Reference:
	 [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003
	Adapted from https://github.com/andrewdyates/quantile_normalize
	"""
	M = df.values.copy()

	Q = M.argsort(axis=0)
	m,n = M.shape

	# compute quantile vector
	quantiles = np.zeros(m)
	for i in range(n):
		quantiles += M[Q[:,i],i]
	quantiles = quantiles / n

	for i in range(n):
		# Get equivalence classes; unique values == 0
		dupes = np.zeros(m, dtype=np.int)
		for j in range(m-1):
			if M[Q[j,i],i]==M[Q[j+1,i],i]:
				dupes[j+1] = dupes[j]+1

		# Replace column with quantile ranks
		M[Q[:,i],i] = quantiles

		# Average together equivalence classes
		j = m-1
		while j >= 0:
			if dupes[j] == 0:
				j -= 1
			else:
				idxs = Q[j-dupes[j]:j+1,i]
				M[idxs,i] = np.median(M[idxs,i])
				j -= 1 + dupes[j]
		assert j == -1

	return pd.DataFrame(M, index=df.index, columns=df.columns)



def edgeR_calcNormFactors(counts, ref=None, logratio_trim=0.3, sum_trim=0.05, acutoff=-1e10, verbose=False):
	"""
	Calculate TMM (Trimmed Mean of M values) normalization.
	Reproduces edgeR::calcNormFactors.default
	Scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes.
	Effective library size: TMM scaling factor * library size
	References:
	 [1] Robinson & Oshlack, 2010
	 [2] R functions:
		  edgeR::calcNormFactors.default
		  edgeR:::.calcFactorWeighted
		  edgeR:::.calcFactorQuantile
	"""
	from scipy import stats

	# discard genes with all-zero counts
	Y = counts.values.copy()
	allzero = np.sum(Y>0,axis=1)==0
	if np.any(allzero):
		Y = Y[~allzero,:]

	# select reference sample
	if ref is None:  # reference sample index
		# f75 = np.percentile(Y/np.sum(Y,axis=0), 75, axis=0)
		f75 = np.array([np.percentile(x,75) for x in (Y/np.sum(Y,axis=0)).T])
		ref = np.argmin(np.abs(f75-np.mean(f75)))
		if verbose:
			logger.write('Reference sample index: '+str(ref))

	N = np.sum(Y, axis=0)  # total reads in each library

	# with np.errstate(divide='ignore'):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		logR = np.log2((Y/N).T / (Y[:,ref]/N[ref])).T  # log fold change; Mg in [1]
		absE = 0.5*(np.log2(Y/N).T + np.log2(Y[:,ref]/N[ref])).T  # average log relative expression; Ag in [1]
		v = (N-Y)/N/Y
		v = (v.T + v[:,ref]).T  # w in [1]

	ns = Y.shape[1]
	tmm = np.zeros(ns)
	for i in range(ns):
		fin = np.isfinite(logR[:,i]) & np.isfinite(absE[:,i]) & (absE[:,i] > acutoff)
		n = np.sum(fin)

		loL = np.floor(n*logratio_trim)+1
		hiL = n + 1 - loL
		loS = np.floor(n*sum_trim)+1
		hiS = n + 1 - loS
		rankR = stats.rankdata(logR[fin,i])
		rankE = stats.rankdata(absE[fin,i])
		keep = (rankR >= loL) & (rankR <= hiL) & (rankE >= loS) & (rankE <= hiS)
		# in [1], w erroneously defined as 1/v ?
		tmm[i] = 2**(np.nansum(logR[fin,i][keep]/v[fin,i][keep]) / np.nansum(1/v[fin,i][keep]))

	tmm = tmm / np.exp(np.mean(np.log(tmm)))
	return tmm



