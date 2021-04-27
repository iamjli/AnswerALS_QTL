#!/usr/bin/env python3

import numpy as np
import pandas as pd

from pathlib import Path
import warnings

__all__ = ["counts_to_tpm", "inverse_normal_transform", "normalize_quantiles", "edgeR_calcNormFactors"]

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
