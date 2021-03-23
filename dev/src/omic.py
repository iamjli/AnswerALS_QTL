#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scipy.stats as stats
import warnings

from pathlib import Path

from . import BASE_DIR, CHROMS, logger
from .bed import Regions


def load_omics(path, initialize=False): 

	df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")

	if df.index.name == "gene_id": 
		omic = "rna"
		regions = df[["chrom", "start", "end", "strand"]]
		lengths = df["lengths"]
		counts = df.iloc[:,5:]

	if df.index.name == "peak_id": 
		omic = "atac"
		regions = df[["chrom", "start", "end"]]
		lengths = df["lengths"]
		counts = df.iloc[:,4:]

	return Omic(omic, counts, regions, lengths, initialize)


class Omic: 

	def __init__(self, omic, counts, regions, lengths=None, initialize=False): 
		"""Accept either raw counts, or phenotype file."""

		assert (omic == "rna") or (omic == "atac")
		self.omic    = omic
		self.counts  = counts.copy()
		self.regions = regions.copy()

		# prepare regions and use for indexing
		if not self.regions.iloc[0,0].startswith("chr"): 
			self.regions.iloc[:,0] = "chr" + self.regions.iloc[:,0]
		self.regions.iloc[:,1:3] = self.regions.iloc[:,1:3].astype(int)

		if self.omic == "rna": 
			self.regions.columns = ["chrom", "start", "end", "strand"]
			self.regions.index.name = "gene_id"
		else: 
			self.regions.columns = ["chrom", "start", "end"]
			self.regions.index = self.regions["chrom"] + ":" + self.regions["start"].astype(str) + "-" + self.regions["end"].astype(str)
			self.regions.index.name = "peak_id"

		self.regions = Regions(self.regions)

		# compute lengths for TPM
		if lengths is None: 
			if self.omic == "rna":
				logger.write("Warning: lengths should be specified for RNA.")
			self.lengths = (self.regions["end"] - self.regions["start"]).astype(int).rename("lengths")
		else: 
			self.lengths = lengths.copy().rename("lengths")
		
		# get everything onto the same index
		self.counts.index  = self.regions.index
		self.lengths.index = self.regions.index

		# Normalization
		self.sample_frac_threshold = 0.2
		self.count_threshold = 6
		self.tpm_threshold = 0.1 
		self._mask = None

		self._lib_norm_factors = None
		self._tmm_norm_factors = None

		self._tpm = None
		self._tmm = None
		self._gtex = None

		if initialize: 
			self.tmm
			self.tpm
			self.gtex

	@property
	def mask(self):
		n_threshold = self.counts.shape[1] * self.sample_frac_threshold

		if self._mask is None: 
			if self.omic == "rna": 
				self._mask = (
					((self.counts >= self.count_threshold).sum(axis=1) >= n_threshold) & 
					((self.tpm >= self.tpm_threshold).sum(axis=1) >= n_threshold) & 
					(self.regions["chrom"].isin(CHROMS))
				)
			else: 
				self._mask = (
					((self.counts >= self.count_threshold).sum(axis=1) >= n_threshold)
				)
		return self._mask
	
	@property
	def tpm(self):
		if self._tpm is None: 
			tpm = self.counts.div(self.lengths, axis=0)
			self._tpm = tpm.div(tpm.sum() / 1e6, axis=1)
		return self._tpm
	
	@property
	def tmm(self):
		"""edgeR normalization: normalize by library size and TMM factors."""
		if self._tmm is None: 
			logger.write("Computing TMM matrix...")
			norm_factors = self.lib_norm_factors * self.tmm_norm_factors
			self._tmm = self.counts / norm_factors
		return self._tmm

	@property
	def gtex(self):
		if self._gtex is None: 
			self._gtex = inverse_normal_transform(self.tmm[self.mask])
		return self._gtex

	@property
	def tmm_norm_factors(self):
		if self._tmm_norm_factors is None: 
			self._tmm_norm_factors = edgeR_calcNormFactors(self.counts)
		return self._tmm_norm_factors

	@property
	def lib_norm_factors(self):
		if self._lib_norm_factors is None: 
			if self.omic == "rna": 
				self._lib_norm_factors = self.counts.sum(axis=0)
			elif self.omic == "atac": 
				self._lib_norm_factors = 1  # ATAC data has already been normalized by library size
		return self._lib_norm_factors

	def write_tensorqtl_phenotypes(self, output_dir): 

		output_dir = Path(output_dir)
		prefix = "{}_gtex".format(self.omic)

		df = pd.concat([self.regions, self.gtex], axis=1, join="inner")
		df = df.reset_index().rename(columns={"chrom": "#chr", self.regions.index.name:"gene_id"})
		df = df[["#chr", "start", "end", "gene_id"] + list(df.columns[4:])]
		df["#chr"] = df["#chr"].str[3:] # remove chromosome prefix

		# Write phenotype files for tensorqtl
		df.to_csv(output_dir / "{}.bed.gz".format(prefix), sep="\t", index=False, compression="gzip")

		# Write per-chromosome phenotype files for parallel processing
		for chrom in df["#chr"].unique(): 
			df_chrom = df[df["#chr"] == chrom]
			df_chrom.to_csv(output_dir / "{}.chr{}.bed.gz".format(prefix, chrom), sep="\t", index=False, compression="gzip")

		# Write phenotype file for gtex PEER pipeline. Remove strand column
		df.drop(columns=["strand"], errors="ignore").to_csv(output_dir / "{}.for_PEER.bed.gz".format(prefix), sep="\t", index=False, compression="gzip")

	def dump(self, output_dir): 
		"""Dump contents to load later."""
		output_dir = Path(output_dir)
		file_name = "_{}_counts.omics_dump.txt.gz".format(self.omic)

		df = pd.concat([
			self.regions, 
			self.lengths.rename("lengths"), 
			self.counts
		], axis=1)

		df.to_csv(output_dir / file_name, sep="\t", index=True, compression="gzip")
		logger.write("Dumped counts and metadata to", output_dir / file_name)







######################################
## 			NORMALIZATIONS			##
######################################

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

def inverse_normal_transform(M):
    """
    Transform rows to a standard normal distribution
    """
    R = stats.mstats.rankdata(M, axis=1)  # ties are averaged
    if isinstance(M, pd.DataFrame):
        Q = pd.DataFrame(stats.norm.ppf(R/(M.shape[1]+1)), index=M.index, columns=M.columns)
    else:
        Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q