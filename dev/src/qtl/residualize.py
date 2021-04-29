#!/usr/bin/env python3

import numpy as np
import pandas as pd
import torch

from src import logger, base_dir


class Residualizer(object):
	"""
	Based on `tensorqtl.core.Residualizer` but added support for dataframes.
	"""
	def __init__(self, C):
		"""
		C: samples x covariates
		"""
		self.C = C
		if isinstance(C, pd.DataFrame): 
			if not C.index.str.startswith("NEU").all(): 
				logger.write("Warning: check that input is indexed by sample name")
			C_t = torch.tensor(C.values, dtype=torch.float32).to("cpu") 
		elif isinstance(C, torch.Tensor): 
			C_t = C
		else: 
			logger.write("Must provide as dataframe or tensor.")

		# center and orthogonalize
		self.Q_t, _ = torch.qr(C_t - C_t.mean(0))
		self.dof = C_t.shape[0] - 2 - C_t.shape[1]

		self.n_samples = self.C.shape[0]

	@classmethod
	def load_covariates(cls, path=None, prefix=None): 
		if path is None: 
			path = base_dir / "tensorqtl_runs/_phenotypes" / f"{prefix}_gtex.PEER_10.condition_sex.PEER_covariates.txt"
		covariates_df = pd.read_csv(path, sep="\t", index_col=0)
		return cls(covariates_df.T)

	def transform(self, M, center=True):
		"""Residualize rows of M wrt columns of C. Does not necessarily need to be normalized."""

		if isinstance(M, pd.DataFrame): 
			input_format = "dataframe"
			M_t = torch.tensor(M.values, dtype=torch.float).to("cpu")
		elif isinstance(M, np.ndarray): 
			input_format = "array"

			M_input_shape = M.shape
			if len(M_input_shape) > 2: 
				assert M_input_shape[-1] == self.n_samples # check that last axis in M corresponds to samples
				M = M.reshape((-1, self.n_samples)) # stack along the first axes

			M_t = torch.tensor(M, dtype=torch.float).to("cpu")
		else: 
			input_format = "tensor"
			M_t = M

		# center row means
		M0_t = M_t - M_t.mean(1, keepdim=True) if center else M_t

		# the second term is the components of M that are explainable by Q. First projects into covariate space, then projects back
		# Note that normally the projection back would be Q_inverse, but because it's orthonormal, that's equal to Q^T
		M_t_transformed = M_t - torch.mm(torch.mm(M0_t, self.Q_t), self.Q_t.t())  # keep original mean

		if input_format == "dataframe": 
			return pd.DataFrame(M_t_transformed.numpy(), index=M.index, columns=M.columns)
		elif input_format == "array": 
			M_t_transformed = M_t_transformed.numpy()
			# return M_t_transformed
			if len(M_input_shape) > 2: 
				M_t_transformed = M_t_transformed.reshape(M_input_shape)
			return M_t_transformed
		else: 
			return M_t_transformed