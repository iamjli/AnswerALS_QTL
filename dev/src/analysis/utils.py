#!/usr/bin/env python3

import numpy as np
import pandas as pd

from src import logger


def unify_dfs(df1, df2): 

	assert df1.index.is_unique and df2.index.is_unique

	common_idx = df1.index & df2.index
	logger.write("Number of rows (df1, df2, common): {}, {}, {}".format(len(df1), len(df2), len(common_idx)))

	return df1.reindex(common_idx), df2.reindex(common_idx)

def outer_product(s1, s2): 

	return pd.DataFrame(np.outer(s2, s1), index = s2.index, columns = s1.index)


def explode(df, lst_cols, fill_value='', preserve_index=False):
	# make sure `lst_cols` is list-alike
	if (lst_cols is not None
		and len(lst_cols) > 0
		and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
		lst_cols = [lst_cols]
	# all columns except `lst_cols`
	idx_cols = df.columns.difference(lst_cols)
	# calculate lengths of lists
	lens = df[lst_cols[0]].str.len()
	# preserve original index values	
	idx = np.repeat(df.index.values, lens)
	# create "exploded" DF
	res = (pd.DataFrame({
				col:np.repeat(df[col].values, lens)
				for col in idx_cols},
				index=idx)
			 .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
							for col in lst_cols}))
	# append those rows that have empty lists
	if (lens == 0).any():
		# at least one list in cells is empty
		res = (res.append(df.loc[lens==0, idx_cols], sort=False)
				  .fillna(fill_value))
	# revert the original index order
	res = res.sort_index()
	# reset index if requested
	if not preserve_index:		
		res = res.reset_index(drop=True)
	return res
