#!/usr/bin/env python3

import numpy as np
import pandas as pd

from src import logger


def unify_dfs(df1, df2): 

	assert df1.index.is_unique and df2.index.is_unique

	common_idx = df1.index & df2.index
	logger.write("Number of rows (df1, df2, common): {}, {}, {}".format(len(df1), len(df2), len(common_idx)))

	return df1.reindex(common_idx), df2.reindex(common_idx)
