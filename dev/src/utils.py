#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy import interpolate, stats


flatten = lambda l: [item for sublist in l for item in sublist]

# CHROMS = ["chr"+str(i) for i in list(range(1,23))] + ["chrX", "chrY"]


## Mapping
# ENSG = load_ENSG_to_symbol()

# def symbol_to_ENSG(df): 
# 	symbols = pd.Series(dict((v,k) for k,v in ENSG.drop_duplicates().iteritems()))
# 	return df.groupby(symbols.reindex(df.index)).mean()


# 	# ENSGs = ENSG[ENSG == symbol]
# 	# if len(ENSGs) != 1: print("{} entries found for {}".format(symbol))
# 	# return ENSGs.index[0]

# def ENSG_to_symbol(df): 
# 	return df.groupby(ENSG.reindex(df.index)).mean()


=







