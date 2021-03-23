#!/usr/bin/env python3
import sys
import traceback 

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

def print_update(s): 

	sys.stdout.write("\r{}".format(s))
	sys.stdout.flush()

class SimpleLogger(object):
	def __init__(self, logfile=None, verbose=True):
		self.console = sys.stdout
		self.verbose = verbose
		if logfile is not None:
			self.log = open(logfile, 'w')
		else:
			self.log = None

	def __enter__(self): 
		return self.update

	def __exit__(self, exc_type, exc_value, tb):
		if exc_type is not None:
			traceback.print_exception(exc_type, exc_value, tb)
			return False # uncomment to pass exception through
		# self.flush()
		return True

	def write(self, message):
		if self.verbose:
			self.console.write(message+"\n")
		if self.log is not None:
			self.log.write(message+'\n')
			self.log.flush()

	def update(self, message): 
		if self.verbose: 
			self.console.write("\r{}".format(message))
			self.console.flush()
		if self.log is not None: 
			self.log.write(message+'\n')
			self.log.flush()

	def flush(self): 
		if self.verbose: 
			self.console.write('\n')






