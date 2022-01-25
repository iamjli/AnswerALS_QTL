#!/usr/bin/env python3

import numpy as np
from pathos import multiprocessing


def chunk_list(l, n): 
	"""Yield successive n-sized chunks from lst."""
	for i in range(0, len(l), n):
	    yield l[i:i + n]


def multiprocess_by_chunk(func, chunk_arg, kwargs=dict(), n_cpus=24): 
	"""Useful when batch arg is large and would normally incur a large overhead."""
	arg_type = type(chunk_arg)
	if arg_type == dict: 
		keys, chunk_arg = zip(*chunk_arg.items())
		# chunk_arg = list(chunk_arg)
	else: 
		assert arg_type == list

	chunked_args = chunk_list(chunk_arg, n_cpus)
	chunk_func = lambda args: [func(arg, **kwargs) for arg in args] 

	with multiprocessing.ProcessingPool(n_cpus) as pool: 
		results = pool.map(chunk_func, chunked_args)
	results = [item for sublist in results for item in sublist]

	if arg_type == dict:
		return dict(zip(keys, results))
	else: 
		return results


	# batch_arg_split = np.array_split(batch_arg, n_cpus)

	# batch_func = lambda args: [func(arg, *other_args) for arg in args] 

	# with multiprocessing.ProcessingPool(n_cpus) as pool: 
	# 	results = pool.map(batch_func, batch_arg_split)
	# results = [item for sublist in results for item in sublist]

	# return results