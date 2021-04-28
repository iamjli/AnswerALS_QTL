#!/usr/bin/env python3

import sys
import traceback 


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

	def write(self, *message, verbose=True):
		message = " ".join([str(m) for m in message])
		if verbose:
			self.console.write(message+"\n")
		if self.log is not None:
			self.log.write(message+'\n')
			self.log.flush()

	def update(self, *message, verbose=True): 
		message = " ".join([str(m) for m in message])
		if verbose: 
			self.console.write("\r{}".format(message))
			self.console.flush()
		if self.log is not None: 
			self.log.write(message+'\n')
			self.log.flush()

	def flush(self, verbose=True): 
		if verbose: 
			self.console.write('\n')

logger = SimpleLogger()