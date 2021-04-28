from setuptools import find_packages, setup
from Cython.Build import cythonize

import numpy

import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

setup(
	name='src',
	packages=find_packages(),
    ext_modules=cythonize("src/query/pysam_utils.pyx", annotate=True),
    include_dirs=[numpy.get_include()],
    zip_safe=False,
	version='0.0.1',
	description='QTL analysis',
	author='iamjli',
	license='MIT',
)