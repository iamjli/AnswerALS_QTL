#!/usr/bin/env python3

# These need to be imported at the base level in this specific order 
# to avoid static TLS import errors
import torch
import sklearn
import scipy

from pathlib import Path
base_dir = Path("/nfs/latdata/iamjli/ALS/analysis/210211_250_lines")

from src.utils.logger import logger