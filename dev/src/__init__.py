#!/usr/bin/env python3

from src.utils import set_imports
set_imports()


from pathlib import Path
base_dir = Path("/nfs/latdata/iamjli/ALS/analysis/210211_250_lines")


from src.logger import logger

# module-wide singletons
from src.data_repo.aals_data import aals
from src.data_repo.external_data import data
from src.data_repo.genome_data import hg38