#!/usr/bin/env python3

from pathlib import Path

base_dir = Path("/nfs/latdata/iamjli/ALS/analysis/210211_250_lines")

from src.logger import logger

# module-wide singletons
from src.data_repo.reference_genomes import hg38
from src.data_repo.aals_data import aals
from src.data_repo.external_data import data