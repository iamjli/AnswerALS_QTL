#!/usr/bin/env python3

__all__ = ["aals", "data", "hg38"]

# module-wide singletons for accessing data
from src.load.aals_data import aals
from src.load.external_data import data
from src.load.genome_data import hg38