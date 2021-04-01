#!/usr/bin/env python3
from pathlib import Path

BASE_DIR = Path("/nfs/latdata/iamjli/ALS/analysis/210211_250_lines")
CHROMS = ["chr"+str(i) for i in list(range(1,23))] + ["chrX", "chrY"]

from .utils import SimpleLogger
logger = SimpleLogger()

from .repository import Repository, RESULTS_PATHS 
DATA = Repository(BASE_DIR)

from .omic import Omic
from .bed import Regions
from .genomic import BamQuery, TrackPlotter, Interval
from .qtl import QTL

# from .bed import pr_to_df, get_point, join_regions_by_window, phenotype_snp_distance
# from .annotate import get_variant_annos_as_df, annotation_counts_per_region

from .graph import Forest, get_node_attributes_as_df, get_edge_attributes_as_df


# sns.set(font='arial', font_scale=1.5)