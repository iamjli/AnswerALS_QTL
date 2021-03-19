#!/usr/bin/env python3
from pathlib import Path

BASE_DIR = Path("/nfs/latdata/iamjli/ALS/analysis/210211_250_lines")
CHROMS = ["chr"+str(i) for i in list(range(1,23))] + ["chrX", "chrY"]

from .repository import Repository
from .omic import Omic
from .genomic import Genomic
# from .data import Genomic, Omic, External

# from .utils import symbol_to_ENSG, ENSG_to_symbol

# from .plot import get_coverage, get_snps

from .bed import df_to_pr, pr_to_df, get_point, join_regions_by_window, phenotype_snp_distance
from .annotate import get_variant_annos_as_df, annotation_counts_per_region
from .qtl import *


from .graph import Forest, get_node_attributes_as_df, get_edge_attributes_as_df



# sns.set(font='arial', font_scale=1.5)