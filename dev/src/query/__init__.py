#!/usr/bin/env python3

from src.load import aals
from src.query import bam, vcf

atac_bams = bam.QueryBams(omic="atac", bam_paths=aals.atac_bams)
rna_bams = bam.QueryBams(omic="rna", bam_paths=aals.rna_bams)

vcf = vcf.QueryVCF(vcf_path=aals.paths["vcf"])

# Validations
assert vcf.sample_names.equals(aals.sample_names), "VCF sample name mismatch"
