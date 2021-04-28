#!/usr/bin/env python3

from src.load import aals
from src.query import bam, vcf

atac_bam = bam.QueryBam(bam_files=aals.atac_bams)
rna_bam = bam.QueryBam(bam_files=aals.rna_bams)

vcf = vcf.QueryVCF(vcf_path=aals.paths["vcf"])

# Validations
assert vcf.sample_names.equals(aals.sample_names), "VCF sample name mismatch"
