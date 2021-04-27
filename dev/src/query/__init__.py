#!/usr/bin/env python3


from src import aals
from src.query.bam import QueryBam
from src.query.vcf import QueryVCF


atac_bam = QueryBam(bam_files=aals.atac_bams)
rna_bam = QueryBam(bam_files=aals.rna_bams)

vcf = QueryVCF(vcf_path=aals.paths["vcf"])

# Validations
assert vcf.sample_names.equals(aals.sample_names), "VCF sample name mismatch"
