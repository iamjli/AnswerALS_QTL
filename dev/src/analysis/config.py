#!/usr/bin/env python3

# from src.analysis.config import _omic_types, _genomic_features, _molecular_qtls


_feature_tag = ["atac", "rna", "erna", "srna"]

# Same omics type should share: counts normalization, residualizer, prefix
_omic_tag = {
	"atac": "atac", 
	"rna":  "rna", 
	"erna": "rna", 
	"srna": "rna",
}

# Same genomic feature type should share: regions metadata, 
_genomic_feature_tag = {
	"atac": "peak", 
	"rna":  "gene", 
	"erna": "peak", 
	"srna": "gene",
}

_mol_qtl_tag = {
	"atac": "caQTL", 
	"rna":  "eQTL", 
	"erna": "erQTL", 
	"srna": "sQTL",
}

# Which version of processed counts to use? 
_counts_prefix = {
	"atac": "atac_pysam_200bp",
	"rna":  "rna",
}