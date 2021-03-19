

# 3/5/21

ATAC diffBind, counts normalization and PCA
http://localhost:8891/notebooks/210211_250_lines/atac/diffbind3/prep_diffbind3.ipynb


Harmonization of sample IDs, generation of FastQTL inputs, various normalizations for RNA & ATAC
http://localhost:8894/notebooks/prep_QTLs/harmonization.ipynb


Commands for tensorQTL
http://localhost:8894/notebooks/prep_QTLs/tensorqtl/run_tensorqtl.ipynb


Parse tensorQTL results, adjusted p-values
http://localhost:8894/notebooks/prep_QTLs/tensorqtl/parse_tensorqtl.ipynb






# 3/13/21 - VQSR

```bash

bcftools=/nfs/apps/bcftools-1.11/bin/bcftools
wd=/globus/genomics_base/Sentieon_JG_Fall_2020


# 1a. split into individual chromosomes: `./final_merged_individual_chromosomes`
$bcftools view $wd/final_merged_chromosomes/final_chr4_and_chr5.vcf.gz -r chr4 -o $wd/final_merged_individual_chromosomes/final_chr4.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr6_and_chr7.vcf.gz -r chr6 -o $wd/final_merged_individual_chromosomes/final_chr6.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr8_and_chr9.vcf.gz -r chr8 -o $wd/final_merged_individual_chromosomes/final_chr8.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr10_and_chr11.vcf.gz -r chr10 -o $wd/final_merged_individual_chromosomes/final_chr10.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr12_and_chr13.vcf.gz -r chr12 -o $wd/final_merged_individual_chromosomes/final_chr12.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr14_and_chr15.vcf.gz -r chr14 -o $wd/final_merged_individual_chromosomes/final_chr14.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr16_and_chr17.vcf.gz -r chr16 -o $wd/final_merged_individual_chromosomes/final_chr16.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr18_and_chr19.vcf.gz -r chr18 -o $wd/final_merged_individual_chromosomes/final_chr18.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr20_and_chr21.vcf.gz -r chr20 -o $wd/final_merged_individual_chromosomes/final_chr20.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr22_and_X_and_Y.vcf.gz -r chr22 -o $wd/final_merged_individual_chromosomes/final_chr22.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr22_and_X_and_Y.vcf.gz -r chrX -o $wd/final_merged_individual_chromosomes/final_chrX.vcf.gz -O z --threads 12 &

$bcftools view $wd/final_merged_chromosomes/final_chr4_and_chr5.vcf.gz -r chr5 -o $wd/final_merged_individual_chromosomes/final_chr5.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr6_and_chr7.vcf.gz -r chr7 -o $wd/final_merged_individual_chromosomes/final_chr7.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr8_and_chr9.vcf.gz -r chr9 -o $wd/final_merged_individual_chromosomes/final_chr9.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr10_and_chr11.vcf.gz -r chr11 -o $wd/final_merged_individual_chromosomes/final_chr11.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr12_and_chr13.vcf.gz -r chr13 -o $wd/final_merged_individual_chromosomes/final_chr13.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr14_and_chr15.vcf.gz -r chr15 -o $wd/final_merged_individual_chromosomes/final_chr15.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr16_and_chr17.vcf.gz -r chr17 -o $wd/final_merged_individual_chromosomes/final_chr17.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr18_and_chr19.vcf.gz -r chr19 -o $wd/final_merged_individual_chromosomes/final_chr19.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr20_and_chr21.vcf.gz -r chr21 -o $wd/final_merged_individual_chromosomes/final_chr21.vcf.gz -O z --threads 12 &
$bcftools view $wd/final_merged_chromosomes/final_chr22_and_X_and_Y.vcf.gz -r chrY -o $wd/final_merged_individual_chromosomes/final_chrY.vcf.gz -O z --threads 12 &

# 1b. index
for f in final*.vcf.gz; 
do 
	tabix -p vcf $f
done



# 2. Create small files with just genotype field
for raw_vcf in final*.vcf.gz; 
do 
	gt_only_vcf="GT_only.$(basename "$raw_vcf" .vcf.gz).vcf.gz"
	cmd="$bcftools annotate -x FORMAT $raw_vcf -o $gt_only_vcf -O z --threads 12 && tabix -p vcf $gt_only_vcf &"
	echo $cmd
done



# 3. Generate full list of SNPs. Output: `/globus/genomics_base/Sentieon_JG_Fall_2020/VQSR_210311/GVCFtyper_annotations.vcf.gz`




# 4a. Run VQSR, apply to `VQSR/GVCFtyper_annotations.vcf.gz`. This might be a good opportunity to annotate with rsIDs. 
bash $wd/VQSR/210225_varcal/run_vqsr.sh

# 4b. Add annotations (rsID)
rsid_anno_file=/nfs/latdata/iamjli/ALS/analysis/210211_250_lines/data/external/rsID_hg38/All_20180418.chr_anno.vcf.gz
vqsr_output=$wd/VQSR_210311/varcal/vqsrSNPINDEL.hc.sensitivity_99.recaled.vcf.gz
rsID_vqsr_output=$wd/VQSR_210311/varcal/vqsrSNPINDEL.hc.sensitivity_99.recaled.rsID.vcf.gz
$bcftools annotate -a $rsid_anno_file -c ID -o $rsID_vqsr_output -O z $vqsr_output --threads 12 && tabix -p vcf $rsID_vqsr_output

# 4c. Split output file by chromosome for later speed increases
tabix -l $rsID_vqsr_output | while read -r chrom; 
do 
	out_filename="${rsID_vqsr_output%.vcf.gz}.$chrom.vcf.gz"
	cmd="$bcftools view $rsID_vqsr_output -r $chrom -o $out_filename -O z --threads 12 && tabix -p vcf $out_filename"
	echo $cmd
done



# 5. Update ID,FILTER columns with rsID,VQSR status (PASS/tranche), then report only variants that passed. 
for filename in $wd/final_merged_individual_chromosomes/GT_only.final_chr*.vcf.gz; 
do
	base=$(basename "$filename" .vcf.gz) # get basename
	chrom=${base#"GT_only.final_"} # get chromosome
	anno_file=$wd/VQSR_210311/varcal/vqsrSNPINDEL.hc.sensitivity_99.recaled.rsID.${chrom}.vcf.gz
	out_filename="$wd/final_merged_individual_chromosomes/VQSR_filtered_99.rsID.$base.vcf.gz"

	cmd="""
	$bcftools annotate -a $anno_file -c ID,FILTER -r $chrom $filename -O u | \
	$bcftools view -f PASS -O u | \
	$bcftools annotate -x FILTER -o $out_filename -O z --threads 12 && \
	tabix -p vcf $out_filename"""
	echo $cmd
done
```


## Generate variant list
```bash
wd=/nfs/latdata/iamjli/ALS/analysis/210211_250_lines
out=$wd/tensorqtl/genomes/snp_list.biallelic.harmonized.VQSR_filtered_99.rsID.txt
echo -e "chrom\tpos\tvariant_id\tref\talt" > $out

while read chrom; 
do
    vcf=$wd/tensorqtl/genomes/biallelic.harmonized.VQSR_filtered_99.rsID.${chrom}.vcf.gz
    bcftools view -H $vcf | cut -f 1-5 >> $out
done < $wd/tensorqtl/genomes/chrom_list.txt
```
```python
snp_df = pd.read_csv(BASE_DIR / "tensorqtl/genomes/snp_list.biallelic.harmonized.VQSR_filtered_99.rsID.txt", sep="\t")
snp_df = snp_df[snp_df.variant_id != "."].set_index("variant_id")
snp_df.to_parquet(BASE_DIR / "tensorqtl/genomes/snp_list_filtered.biallelic.harmonized.VQSR_filtered_99.rsID.parquet", engine="pyarrow")
```

## Recode Project MINE
```python
# positive b means enriched in ALS
projectmine = pd.read_csv(BASE_DIR / "data/external/Summary_Statistics_GWAS_2016/als.sumstats.lmm.txt", sep=' ', index_col=1, low_memory=False)
projectmine.to_parquet(BASE_DIR / "data/external/Summary_Statistics_GWAS_2016/als.sumstats.lmm.parquet")
projectmine.head()
```
