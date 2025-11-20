#!/usr/bin/env bash

#SBATCH --job-name=gwas_qc_struct
#SBATCH --partition=msibigmem,msilarge,msismall
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=balay011@umn.edu

set -euo pipefail
echo "[$(date)] Job started on ${HOSTNAME} with ${SLURM_CPUS_PER_TASK:-unknown} CPUs"

######################## This script takes directories with input vcf SNP files and runs:
######################## I) QC on all provided germline SNPs
######################## II) Processing of 1KG phase 3 SNPs with vcf SNPs to make PCs
######################## III) Inference of per-sample genome-wide global ancestry

# Threads from SLURM (fallback 8)
THREADS=${SLURM_CPUS_PER_TASK:-8}

THREADS=8
### Input vcf SNP directory
SNPS_DIR="/scratch.global/balay011/germline_calls/vqsr/sample_stats/final_merged/snps"

### Tools/modules
module load bcftools/1.16-gcc-8.2.0-5d4xg4y
module load conda
source activate gwas_env

BCFTOOLS="bcftools"
PLINK2="plink2"
PLINK="plink"
ADMIXTURE="admixture"

### Output directory
OUTDIR="/scratch.global/balay011/gwas_qc_struct"
mkdir -p "$OUTDIR"/{logs,tmp,ref,admixture,pca}
cd "$OUTDIR"

######################## I) QC on all provided germline SNPs
# QC thresholds
MIND=0.1 # sample missingness
GENO=0.05 # variant missingness
MAF_CUTOFF=0.01 # MAF cutoff
KING_CUTOFF=0.354 # duplicate cutoff

# LD pruning
PRUNE_WINDOW_KB=200
PRUNE_STEP=50
PRUNE_R2=0.1

# ADMIXTURE
POPS=(AFR AMR EAS EUR SAS)
ADMIXED_SECONDARY_THRESHOLD=0.20

echo ">>> Starting pipeline at $(date). Output: $OUTDIR"

MERGED_FINAL_SNPS_VCF="/scratch.global/balay011/germline_calls/vqsr/sample_stats/final_merged/cohort.snps.merged.vcf.gz"

### Convert merged VCF to PLINK2 PGEN by acquiring autosomal biallelic SNPs
echo ">>> Step 1: Import autosomal biallelic SNPs directly to PGEN"
$PLINK2 --threads $THREADS \
  --vcf $MERGED_FINAL_SNPS_VCF \
  --autosome \
  --snps-only just-acgt \
  --max-alleles 2 \
  --set-all-var-ids @:#:\$r:\$a \
  --rm-dup force-first \
  --make-pgen --out cohort1_auto
# --vcf: 23062756 variants scanned (745012 skipped).
# --vcf: cohort1_auto-temporary.pgen + cohort1_auto-temporary.pvar.zst +
# cohort1_auto-temporary.psam written.
# 215 samples (0 females, 0 males, 215 ambiguous; 215 founders) loaded from
# cohort1_auto-temporary.psam.
# 22919404 out of 23062756 variants loaded from cohort1_auto-temporary.pvar.zst.
# Note: No phenotype data present.
# Note: Skipping --rm-dup since no duplicate IDs are present.
# 22919404 variants remaining after main filters.
# Writing cohort1_auto.psam ... done.
# Writing cohort1_auto.pvar ... done.
# Writing cohort1_auto.pgen ... done.

### Missingness, MAF, relatedness
echo ">>> Step 2: Filtering (mind=$MIND, geno=$GENO, maf=$MAF_CUTOFF), KING>$KING_CUTOFF"
$PLINK2 --threads $THREADS --pfile cohort1_auto \
  --geno $GENO \
  --make-pgen --out cohort2_auto
# 215 samples (0 females, 0 males, 215 ambiguous; 215 founders) loaded from
# cohort1_auto.psam.
# 22919404 variants loaded from cohort1_auto.pvar.
# Note: No phenotype data present.
# Calculating allele frequencies... done.
# --geno: 10724690 variants removed due to missing genotype data.
# 12194714 variants remaining after main filters.
# Writing cohort2_auto.psam ... done.
# Writing cohort2_auto.pvar ... done.
# Writing cohort2_auto.pgen ... done.

$PLINK2 --threads $THREADS --pfile cohort2_auto \
  --mind $MIND --make-pgen --out cohort3_auto
# 215 samples (0 females, 0 males, 215 ambiguous; 215 founders) loaded from
# cohort2_auto.psam.
# 12194714 variants loaded from cohort2_auto.pvar.
# Note: No phenotype data present.
# Calculating sample missingness rates... done.
# 15 samples removed due to missing genotype data (--mind).
# IDs written to cohort3_auto.mindrem.id .
# 200 samples (0 females, 0 males, 200 ambiguous; 200 founders) remaining after
# main filters.
# Writing cohort3_auto.psam ... done.
# Writing cohort3_auto.pvar ... done.
# Writing cohort3_auto.pgen ... done.

##### Samples removed after MIND=0.1
# cat cohort3_auto.mindrem.id
# #IID
# SRR1018300_sortedbycoords_MD_BQ
# SRR606394_sorted_MD_BQ
# SRR606396_sorted_MD_BQ
# SRR619138_sortedbycoords_MD_BQ
# SRR619144_sortedbycoords_MD_BQ
# SRR619150_sortedbycoords_MD_BQ
# SRR619156_sortedbycoords_MD_BQ
# SRR619160_sortedbycoords_MD_BQ
# SRR619166_sortedbycoords_MD_BQ
# SRR619171_sortedbycoords_MD_BQ
# SRR619178_sortedbycoords_MD_BQ
# SRR619184_sortedbycoords_MD_BQ
# SRR619189_sortedbycoords_MD_BQ
# SRR619190_sortedbycoords_MD_BQ
# SRR619196_sortedbycoords_MD_BQ

$PLINK2 --threads $THREADS --pfile cohort3_auto \
  --maf $MAF_CUTOFF --make-pgen --out cohort4_auto
# 200 samples (0 females, 0 males, 200 ambiguous; 200 founders) loaded from
# cohort3_auto.psam.
# 12194714 variants loaded from cohort3_auto.pvar.
# Note: No phenotype data present.
# Calculating allele frequencies... done.
# 5735237 variants removed due to allele frequency threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 6459477 variants remaining after main filters.
# Writing cohort4_auto.psam ... done.
# Writing cohort4_auto.pvar ... done.
# Writing cohort4_auto.pgen ... done.

### LD pruning
echo ">>> Step 4: LD pruning with r^2 <= $PRUNE_R2"
$PLINK2 --threads $THREADS --pfile cohort4_auto \
  --indep-pairwise ${PRUNE_WINDOW_KB} ${PRUNE_STEP} ${PRUNE_R2} \
  --out tmp/pruned

$PLINK2 --threads $THREADS --pfile cohort4_auto --extract tmp/pruned.prune.in \
  --make-bed --out cohort5_pruned
# 200 samples (0 females, 0 males, 200 ambiguous; 200 founders) loaded from
# cohort4_auto.psam.
# 6459477 variants loaded from cohort4_auto.pvar.
# Note: No phenotype data present.
# --extract: 433106 variants remaining.
# 433106 variants remaining after main filters.
# Writing cohort5_pruned.fam ... done.
# Writing cohort5_pruned.bim ... done.
# Writing cohort5_pruned.bed ... done.

### Remove duplicates with KING
echo ">>> Step 5: Remove duplicates with KING cutoff > $KING_CUTOFF"
$PLINK2 --threads $THREADS --bfile cohort5_pruned --make-king-table --king-cutoff $KING_CUTOFF --out logs/king_cohort_duplicates
##### Excluded samples after KING cutoff 0.354
# cat logs/king_cohort_duplicates.king.cutoff.out.id
# #FID	IID
# 0	chord167-blood_S19_sorted_MD_BQ - matches chord158-blood_S11_sorted_MD_BQ with 50% kinship
# 0	SRR619138_sortedbycoords_MD_BQ - matches SRR619137_sortedbycoords_MD_BQ with 50% kinship
# 0	SRR619166_sortedbycoords_MD_BQ - matches SRR619165_sortedbycoords_MD_BQ with 50% kinship
# 0	SRR619172_sortedbycoords_MD_BQ - matches SRR619171_sortedbycoords_MD_BQ with 50% kinship

$PLINK2 --threads $THREADS --bfile cohort5_pruned --keep logs/king_cohort_duplicates.king.cutoff.in.id --make-pgen --out cohort5_unrel
# 200 samples (0 females, 0 males, 200 ambiguous; 200 founders) loaded from
# cohort5_pruned.fam.
# 433106 variants loaded from cohort5_pruned.bim.
# Note: No phenotype data present.
# --keep: 199 samples remaining.
# 199 samples (0 females, 0 males, 199 ambiguous; 199 founders) remaining after
# main filters.
# Writing cohort5_unrel.psam ... done.
# Writing cohort5_unrel.pvar ... done.
# Writing cohort5_unrel.pgen ... done.


######################## II) Processing of 1KG phase 3 SNPs with vcf SNPs to make PCs
# Reference basename (BED/BIM/FAM or PGEN/PVAR/PSAM)
REF1KG_PLINK_PREFIX="$MSIPROJECT/balay011/references/plink_1kg_files/all_hg38"
# 3202-sample panel: FamilyID SampleID FatherID MotherID Sex Population Superpopulation
REF1KG_PANEL="$MSIPROJECT/balay011/references/plink_1kg_files/20130606_g1k_3202_samples_ped_population.txt"

# Convert PGEN to BED if needed
$PLINK2 --threads $THREADS \
  --pfile $MSIPROJECT/balay011/references/plink_1kg_files/all_hg38 \
  --autosome \
  --snps-only just-acgt \
  --max-alleles 2 \
  --geno $GENO \
  --make-bed --out ref/1kg_auto_1
# 3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) loaded from
# /projects/standard/aventeic/balay011/references/plink_1kg_files/all_hg38.psam.
# Note: 2585 nonstandard chromosome codes present.
# 61599150 out of 75193455 variants loaded from
# /projects/standard/aventeic/balay011/references/plink_1kg_files/all_hg38.pvar.
# 2 categorical phenotypes loaded.
# Calculating allele frequencies... done.
# --geno: 0 variants removed due to missing genotype data.
# 61599150 variants remaining after main filters.
# Writing ref/1kg_auto_1.fam ... done.
# Writing ref/1kg_auto_1.bim ... done.
# Writing ref/1kg_auto_1.bed ... done.

$PLINK2 --threads $THREADS \
  --bfile ref/1kg_auto_1 \
  --mind $MIND \
  --make-bed --out ref/1kg_auto_2
# 3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) loaded from
# ref/1kg_auto_1.fam.

# 61599150 variants loaded from ref/1kg_auto_1.bim.
# Note: No phenotype data present.
# Calculating sample missingness rates... done.
# 0 samples removed due to missing genotype data (--mind).
# 3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) remaining
# after main filters.
# Writing ref/1kg_auto_2.fam ... done.
# Writing ref/1kg_auto_2.bim ... done.
# Writing ref/1kg_auto_2.bed ... done.

$PLINK2 --threads $THREADS \
  --bfile ref/1kg_auto_2 \
  --maf $MAF_CUTOFF \
  --make-bed --out ref/1kg_auto_3
# 3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) loaded from
# ref/1kg_auto_2.fam.
# 61599150 variants loaded from ref/1kg_auto_2.bim.
# Note: No phenotype data present.
# Calculating allele frequencies... done.
# 49012263 variants removed due to allele frequency threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 12586887 variants remaining after main filters.
# Writing ref/1kg_auto_3.fam ... done.
# Writing ref/1kg_auto_3.bim ... done.
# Writing ref/1kg_auto_3.bed ... done.

# ### LD pruning
# $PLINK2 --threads $THREADS --bfile ref/1kg_auto_3 \
#   --set-all-var-ids @:#:\$r:\$a \
#   --indep-pairwise ${PRUNE_WINDOW_KB} ${PRUNE_STEP} ${PRUNE_R2} \
#   --out tmp/pruned

# $PLINK2 --threads $THREADS --bfile ref/1kg_auto_3 --extract tmp/pruned.prune.in \
#   --make-bed --out ref/1kg_auto_maf01

### Remove duplicates with KING
$PLINK2 --threads $THREADS --bfile ref/1kg_auto_3 --make-king-table --king-cutoff $KING_CUTOFF --out logs/king_1kg_duplicates
$PLINK2 --threads $THREADS --bfile ref/1kg_auto_3 --keep logs/king_1kg_duplicates.king.cutoff.in.id --make-pgen --out ref/1kg_auto_unrel
# 3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) loaded from
# ref/1kg_auto_3.fam.
# 12586887 variants loaded from ref/1kg_auto_3.bim.
# Note: No phenotype data present.
# --keep: 3202 samples remaining.
# 3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) remaining
# after main filters.
# Writing ref/1kg_auto_unrel.psam ... done.
# Writing ref/1kg_auto_unrel.pvar ... done.
# Writing ref/1kg_auto_unrel.pgen ... done.

### Re-ID each dataset to CHR:POS:REF:ALT
$PLINK2 --threads $THREADS --pfile ref/1kg_auto_unrel \
  --snps-only just-acgt --max-alleles 2 \
  --set-all-var-ids @:#:\$r:\$a --rm-dup force-first --new-id-max-allele-len 200 \
  --make-pgen --out tmp/1kg_reid
# 3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) loaded from
# ref/1kg_auto_unrel.psam.
# 12586887 variants loaded from ref/1kg_auto_unrel.pvar.
# Note: No phenotype data present.
# Note: Skipping --rm-dup since no duplicate IDs are present.
# 12586887 variants remaining after main filters.
# Writing tmp/1kg_reid.psam ... done.
# Writing tmp/1kg_reid.pvar ... done.
# Writing tmp/1kg_reid.pgen ... done.

$PLINK2 --threads $THREADS --pfile cohort5_unrel \
  --snps-only just-acgt --max-alleles 2 \
  --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 200 \
  --make-pgen --out tmp/cohort_reid
# 199 samples (0 females, 0 males, 199 ambiguous; 199 founders) loaded from
# cohort5_unrel.psam.
# 433106 variants loaded from cohort5_unrel.pvar.
# Note: No phenotype data present.
# 433106 variants remaining after main filters.
# Writing tmp/cohort_reid.psam ... done.
# Writing tmp/cohort_reid.pvar ... done.
# Writing tmp/cohort_reid.pgen ... done.

### Acquire common SNPs between 1KG and cohort
awk 'NR>1{print $3}' tmp/1kg_reid.pvar    | sort -u > tmp/1kg.ids
awk 'NR>1{print $3}' tmp/cohort_reid.pvar | sort -u > tmp/cohort.ids
comm -12 tmp/1kg.ids tmp/cohort.ids > tmp/common.ids
wc -l tmp/1kg.ids tmp/cohort.ids tmp/common.ids
#  12586887 tmp/1kg.ids
#    433106 tmp/cohort.ids
#    240254 tmp/common.ids
#  13260247 total

### Extract common SNPs between 1KG and cohort
$PLINK2 --threads $THREADS --pfile tmp/1kg_reid \
  --extract tmp/common.ids \
  --make-bed --out pca/1kg_aligned
# 3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) loaded from
# tmp/1kg_reid.psam.
# 12586887 variants loaded from tmp/1kg_reid.pvar.
# Note: No phenotype data present.
# --extract: 240254 variants remaining.
# 240254 variants remaining after main filters.
# Writing pca/1kg_aligned.fam ... done.
# Writing pca/1kg_aligned.bim ... done.
# Writing pca/1kg_aligned.bed ... done.

$PLINK2 --threads $THREADS --pfile tmp/cohort_reid \
  --extract tmp/common.ids \
  --make-bed --out pca/cohort_aligned
# 199 samples (0 females, 0 males, 199 ambiguous; 199 founders) loaded from
# tmp/cohort_reid.psam.
# 433106 variants loaded from tmp/cohort_reid.pvar.
# Note: No phenotype data present.
# --extract: 240254 variants remaining.
# 240254 variants remaining after main filters.
# Writing pca/cohort_aligned.fam ... done.
# Writing pca/cohort_aligned.bim ... done.
# Writing pca/cohort_aligned.bed ... done.

diff -s pca/1kg_aligned.bim pca/cohort_aligned.bim 

### Merge cohort with 1KG
printf "%s\n" \
  "${OUTDIR}/pca/1kg_aligned" \
  > "${OUTDIR}/pca/merge_list_bed.txt"

$PLINK --bfile "${OUTDIR}/pca/cohort_aligned" \
  --merge-list "${OUTDIR}/pca/merge_list_bed.txt" \
  --make-bed \
  --memory 120000 \
  --out "${OUTDIR}/pca/merged_all_autosomes_bed"

$PLINK2 --bfile "${OUTDIR}/pca/merged_all_autosomes_bed" \
  --make-pgen \
  --out "${OUTDIR}/pca/merged_all_autosomes"
# 3401 samples (1603 females, 1598 males, 200 ambiguous; 2782 founders) loaded
# from /scratch.global/balay011/gwas_qc_struct/pca/merged_all_autosomes_bed.fam.
# 240254 variants loaded from
# /scratch.global/balay011/gwas_qc_struct/pca/merged_all_autosomes_bed.bim.
# Note: No phenotype data present.
# Writing /scratch.global/balay011/gwas_qc_struct/pca/merged_all_autosomes.psam
# ... done.
# Writing /scratch.global/balay011/gwas_qc_struct/pca/merged_all_autosomes.pvar
# ... done.
# Writing /scratch.global/balay011/gwas_qc_struct/pca/merged_all_autosomes.pgen
# ... done.

### PCA 
$PLINK2 --pfile "${OUTDIR}/pca/merged_all_autosomes" \
  --pca approx 20 allele-wts \
  --out "${OUTDIR}/pca/pca_all"
# 3401 samples (1603 females, 1598 males, 200 ambiguous; 2782 founders) loaded
# from /scratch.global/balay011/gwas_qc_struct/pca/merged_all_autosomes.psam.
# 240254 variants loaded from
# /scratch.global/balay011/gwas_qc_struct/pca/merged_all_autosomes.pvar.
# Note: No phenotype data present.
# Calculating allele frequencies... done.
# Warning: "--pca approx" is only recommended for analysis of >5000 samples.
# Projecting random vectors (8 compute threads)... 
# Projecting random vectors (8 compute threads)... 21/21.
# Computing SVD of Krylov matrix... done.
# Recovering top PCs from range approximation... done.
# --pca approx: Allele weights written to
# /scratch.global/balay011/gwas_qc_struct/pca/pca_all.eigenvec.allele .
# --pca approx: Eigenvectors written to
# /scratch.global/balay011/gwas_qc_struct/pca/pca_all.eigenvec , and eigenvalues
# written to /scratch.global/balay011/gwas_qc_struct/pca/pca_all.eigenval .

echo "PCA done:"
echo "  PCs:     ${OUTDIR}/pca/pca_all.eigenvec"
echo "  Eigenval ${OUTDIR}/pca/pca_all.eigenval"
echo "  Weights: ${OUTDIR}/pca/pca_all.eigenvec.allele_wts"

######################## III) Global ancestry via ADMIXTURE K=5 (unsupervised)
echo ">>> Run unsupervised ADMIXTURE on 1KG panel with k=5"
admixture pca/1kg_aligned.bed 5
echo ">>> Use learned allele frequencies as (fixed) input to next step"
cp 1kg_aligned.5.P cohort_aligned.5.P.in
cp pca/cohort_aligned.bed cohort_aligned.bed
cp pca/cohort_aligned.bim cohort_aligned.bim
cp pca/cohort_aligned.fam cohort_aligned.fam
echo ">>> Run projection ADMIXTURE with k=5"
admixture -P cohort_aligned.bed 5

echo ">>> DONE at $(date)"
