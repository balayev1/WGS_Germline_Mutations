#!/usr/bin/env bash

#SBATCH --job-name=gwas_prep
#SBATCH --partition=msibigmem,msilarge,msismall
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=balay011@umn.edu

set -euo pipefail
echo "[$(date)] Job started on ${HOSTNAME} with ${SLURM_CPUS_PER_TASK:-unknown} CPUs"

######################## This script takes directories with input vcf SNP and INDEL files and runs GWAS

# Threads from SLURM (fallback 8)
THREADS=16

### Input vcf SNP and INDEL directory
SNPS_DIR="/scratch.global/balay011/germline_calls/vqsr/sample_stats/final_merged/snps"
INDELS_DIR="/scratch.global/balay011/germline_calls/vqsr/sample_stats/final_merged/indels"
ONEKG_DIR="/scratch.global/balay011/germline_calls/1kg_samples/vcf_files"

### Tools/modules
module load bcftools/1.16-gcc-8.2.0-5d4xg4y
module load conda
source activate gwas_env

BCFTOOLS="bcftools"
PLINK2="plink2"
PLINK="plink"

### Output directory
OUTDIR="/scratch.global/balay011/gwas_analysis"
mkdir -p "${OUTDIR}"/{qc,tmp,eur,eas}
cd "$OUTDIR"

### Other files
REF="$MSIPROJECT/balay011/references/reference_genome/GRCh38.primary_assembly.genome.fa"
ALL_EUR_SAMPLE_LIST="$OUTDIR/ctrl_eur_samples.txt"         # one IID per line (controls with european ancestry)
ALL_EAS_SAMPLE_LIST="$OUTDIR/ctrl_eas_samples.txt"         # one IID per line (controls with eastern asian ancestry)
CHORD_EUR_SAMPLE_LIST="$OUTDIR/chord_eur_samples.txt"    # one IID per line (chordoma with > 50% european ancestry)
CHORD_EAS_SAMPLE_LIST="$OUTDIR/chord_eas_samples.txt"    # one IID per line (chordoma with > 50% eastern asian ancestry)
MERGED_FINAL_SNPS_VCF="/scratch.global/balay011/germline_calls/vqsr/sample_stats/final_merged/cohort.snps.merged.vcf.gz"
MERGED_FINAL_INDELS_VCF="/scratch.global/balay011/germline_calls/vqsr/sample_stats/final_merged/cohort.indels.merged.vcf.gz"

######################## Subset 1KG panel by ancestry (EUR/EAS)
echo ">>> Subsetting 1KG per-chrom VCFs by ancestry from ${ONEKG_DIR}"
mkdir -p "${ONEKG_DIR}/eur" "${ONEKG_DIR}/eas"

for ANC in eur eas; do
  if [[ "$ANC" == "eur" ]]; then
    SAMPLE_LIST="$ALL_EUR_SAMPLE_LIST"
  else
    SAMPLE_LIST="$ALL_EAS_SAMPLE_LIST"
  fi

  [[ -s "$SAMPLE_LIST" ]] || { echo ">>> Skipping ${ANC} (empty or missing sample list: $SAMPLE_LIST)"; continue; }
  outdir="${ONEKG_DIR}/${ANC}"
  mkdir -p "${outdir}"

  rm -f "${outdir}/${ANC}.concat.list"
  for CHR in {1..22} X; do
    IN="${ONEKG_DIR}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.recalibrated_variants.vcf.gz"
    if [[ ! -f "$IN" ]]; then
      echo "WARN: Missing ${IN}; skipping chr${CHR}"
      continue
    fi

    echo ">>> ${ANC} chr${CHR}: subset → normalize (left-normalize + split multiallelic)"
    ${BCFTOOLS} view --threads ${THREADS} -S "$SAMPLE_LIST" -Ou "$IN" --force-samples | \
    ${BCFTOOLS} norm  --threads ${THREADS} -f "$REF" -m -both -Oz \
      -o "${outdir}/${ANC}.chr${CHR}.norm.vcf.gz"
    tabix -p vcf "${outdir}/${ANC}.chr${CHR}.norm.vcf.gz"

    echo "${outdir}/${ANC}.chr${CHR}.norm.vcf.gz" >> "${outdir}/${ANC}.concat.list"
  done

  [[ -s "${outdir}/${ANC}.concat.list" ]] || { echo "ERROR: No normalized chr files for ${ANC}"; continue; }

  echo ">>> ${ANC}: concatenating chromosomes"
  ${BCFTOOLS} concat --threads ${THREADS} -f "${outdir}/${ANC}.concat.list" -Oz \
    -o "${outdir}/${ANC}.1KG.norm.vcf.gz"
  tabix -p vcf "${outdir}/${ANC}.1KG.norm.vcf.gz"

  echo ">>> ${ANC}: split SNPs / INDELs"
  ${BCFTOOLS} view --threads ${THREADS} -v snps   -Oz -o "${outdir}/${ANC}.1KG.snps.vcf.gz"   "${outdir}/${ANC}.1KG.norm.vcf.gz"
  tabix -p vcf "${outdir}/${ANC}.1KG.snps.vcf.gz"
  ${BCFTOOLS} view --threads ${THREADS} -v indels -Oz -o "${outdir}/${ANC}.1KG.indels.vcf.gz" "${outdir}/${ANC}.1KG.norm.vcf.gz"
  tabix -p vcf "${outdir}/${ANC}.1KG.indels.vcf.gz"

  # optional cleanup
  rm -f ${outdir}/${ANC}.chr*.norm.vcf.gz* "${outdir}/${ANC}.concat.list"
done

echo "[$(date)] Done. Final 1KG files:"
echo "  EUR: ${ONEKG_DIR}/eur/eur.1KG.snps.vcf.gz   and   ${ONEKG_DIR}/eur/eur.1KG.indels.vcf.gz"
echo "  EAS: ${ONEKG_DIR}/eas/eas.1KG.snps.vcf.gz   and   ${ONEKG_DIR}/eas/eas.1KG.indels.vcf.gz"


######################## Convert 1KG panel VCFs to PLINK PGEN format using sex information
$PLINK2 --vcf ${ONEKG_DIR}/eur/eur.1KG.snps.vcf.gz --psam $OUTDIR/1kg.ctrl.eur.sex.psam --split-par b38 --make-pgen --out $OUTDIR/eur/eur_ctrls_snps_raw
$PLINK2 --vcf ${ONEKG_DIR}/eur/eur.1KG.indels.vcf.gz --psam $OUTDIR/1kg.ctrl.eur.sex.psam --split-par b38 --make-pgen --out $OUTDIR/eur/eur_ctrls_indels_raw
$PLINK2 --vcf ${ONEKG_DIR}/eas/eas.1KG.snps.vcf.gz --psam $OUTDIR/1kg.ctrl.eas.sex.psam  --split-par b38 --make-pgen --out $OUTDIR/eas/eas_ctrls_snps_raw
$PLINK2 --vcf ${ONEKG_DIR}/eas/eas.1KG.indels.vcf.gz --psam $OUTDIR/1kg.ctrl.eas.sex.psam --split-par b38 --make-pgen --out $OUTDIR/eas/eas_ctrls_indels_raw

######################## Infex sex of 1KG panel samples for verification purposes
$PLINK2 --pfile $OUTDIR/eur/eur_ctrls_snps_raw \
  --check-sex max-female-xf=0.2 min-male-xf=0.8 \
  --threads $THREADS --out eur/eur_ctrls_snps_sexcheck
$PLINK2 --pfile $OUTDIR/eas/eas_ctrls_snps_raw \
  --check-sex max-female-xf=0.2 min-male-xf=0.8 \
  --threads $THREADS --out eas/eas_ctrls_snps_sexcheck

# This writes: ${OUT}.sexcheck with columns:
# FID IID PEDSEX SNPSEX STATUS F ...
#   PEDSEX = your provided sex (1/2/0)
#   SNPSEX = inferred (1/2/0=ambiguous)
#   F = X inbreeding coefficient (≈ males → ~1, females → ~0)

######################## Create PLINK .psam file for tumor samples sex information
zgrep -m 1 "^#CHROM" $MERGED_FINAL_SNPS_VCF | cut -f10- | tr '\t' '\n' \
 | awk 'BEGIN{print "#FID IID SEX"} {print $1, $1, 0}' OFS='\t' \
 > tumors.sex.psam

######################## Convert merged tumor VCFs to PLINK2 PGEN
echo ">>> Import SNPs directly to PGEN"
$PLINK2 --threads $THREADS --psam $OUTDIR/tumors.sex.psam --vcf $MERGED_FINAL_SNPS_VCF --split-par b38 --make-pgen --out cohort1_cases_snps

echo ">>> Import INDELs directly to PGEN"
$PLINK2 --threads $THREADS --psam $OUTDIR/tumors.sex.psam --vcf $MERGED_FINAL_INDELS_VCF --split-par b38 --make-pgen --out cohort1_cases_indels

######################## Sex inference by chrX chromosome
echo ">>> Inferring sex from chrX (split PAR)"
$PLINK2 --pfile cohort1_cases_snps --check-sex max-female-xf=0.2 min-male-xf=0.8 --threads ${THREADS} --out qc/tumors.sex
# Comments: despite extremely low F values ~ -1 for females, most of them are actual females except Chord_55b which may be outlier

######################## SNP call rate stats
echo ">>> Estimate per-sample SNP call rates"
$PLINK2 --pfile cohort1_cases_snps \
       --missing \
       --out sample_snp_callrate

######################## INDEL call rate stats
echo ">>> Estimate per-sample INDEL call rates"
$PLINK2 --pfile cohort1_cases_indels \
       --missing \
       --out sample_indel_callrate

######################## Split SNPs and INDEL VCFs by ancestry 
echo ">>> Splitting cases by ancestry (EUR/EAS)"
${BCFTOOLS} view --threads "$THREADS" \
  -S "$CHORD_EUR_SAMPLE_LIST" -Ou $MERGED_FINAL_SNPS_VCF | \
${BCFTOOLS} norm --threads "$THREADS" \
  -f "$REF" -m -both -Ou | \
${BCFTOOLS} sort -Oz -o cohort.eur.cases.snps.vcf.gz
tabix -p vcf cohort.eur.cases.snps.vcf.gz

${BCFTOOLS} view --threads "$THREADS" \
  -S "$CHORD_EAS_SAMPLE_LIST" -Ou $MERGED_FINAL_SNPS_VCF | \
${BCFTOOLS} norm --threads "$THREADS" \
  -f "$REF" -m -both -Ou | \
${BCFTOOLS} sort -Oz -o cohort.eas.cases.snps.vcf.gz
tabix -p vcf cohort.eas.cases.snps.vcf.gz

${BCFTOOLS} view --threads "$THREADS" \
  -S "$CHORD_EUR_SAMPLE_LIST" -Ou $MERGED_FINAL_INDELS_VCF | \
${BCFTOOLS} norm --threads "$THREADS" \
  -f "$REF" -m -both -Ou | \
${BCFTOOLS} sort -Oz -o cohort.eur.cases.indels.vcf.gz
tabix -p vcf cohort.eur.cases.indels.vcf.gz

${BCFTOOLS} view --threads "$THREADS" \
  -S "$CHORD_EAS_SAMPLE_LIST" -Ou $MERGED_FINAL_INDELS_VCF | \
${BCFTOOLS} norm --threads "$THREADS" \
  -f "$REF" -m -both -Ou | \
${BCFTOOLS} sort -Oz -o cohort.eas.cases.indels.vcf.gz
tabix -p vcf cohort.eas.cases.indels.vcf.gz

######################## Split sex inference files by ancestry
grep -F -f "$CHORD_EUR_SAMPLE_LIST" qc/tumors.sex.sexcheck | cut -f1,2,4 > qc/tumors.eur.sex.psam
grep -F -f "$CHORD_EAS_SAMPLE_LIST" qc/tumors.sex.sexcheck | cut -f1,2,4 > qc/tumors.eas.sex.psam
grep -F -f "$ALL_EUR_SAMPLE_LIST" eur/eur_ctrls_snps_sexcheck.sexcheck | cut -f1,2,4 > qc/1kg.ctrl.eur.sex.psam
grep -F -f "$ALL_EAS_SAMPLE_LIST" eas/eas_ctrls_snps_sexcheck.sexcheck | cut -f1,2,4 > qc/1kg.ctrl.eas.sex.psam

######################## Generate phenotype files for GWAS analysis
echo ">>> Generating phenotype files for GWAS analysis"
echo -e "FID\tIID\tTUMOR" > eur/pheno_eur.txt
awk '{print $1, $2, 2}' OFS='\t' qc/tumors.eur.sex.psam >> eur/pheno_eur.txt
awk '{print $1, $2, 1}' OFS='\t' qc/1kg.ctrl.eur.sex.psam >> eur/pheno_eur.txt

echo -e "FID\tIID\tTUMOR" > eas/pheno_eas.txt
awk '{print $1, $2, 2}' OFS='\t' qc/tumors.eas.sex.psam >> eas/pheno_eas.txt
awk '{print $1, $2, 1}' OFS='\t' qc/1kg.ctrl.eas.sex.psam >> eas/pheno_eas.txt

echo "[$(date)] Done."
