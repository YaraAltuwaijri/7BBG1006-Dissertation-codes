#!/bin/bash -l
#SBATCH -p cpu,simpson_cpu
#SBATCH --output=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/merge_vcf/%j.out
#SBATCH --time=0-4:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=4
#SBATCH --nodes=1

# Tools
bcftools() {
  /scratch/prj/mmg_data/software/bcftools/bcftools-1.20/bcftools "$@"
}

# Paths
VCF_DIR=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/region_subset
OUT_DIR=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/merge_vcf
OUT_PREFIX=ALDH2_merged
OUT_VCF=$OUT_DIR/${OUT_PREFIX}.vcf.gz

mkdir -p "$OUT_DIR"
cd "$VCF_DIR" || exit 1

echo "Merging region-extracted VCFs..."
bcftools merge -m none -Oz -o "$OUT_VCF" *.vcf.gz

echo "Indexing merged VCF..."
bcftools index "$OUT_VCF"

echo "Merged and indexed: $OUT_VCF"
