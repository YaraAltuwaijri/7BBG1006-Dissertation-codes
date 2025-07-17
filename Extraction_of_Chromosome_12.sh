#!/bin/bash -l
#SBATCH --job-name=extract_chr12
#SBATCH --output=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/chr12_only/extract_chr12_%j.out
#SBATCH --error=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/chr12_only/extract_chr12_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

# Set input and output directories
INPUT_DIR="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup"
OUTPUT_DIR="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/chr12_only"

# Extract chromosome 12 rows (including header)
for f in "$INPUT_DIR"/*_ready.tsv; do
  base=$(basename "$f" .tsv)
  awk 'NR==1 || $1=="12"' "$f" > "$OUTPUT_DIR/${base}_chr12.tsv"
done

