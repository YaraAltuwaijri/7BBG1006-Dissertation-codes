#!/bin/bash -l
#SBATCH --job-name=liftover_all_chr12
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/liftover/%j.out
#SBATCH --error=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/liftover/%j.err

# === Tools ===
export BCFTOOLS_PLUGINS=/scratch/prj/mmg_data/software/score
bcftools() {
  /scratch/prj/mmg_data/software/bcftools/bcftools-1.20/bcftools "$@"
}
bgzip() {
  /scratch/prj/mmg_data/software/bcftools/bcftools-1.20/htslib-1.20/htslib/bin/bgzip "$@"
}
tabix() {
  /scratch/prj/mmg_data/software/bcftools/bcftools-1.20/htslib-1.20/htslib/bin/tabix "$@"
}

# === References ===
REF37="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/FASTA_ref_37_38/ref37/hs37d5.fa"
REF38="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/FASTA_ref_37_38/ref38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
CHAIN="/scratch/prj/mmg_data/reference/GRCh37_to_GRCh38.chain.gz"

# === Directories ===
IN_DIR="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/vcf_ready"
OUT_DIR="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/liftover"

cd "$IN_DIR" || exit 1

FILES=(*_buildGRCh37_*.vcf.gz)

for FILE in "${FILES[@]}"; do
  BASE=$(basename "$FILE" .vcf.gz)
  echo "▶ Processing $BASE"

  # Step 1: Remove REF mismatches
  bcftools view -e 'FILTER="REF_MISMATCH"' "$FILE" -Oz -o tmp0.vcf.gz
  bcftools index -t tmp0.vcf.gz

  # Step 2: Normalize multiallelics (split)
  bcftools norm --no-version -Oz -m+ tmp0.vcf.gz -f "$REF37" -o tmp1.vcf.gz
  bcftools index -t tmp1.vcf.gz

  # Step 3: Liftover
  bcftools +liftover --no-version -Ou tmp1.vcf.gz -- \
    -f "$REF38" \
    -s "$REF37" \
    -c "$CHAIN" \
    --reject "${OUT_DIR}/${BASE}_rejected_liftover.vcf.gz" \
    --reject-type z \
    --write-src | \
  bcftools sort -Oz -o tmp2.vcf.gz
  bcftools index -t tmp2.vcf.gz

  # Step 4: Final normalization (merge back)
  bcftools norm --no-version -Oz -m- tmp2.vcf.gz -f "$REF38" \
    -o "${OUT_DIR}/${BASE/_buildGRCh37_/_buildGRCh38_}.vcf.gz"
  bcftools index -t "${OUT_DIR}/${BASE/_buildGRCh37_/_buildGRCh38_}.vcf.gz"

  echo "✅ Done: ${BASE/_buildGRCh37_/_buildGRCh38_}.vcf.gz"
done

# Clean up temp files
rm -f tmp0.vcf.gz* tmp1.vcf.gz* tmp2.vcf.gz*

echo "Liftover complete for ${#FILES[@]} files."

