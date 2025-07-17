#!/bin/bash -l
#SBATCH --job-name=convert_to_vcf
#SBATCH --partition=cpu
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/vcf_ready/convert_%j.out
#SBATCH --error=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/vcf_ready/convert_%j.err

echo "Starting VCF conversion at $(date)"

OUT_DIR="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/vcf_ready"
module load samtools

# Tool paths
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

HEADER_FIX="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/vcf_ready/b.txt"
COLHEADERS="/scratch/prj/mmg_data/software/useful_scripts/combined_colheaders2.tsv"
REF37="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/FASTA_ref_37_38/ref37/hs37d5.fa"
REF38="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/FASTA_ref_37_38/ref38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INPUT_DIR="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/chr12_only_withLP"

[ ! -f "${REF37}.fai" ] && samtools faidx "$REF37"
[ ! -f "${REF38}.fai" ] && samtools faidx "$REF38"

for FILE in "$INPUT_DIR"/*.tsv; do
  BASENAME=$(basename "$FILE" .tsv)
  LABEL=${BASENAME%_ready} 
  if [[ "$BASENAME" == *_buildGRCh37_* ]]; then
    REF="$REF37"
    echo "→ Using GRCh37 for $BASENAME"
  elif [[ "$BASENAME" == *_buildGRCh38_* ]]; then
    REF="$REF38"
    echo "→ Using GRCh38 for $BASENAME"
  else
    echo "Skipping $BASENAME — genome build not specified"
    continue
  fi

  echo "Filtering and converting: $FILE → $OUT_DIR/${LABEL}.vcf.gz"

  awk 'BEGIN {FS=OFS="\t"} NR==1 || ($2 != "NA" && $4 != "NA" && $5 != "NA")' "$FILE" > filtered.tsv

  # Step 1: Convert to basic VCF using bcftools +munge
  bcftools +munge --no-version \
    -f "$REF" \
    -s "$LABEL" \
    -C "$COLHEADERS" \
    -Ov filtered.tsv > temp.vcf

  # Step 2: Inject correct FORMAT metadata
  bcftools annotate --header-lines "$HEADER_FIX" temp.vcf -Ov > temp_fixed.vcf

  # Step 3: Sort and compress
  bcftools sort temp_fixed.vcf -Oz -o "$OUT_DIR/${LABEL}#.vcf.gz"

  # Step 4: Index
  bcftools index -t "$OUT_DIR/${LABEL}#.vcf.gz"

  # Step 5: Cleanup
  rm -f filtered.tsv temp.vcf temp_fixed.vcf
done

echo "VCF conversion completed at $(date)"

