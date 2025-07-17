#!/bin/bash

# ALDH2 ± 500,000 bp
REGION="12:111266887-112317532"

INDIR=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/liftover
OUTDIR=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/region_subset

for f in "$INDIR"/*_buildGRCh38_ready_chr12_withLP#.vcf.gz ; do
  fname=$(basename "$f")
  out="$OUTDIR/${fname/.vcf.gz/_region.vcf.gz}"
  echo "Extracting $REGION from $fname → $out"

  # Extract and preserve FORMAT fields 
  bcftools view -r "$REGION" -Oz -o "$out" "$f"
  bcftools index -t "$out"
done
