#!/bin/bash -l
#SBATCH -p cpu,simpson_cpu
#SBATCH --time=0-4:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --output=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/merge_vcf/%j.out

# Load tools
export BCFTOOLS_PLUGINS=/scratch/prj/mmg_data/software/score
bcftools() {
  /scratch/prj/mmg_data/software/bcftools/bcftools-1.20/bcftools "$@"
}

# Directories and file names
IN_DIR=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/merge_vcf
OUT_DIR=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/merge_vcf

MERGED_VCF="$IN_DIR/ALDH2_merged.vcf.gz"
GWS_VCF="$OUT_DIR/ALDH2_gws.vcf.gz"
TOP_HITS="$OUT_DIR/ALDH2_maxLogP.txt"

echo "Extracting genome-wide significant variants (P < 5e-8)..."
bcftools view -i 'FMT/LP>7.301' "$MERGED_VCF" -Oz -o "$GWS_VCF"
bcftools index -t "$GWS_VCF"

echo "Summarizing top hit per study..."
bcftools annotate -x INFO,^FORMAT/LP -Ov "$MERGED_VCF" |\
grep -v ^## |\
awk '{
  if (NR==1) {
    for (i=10; i<=NF; i++) {
      split($i, a, "_")
      study[i] = a[1]; cpra[i]="none"; maxlp[i]=0
    }
    nc=NF
  } else {
    for (i=10; i<=nc; i++) {
      if ($i != "." && $i > maxlp[i]) {
        cpra[i]=$1":"$2"_"$4"_"$5; maxlp[i]=$i
      }
    }
  }
}
END {
  for (i=10; i<=nc; i++) {
    print study[i], cpra[i], maxlp[i]
  }
}' > "$TOP_HITS"

echo "Done. Outputs:"
echo "  - Genome-wide hits: $GWS_VCF"
echo "  - Top hit per study: $TOP_HITS"
