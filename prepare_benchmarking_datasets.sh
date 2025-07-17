#!/bin/bash
#SBATCH --job-name=prepare_benchmarking
#SBATCH --output=prepare_benchmarking_%j.out
#SBATCH --error=prepare_benchmarking_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

# Load necessary modules
module load bcftools/1.12-gcc-13.2.0-python-3.11.6
module load tabix

# Input master VCF
MASTER_VCF="original_synthetic_10k.vcf"

# Step 1. Master Dataset Preparation
# Standardize variant IDs and simplify the header
echo "Preparing master VCF dataset..."
awk -F"\t" 'NR==1 {print "#CHROM\tPOS\tID\tREF\tALT"; next}
NR>1 {
    id = $1 ":" $2 "_" $4 "_" $5;
    print $1"\t"$2"\t"id"\t"$4"\t"$5
}' "$MASTER_VCF" > master_synthetic_10k.tsv

# Step 2. VCF Subsetting
echo "Subsetting master VCF into smaller datasets..."
SAMPLES=($(grep "^#CHROM" "$MASTER_VCF" | cut -f10- | tr '\t' '\n'))
for subset_size in 500 1000 5000; do
    subset_samples=$(printf ",%s" "${SAMPLES[@]:0:$subset_size}" | cut -c2-)
    echo "Creating subset: ${subset_size} samples"
    bcftools view -s "$subset_samples" "$MASTER_VCF" -Ov -o "synthetic_${subset_size}.vcf"
    bgzip -c "synthetic_${subset_size}.vcf" > "synthetic_${subset_size}.vcf.gz"
    tabix -p vcf "synthetic_${subset_size}.vcf.gz"
done

# Step 3. VCF to TSV Conversion
echo "Converting VCF subsets to TSV..."
for subset_size in 500 1000 5000; do
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "synthetic_${subset_size}.vcf" > "synthetic_${subset_size}.tsv"
done

# Step 4. Compression and Indexing of TSVs
echo "Compressing and indexing TSV files..."
for subset_size in 500 1000 5000; do
    gzip -c "synthetic_${subset_size}.tsv" > "synthetic_${subset_size}.tsv.gz"
    tabix -s 1 -b 2 -e 2 "synthetic_${subset_size}.tsv.gz"
done

echo "Benchmarking datasets prepared successfully."
