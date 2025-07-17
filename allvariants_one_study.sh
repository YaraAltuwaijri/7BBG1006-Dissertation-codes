#!/bin/bash 
#SBATCH --job-name=AllvariantOneStudy.job
#SBATCH --output=/scratch/prj/mmg_data/benchmarking/AllvariantOneStudy/AllvariantOneStudy_%j.out
#SBATCH --error=/scratch/prj/mmg_data/benchmarking/AllvariantOneStudy/AllvariantOneStudy_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

module load tabix
module load bcftools

OUTPUT_DIR="/scratch/prj/mmg_data/benchmarking/AllvariantOneStudy/32"

DATASETS=("synthetic_100" "synthetic_500" "synthetic_1000" "synthetic_5000" "synthetic_10k")
N_STUDIES=(100 500 1000 5000 10000)
N_NAMES=("100" "500" "1000" "5000" "10k")

for idx in "${!DATASETS[@]}"; do
    DATASET=${DATASETS[$idx]}
    N=${N_STUDIES[$idx]}
    N_NAME=${N_NAMES[$idx]}
    HALF=$((N / 2))
    STUDY="study_synth${HALF}"
    SUBSET_SIZE=$N_NAME

    echo "Querying ${DATASET} with ${STUDY}"

    # ----- Plain .tsv -----
    echo "1) Plain .tsv query"
    /usr/bin/time -v bash -c "
    awk -F'\t' -v OFS='\t' -v COLSTART=$((5 + HALF * 4 - 3)) '
    NR==1 {print \"#CHROM\", \"POS\", \"ID\", \"REF\", \"ALT\", \$COLSTART, \$(COLSTART+1), \$(COLSTART+2), \$(COLSTART+3); next}
    {print \$1, \$2, \$3, \$4, \$5, \$(COLSTART), \$(COLSTART+1), \$(COLSTART+2), \$(COLSTART+3)}
    ' \"${DATASET}.tsv\" > \"${OUTPUT_DIR}/${SUBSET_SIZE}_plain_tsv.tsv\""

    # ----- Indexed .tsv.gz -----
    echo "2) Indexed .tsv.gz query"
    /usr/bin/time -v bash -c "
    (echo -e '#CHROM\tPOS\tID\tREF\tALT\t${STUDY}.ES\t${STUDY}.SE\t${STUDY}.LP\t${STUDY}.AF' && \
    tabix ${DATASET}.tsv.gz 22:17000000-52000000 | \
    cut -f1-5,$((5 + HALF * 4 - 3))-$((5 + HALF * 4))) \
    > \"${OUTPUT_DIR}/${SUBSET_SIZE}_indexed_tsv.tsv\""

    # ----- Plain VCF -----
    echo "3) Plain .vcf query"
    /usr/bin/time -v bash -c "
    echo -e '#CHROM\tPOS\tID\tREF\tALT\t${STUDY}.ES\t${STUDY}.SE\t${STUDY}.LP\t${STUDY}.AF' > \"${OUTPUT_DIR}/${SUBSET_SIZE}_plain_vcf.tsv\" && \
    bcftools query -s ${STUDY} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%ES\t%SE\t%LP\t%AF\\n]' ${DATASET}.vcf >> \"${OUTPUT_DIR}/${SUBSET_SIZE}_plain_vcf.tsv\""

    # ----- Indexed VCF.gz -----
    echo "4) Indexed .vcf.gz query"
    /usr/bin/time -v bash -c "
    echo -e '#CHROM\tPOS\tID\tREF\tALT\t${STUDY}.ES\t${STUDY}.SE\t${STUDY}.LP\t${STUDY}.AF' > \"${OUTPUT_DIR}/${SUBSET_SIZE}_indexed_vcf.tsv\" && \
    bcftools query -r 22:17000000-52000000 -s ${STUDY} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%ES\t%SE\t%LP\t%AF\\n]' ${DATASET}.vcf.gz >> \"${OUTPUT_DIR}/${SUBSET_SIZE}_indexed_vcf.tsv\""

    echo "---- Finished ${DATASET} ----"
done
