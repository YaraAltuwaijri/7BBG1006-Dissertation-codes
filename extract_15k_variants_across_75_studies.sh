#!/bin/bash
#SBATCH --job-name=15var75stud
#SBATCH --output=/scratch/prj/mmg_data/benchmarking/15kVar75stud/64G_15var75_%j.out
#SBATCH --error=/scratch/prj/mmg_data/benchmarking/15kVar75stud/64G_15var75_%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

module load tabix
module load bcftools

INPUT_DIR="/scratch/prj/mmg_data/benchmarking"
OUTPUT_DIR="$INPUT_DIR/15kVar75stud"

SIZES=(100 500 1000 5000 10k)

for SIZE in "${SIZES[@]}"; do
  echo "=== STARTING SIZE: $SIZE ==="
  PLAIN_TSV="$INPUT_DIR/synthetic_${SIZE}.tsv"
  INDEXED_TSV="$INPUT_DIR/synthetic_${SIZE}.tsv.gz"
  PLAIN_VCF="$INPUT_DIR/synthetic_${SIZE}.vcf"
  INDEXED_VCF="$INPUT_DIR/synthetic_${SIZE}.vcf.gz"
  TMPDIR="$OUTPUT_DIR/${SIZE}_tmp"
  LOGFILE="$OUTPUT_DIR/64query_timing_${SIZE}.log"

  mkdir -p "$TMPDIR" "$OUTPUT_DIR"

  ###############################
  # Format A: Plain TSV
  ###############################
echo "[Timing] Plain TSV full query" >> "$LOGFILE"
/usr/bin/time -v bash -c "
  head -1 '$PLAIN_TSV' | sed 's/\t/\n/g' | grep -E '\\.ES\$' | sed 's/\\.ES\$//' | shuf -n 75 | sort > '$TMPDIR/studies.list'
  head -1 '$PLAIN_TSV' | sed 's/\t/\n/g' | cat -n > '$TMPDIR/header.rownums'
  awk '{print \$1\".ES\\n\"\$1\".SE\\n\"\$1\".LP\\n\"\$1\".AF\"}' '$TMPDIR/studies.list' > '$TMPDIR/study_fields.list'
  cut -f1 -d\$'\t' <(grep -Ff '$TMPDIR/study_fields.list' '$TMPDIR/header.rownums') | sed '1s/^/3,/' | paste -sd, - > '$TMPDIR/cutcols.string'
  cut -f3 '$PLAIN_TSV' | sed 1d | shuf -n 15000 > '$TMPDIR/varIDs.list'
  awk 'NR==FNR {ids[\$1]; next} \$3 in ids' '$TMPDIR/varIDs.list' '$PLAIN_TSV' > '$TMPDIR/matched.tsv'
  cut -f\$(cat '$TMPDIR/cutcols.string') '$TMPDIR/matched.tsv' > '$OUTPUT_DIR/synthetic_${SIZE}_plain_tsv_15k_75stud.tsv'
  head -1 '$PLAIN_TSV' | cut -f\$(cat '$TMPDIR/cutcols.string') > '$TMPDIR/header_line.tsv'
  cat '$TMPDIR/header_line.tsv' '$OUTPUT_DIR/synthetic_${SIZE}_plain_tsv_15k_75stud.tsv' > '$OUTPUT_DIR/synthetic_${SIZE}_plain_tsv_15k_75stud_with_header.tsv'
  mv '$OUTPUT_DIR/synthetic_${SIZE}_plain_tsv_15k_75stud_with_header.tsv' '$OUTPUT_DIR/synthetic_${SIZE}_plain_tsv_15k_75stud.tsv'
" 2>> "$LOGFILE"

  ###############################
  # Format B: Indexed TSV
  ###############################
  echo "[Timing] Indexed TSV full query" >> "$LOGFILE"
  /usr/bin/time -v bash -c "
    awk -F'[:_]' '{print \$1\":\"\$2\"-\"\$2}' '$TMPDIR/varIDs.list' > '$TMPDIR/tabix_regions.list'
    > '$TMPDIR/matched.tsv'
    while read -r region; do tabix '$INDEXED_TSV' \$region >> '$TMPDIR/matched.tsv'; done < '$TMPDIR/tabix_regions.list'
    cut -f\$(cat '$TMPDIR/cutcols.string') '$TMPDIR/matched.tsv' > '$OUTPUT_DIR/synthetic_${SIZE}_indexed_tsv_15k_75stud.tsv'
    zcat '$INDEXED_TSV' | head -1 | cut -f\$(cat '$TMPDIR/cutcols.string') > '$TMPDIR/header_line.tsv'
    cat '$TMPDIR/header_line.tsv' '$OUTPUT_DIR/synthetic_${SIZE}_indexed_tsv_15k_75stud.tsv' > '$OUTPUT_DIR/synthetic_${SIZE}_indexed_tsv_15k_75stud_with_header.tsv'
    mv '$OUTPUT_DIR/synthetic_${SIZE}_indexed_tsv_15k_75stud_with_header.tsv' '$OUTPUT_DIR/synthetic_${SIZE}_indexed_tsv_15k_75stud.tsv'
  " 2>> "$LOGFILE"

###############################
# Format C: Plain VCF Full Query
###############################
echo "[Timing] Plain VCF full query" >> "$LOGFILE"

/usr/bin/time -v bash -c "
  grep -v '^##' '$PLAIN_VCF' | head -1 > '$TMPDIR/header_plain_vcf.txt'
  cut -f10- '$TMPDIR/header_plain_vcf.txt' | tr '\t' '\n' | shuf -n 75 | sort > '$TMPDIR/studies_vcf.list'
  cut -f3 '$PLAIN_TSV' | sed 1d | shuf -n 15000 > '$TMPDIR/varIDs.list'

  bcftools view -H -i 'ID=@$TMPDIR/varIDs.list' '$PLAIN_VCF' > '$TMPDIR/matched_plain_vcf.txt'
  cp '$TMPDIR/matched_plain_vcf.txt' '$TMPDIR/filtered_plain_vcf.vcf'

  echo -e \"ID\t\$(awk '{printf \"%s.ES\t%s.SE\t%s.LP\t%s.AF\t\", \$1, \$1, \$1, \$1}' '$TMPDIR/studies_vcf.list' | sed 's/\t\$//')\" > '$OUTPUT_DIR/synthetic_${SIZE}_plain_vcf_15k_75stud.tsv'

  awk -v OFS='\t' -v header_file='$TMPDIR/header_plain_vcf.txt' -v studies='$TMPDIR/studies_vcf.list' '
  BEGIN {
    while ((getline s < studies) > 0) selected[s] = 1
    while ((getline h < header_file) > 0) {
      split(h, head, \"\t\")
      for (i = 10; i <= length(head); i++) if (head[i] in selected) selected_cols[i] = head[i]
    }
  }
  {
    id = \$3
    out = id
    split(\$9, fmt, \":\")
    for (i = 10; i <= NF; i++) {
      if (!(i in selected_cols)) continue
      split(\$i, val, \":\")
      es=se=lp=af=\"\"
      for (j in fmt) {
        if (fmt[j] == \"ES\") es = val[j]
        if (fmt[j] == \"SE\") se = val[j]
        if (fmt[j] == \"LP\") lp = (val[j] != \"\" && val[j] != \".\") ? 10^(-val[j]) : \"\"
        if (fmt[j] == \"AF\") af = val[j]
      }
      out = out OFS es OFS se OFS lp OFS af
    }
    print out
  }' '$TMPDIR/filtered_plain_vcf.vcf' >> '$OUTPUT_DIR/synthetic_${SIZE}_plain_vcf_15k_75stud.tsv'
" 2>> "$LOGFILE"

###############################
# Format D: Indexed VCF Full Query
###############################
echo "[Timing] Indexed VCF full query" >> "$LOGFILE"

/usr/bin/time -v bash -c "
  zcat '$INDEXED_VCF' | grep -v '^##' | head -1 > '$TMPDIR/header_indexed_vcf.txt'
  cut -f10- '$TMPDIR/header_indexed_vcf.txt' | tr '\t' '\n' | shuf -n 75 | sort > '$TMPDIR/studies_vcf.list'
  cut -f3 '$PLAIN_TSV' | sed 1d | shuf -n 15000 > '$TMPDIR/varIDs.list'

  bcftools view -H -i 'ID=@$TMPDIR/varIDs.list' '$INDEXED_VCF' > '$TMPDIR/matched_indexed_vcf.txt'
  cp '$TMPDIR/matched_indexed_vcf.txt' '$TMPDIR/filtered_indexed_vcf.vcf'

  echo -e \"ID\t\$(awk '{printf \"%s.ES\t%s.SE\t%s.LP\t%s.AF\t\", \$1, \$1, \$1, \$1}' '$TMPDIR/studies_vcf.list' | sed 's/\t\$//')\" > '$OUTPUT_DIR/synthetic_${SIZE}_indexed_vcf_15k_75stud.tsv'

  awk -v OFS='\t' -v header_file='$TMPDIR/header_indexed_vcf.txt' -v studies='$TMPDIR/studies_vcf.list' '
  BEGIN {
    while ((getline s < studies) > 0) selected[s] = 1
    while ((getline h < header_file) > 0) {
      split(h, head, \"\t\")
      for (i = 10; i <= length(head); i++) if (head[i] in selected) selected_cols[i] = head[i]
    }
  }
  {
    id = \$3
    out = id
    split(\$9, fmt, \":\")
    for (i = 10; i <= NF; i++) {
      if (!(i in selected_cols)) continue
      split(\$i, val, \":\")
      es=se=lp=af=\"\"
      for (j in fmt) {
        if (fmt[j] == \"ES\") es = val[j]
        if (fmt[j] == \"SE\") se = val[j]
        if (fmt[j] == \"LP\") lp = (val[j] != \"\" && val[j] != \".\") ? 10^(-val[j]) : \"\"
        if (fmt[j] == \"AF\") af = val[j]
      }
      out = out OFS es OFS se OFS lp OFS af
    }
    print out
  }' '$TMPDIR/filtered_indexed_vcf.vcf' >> '$OUTPUT_DIR/synthetic_${SIZE}_indexed_vcf_15k_75stud.tsv'
" 2>> "$LOGFILE"

done
