#!/bin/bash
#SBATCH --job-name=region_query
#SBATCH --partition=simpson_cpu,cpu
#SBATCH --output=/scratch/prj/mmg_data/benchmarking/RegionAllstudies/64/64_region_query_%j.out
#SBATCH --error=/scratch/prj/mmg_data/benchmarking/RegionAllstudies/64/64_region_query_%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

module load tabix
module load bcftools

cd /scratch/prj/mmg_data/benchmarking || exit 1

OUTPUT_DIR="/scratch/prj/mmg_data/benchmarking/RegionAllstudies/64"
TIME_LOG="${OUTPUT_DIR}/64_region_query.log"
BASE_NAMES=("synthetic_100" "synthetic_500" "synthetic_1000" "synthetic_5000" "synthetic_10k")

# ========= Warm-up =========
for BASE_NAME in "${BASE_NAMES[@]}"; do
  echo "Warming up files for ${BASE_NAME}..." >> "$TIME_LOG"
  [[ -f "${BASE_NAME}.tsv" ]] && head -n 1000 "${BASE_NAME}.tsv" > /dev/null
  [[ -f "${BASE_NAME}.tsv.gz.tbi" ]] && tabix "${BASE_NAME}.tsv.gz" "22:20000000-20000010" > /dev/null 2>&1
  [[ -f "${BASE_NAME}.vcf" ]] && head -n 1000 "${BASE_NAME}.vcf" > /dev/null
  [[ -f "${BASE_NAME}.vcf.gz.tbi" ]] && tabix "${BASE_NAME}.vcf.gz" "22:20000000-20000010" > /dev/null 2>&1
done

# ========= Header Extract =========
extract_header() {
  local FILE="$1"
  if [[ "$FILE" == *.vcf.gz || "$FILE" == *.vcf ]]; then
    bcftools view -h "$FILE" | grep -v '^##' | head -n 1
  elif [[ "$FILE" == *.tsv.gz ]]; then
    zcat "$FILE" | head -n 1
  else
    head -n 1 "$FILE"
  fi
}

# ========= TSV Parsing =========
parse_tsv() {
  local TMP_OUT="$1"
  local OUT_FILE="$2"

  echo -e "ID\tStudy\tbeta\tse\tpval\tea_freq" > "$OUT_FILE"

  awk '
    BEGIN { FS = OFS = "\t"; }
    NR==1 {
      n = NF;
      for (i = 1; i <= n; i++) {
        colnames[i] = $i;
        if (colnames[i] ~ /\.ES$/) {
          study = substr(colnames[i], 1, length(colnames[i]) - 3);
          es[i] = study;
        } else if (colnames[i] ~ /\.SE$/) se_idx[substr(colnames[i], 1, length(colnames[i]) - 3)] = i;
        else if (colnames[i] ~ /\.LP$/) lp_idx[substr(colnames[i], 1, length(colnames[i]) - 3)] = i;
        else if (colnames[i] ~ /\.AF$/) af_idx[substr(colnames[i], 1, length(colnames[i]) - 3)] = i;
      }
      next;
    }
    {
      for (i in es) {
        study = es[i];
        beta = $i;
        se   = (study in se_idx) ? $(se_idx[study]) : "";
        lp   = (study in lp_idx) ? $(lp_idx[study]) : "";
        af   = (study in af_idx) ? $(af_idx[study]) : "";
        pval = (lp != "" && lp != ".") ? 10^(-lp) : "";
        print $3, study, beta, se, pval, af;
      }
    }
  ' "$TMP_OUT" >> "$OUT_FILE"
}

# ========= VCF Parsing =========
parse_vcf() {
  local TMP_OUT="$1" OUT_FILE="$2" HEADER="$3"
  IFS=$'\t' read -r -a header_cols <<< "$HEADER"
  echo -e "ID\tStudy\tbeta\tse\tpval\tea_freq" > "$OUT_FILE"
  tail -n +1 "$TMP_OUT" | while IFS= read -r row; do
    IFS=$'\t' read -r -a row_vals <<< "$row"
    FORMAT_FIELD="${row_vals[8]}"
    IFS=':' read -r -a format_keys <<< "$FORMAT_FIELD"
    for i in $(seq 9 $((${#row_vals[@]} - 1))); do
      study="${header_cols[$i]}"
      IFS=':' read -r -a values <<< "${row_vals[$i]}"
      for j in "${!format_keys[@]}"; do
        key="${format_keys[$j]}"
        [[ "$key" == "ES" ]] && beta="${values[$j]}"
        [[ "$key" == "SE" ]] && se="${values[$j]}"
        [[ "$key" == "LP" ]] && pval="${values[$j]}"
        [[ "$key" == "AF" ]] && af="${values[$j]}"
      done
      echo -e "${row_vals[2]}\t${study}\t${beta}\t${se}\t${pval}\t${af}" >> "$OUT_FILE"
    done
  done
}

# ========= Main Query Loop =========
for BASE_NAME in "${BASE_NAMES[@]}"; do
  for label in start middle end; do
    case "$label" in
      start) CHR=22; START=17072085; END=17572085 ;;
      middle) CHR=22; START=34611072; END=34886072 ;;
      end) CHR=22; START=50696593; END=51196593 ;;
    esac

    REGION="${CHR}:${START}-${END}"

    ### --- 1. Plain TSV ---
    TMP=$(mktemp)
    OUT="${OUTPUT_DIR}/${label}_${BASE_NAME}_plain_tsv.tsv"
    echo -e "\n[$(date)] File: ${BASE_NAME}.tsv | Region: ${REGION}" >> "$TIME_LOG"
    {
      /usr/bin/time -v awk -F'\t' -v chr="$CHR" -v start="$START" -v end="$END" \
        'NR==1 || ($1 == chr && $2 >= start && $2 <= end)' "${BASE_NAME}.tsv" > "$TMP"

      parse_tsv "$TMP" "$OUT"
    } 2>> "$TIME_LOG"
    rm -f "$TMP"

    ### --- 2. Indexed TSV (.tsv.gz + .tbi) ---
    HEADER=$(extract_header "${BASE_NAME}.tsv.gz")
    TMP=$(mktemp)
    OUT="${OUTPUT_DIR}/${label}_${BASE_NAME}_indexed_tsv.tsv"
    echo -e "\n[$(date)] File: ${BASE_NAME}.tsv.gz | Region: ${REGION}" >> "$TIME_LOG"
    {
      /usr/bin/time -v bash -c "
        { zcat '${BASE_NAME}.tsv.gz' | head -n 1;
          tabix '${BASE_NAME}.tsv.gz' '${REGION}';
        } > '${TMP}'
      "
      parse_tsv "$TMP" "$OUT"
    } 2>> "$TIME_LOG"
    rm -f "$TMP"

    ### --- 3. Plain VCF ---
HEADER=$(extract_header "${BASE_NAME}.vcf")
OUT="${OUTPUT_DIR}/${label}_${BASE_NAME}_plain_vcf.tsv"
echo -e "\n[$(date)] File: ${BASE_NAME}.vcf | Region: ${REGION}" >> "$TIME_LOG"

{
  /usr/bin/time -v bcftools view -r "${REGION}" -Ov "${BASE_NAME}.vcf" > "${OUT}.tmp"
  parse_vcf "${OUT}.tmp" "$OUT" "$HEADER"
} 2>> "$TIME_LOG"

rm -f "${OUT}.tmp"

    ### --- 4. Indexed VCF (.vcf.gz + .tbi) ---
    HEADER=$(extract_header "${BASE_NAME}.vcf.gz")
    TMP=$(mktemp)
    OUT="${OUTPUT_DIR}/${label}_${BASE_NAME}_indexed_vcf.tsv"
    echo -e "\n[$(date)] File: ${BASE_NAME}.vcf.gz | Region: ${REGION}" >> "$TIME_LOG"
    {
      /usr/bin/time -v bcftools view -H -r "$REGION" "${BASE_NAME}.vcf.gz" > "$TMP"
      parse_vcf "$TMP" "$OUT" "$HEADER"
    } 2>> "$TIME_LOG"
    rm -f "$TMP"

  done
done
