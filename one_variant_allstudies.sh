#!/bin/bash
#SBATCH --job-name=onevariant
#SBATCH --output=/scratch/prj/mmg_data/benchmarking/OneVariantAllStudies/64_onevariant_%j.out
#SBATCH --error=/scratch/prj/mmg_data/benchmarking/OneVariantAllStudies/64_onevariant_%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

module load tabix
module load bcftools

cd /scratch/prj/mmg_data/benchmarking || exit 1
OUTPUT_DIR="/scratch/prj/mmg_data/benchmarking/OneVariantAllStudies/64"
TIME_LOG="${OUTPUT_DIR}/64_one_variant.log"

VARIANTS=("22:17072086_A_C" "22:34636072_C_A" "22:51196592_T_C")
DATASET_SIZES=(100 500 1000 5000 10k)

for VARIANT_ID in "${VARIANTS[@]}"; do
  CHR=$(echo "$VARIANT_ID" | cut -d':' -f1)
  POS=$(echo "$VARIANT_ID" | cut -d':' -f2 | cut -d'_' -f1)
  REGION="${CHR}:${POS}-${POS}"
  CLEAN_ID=$(echo "$VARIANT_ID" | tr ':' '_' | tr -d '\r')

  for SIZE in "${DATASET_SIZES[@]}"; do
    if [[ "$SIZE" == "10k" ]]; then
      BASE_NAME="synthetic_10k"
    else
      BASE_NAME="synthetic_${SIZE}"
    fi

    # ========= Warm-up to normalize I/O =========
    [[ -f "${BASE_NAME}.tsv" ]] && head -n 100 "${BASE_NAME}.tsv" > /dev/null 2>&1
    [[ -f "${BASE_NAME}.tsv.gz" ]] && zcat "${BASE_NAME}.tsv.gz" | head -n 100 > /dev/null 2>&1
    [[ -f "${BASE_NAME}.vcf" ]] && head -n 100 "${BASE_NAME}.vcf" > /dev/null 2>&1
    [[ -f "${BASE_NAME}.vcf.gz" ]] && tabix "${BASE_NAME}.vcf.gz" "22:1-100" > /dev/null 2>&1

    # ==== Plain TSV ====
    FILE="${BASE_NAME}.tsv"
    [[ ! -f "$FILE" ]] && echo "Missing file: $FILE" && continue
    EXT_NAME="plain_tsv"
    OUT_FILE="${OUTPUT_DIR}/${CLEAN_ID}_${BASE_NAME}_${EXT_NAME}.tsv"
    TMP_ROW_FILE=$(mktemp)

    header=$(head -n 1 "$FILE")
    IFS=$'\t' read -r -a header_fields <<< "$header"
    for i in "${!header_fields[@]}"; do [[ "${header_fields[$i]}" == "ID" ]] && ID_COL=$((i + 1)) && break; done

    echo -e "\n=== File: $FILE | Query: $VARIANT_ID ===" >> "$TIME_LOG"
    { /usr/bin/time -v grep -P "\t$VARIANT_ID(\t|$)" "$FILE" | tail -n 1 > "$TMP_ROW_FILE"; } 2>> "$TIME_LOG"

    row=$(cat "$TMP_ROW_FILE")
    IFS=$'\t' read -r -a header_cols <<< "$header"
    IFS=$'\t' read -r -a row_vals <<< "$row"
    echo -e "ID\tStudy\tbeta\tse\tpval\tea_freq" > "$OUT_FILE"
    for i in "${!header_cols[@]}"; do
      col="${header_cols[$i]}"
      if [[ "$col" =~ ^(.+)\.ES$ ]]; then
        study="${BASH_REMATCH[1]}"
        beta="${row_vals[$i]}"
        se=""; pval=""; af=""
        for j in "${!header_cols[@]}"; do
          [[ "${header_cols[$j]}" == "${study}.SE" ]] && se="${row_vals[$j]}"
          [[ "${header_cols[$j]}" == "${study}.LP" ]] && pval=$(awk -v lp="${row_vals[$j]}" 'BEGIN{ printf "%.6g", 10^(-lp) }')
          [[ "${header_cols[$j]}" == "${study}.AF" ]] && af="${row_vals[$j]}"
        done
        echo -e "${VARIANT_ID}\t${study}\t${beta}\t${se}\t${pval}\t${af}" >> "$OUT_FILE"
      fi
    done
    rm "$TMP_ROW_FILE"

    # ==== Indexed TSV ====
    FILE="${BASE_NAME}.tsv.gz"
    [[ ! -f "$FILE" ]] && echo "Missing file: $FILE" && continue
    EXT_NAME="indexed_tsv"
    OUT_FILE="${OUTPUT_DIR}/${CLEAN_ID}_${BASE_NAME}_${EXT_NAME}.tsv"
    TMP_ROW_FILE=$(mktemp)

    header=$(zcat "$FILE" | head -n 1)
    IFS=$'\t' read -r -a header_fields <<< "$header"
    for i in "${!header_fields[@]}"; do [[ "${header_fields[$i]}" == "ID" ]] && ID_COL=$((i + 1)) && break; done

    echo -e "\n=== File: $FILE | Query: $VARIANT_ID ===" >> "$TIME_LOG"
    { /usr/bin/time -v tabix "$FILE" "$REGION" > "$TMP_ROW_FILE"; } 2>> "$TIME_LOG"

    row=$(cat "$TMP_ROW_FILE")
    IFS=$'\t' read -r -a header_cols <<< "$header"
    IFS=$'\t' read -r -a row_vals <<< "$row"
    echo -e "ID\tStudy\tbeta\tse\tpval\tea_freq" > "$OUT_FILE"
    for i in "${!header_cols[@]}"; do
      col="${header_cols[$i]}"
      if [[ "$col" =~ ^(.+)\.ES$ ]]; then
        study="${BASH_REMATCH[1]}"
        beta="${row_vals[$i]}"
        se=""; pval=""; af=""
        for j in "${!header_cols[@]}"; do
          [[ "${header_cols[$j]}" == "${study}.SE" ]] && se="${row_vals[$j]}"
          [[ "${header_cols[$j]}" == "${study}.LP" ]] && pval=$(awk -v lp="${row_vals[$j]}" 'BEGIN{ printf "%.6g", 10^(-lp) }')
          [[ "${header_cols[$j]}" == "${study}.AF" ]] && af="${row_vals[$j]}"
        done
        echo -e "${VARIANT_ID}\t${study}\t${beta}\t${se}\t${pval}\t${af}" >> "$OUT_FILE"
      fi
    done
    rm "$TMP_ROW_FILE"

    # ==== Plain VCF ====
    FILE="${BASE_NAME}.vcf"
    [[ ! -f "$FILE" ]] && echo "Missing file: $FILE" && continue
    EXT_NAME="plain_vcf"
    OUT_FILE="${OUTPUT_DIR}/${CLEAN_ID}_${BASE_NAME}_${EXT_NAME}.tsv"
    TMP_ROW_FILE=$(mktemp)

    header=$(grep -v '^##' "$FILE" | head -n 1)

    echo -e "\n=== File: $FILE | Query: $VARIANT_ID ===" >> "$TIME_LOG"
    { /usr/bin/time -v grep -P "^${CHR}\t${POS}\t" "$FILE" | grep -v '^##' > "$TMP_ROW_FILE"; } 2>> "$TIME_LOG"

    if [[ ! -s "$TMP_ROW_FILE" ]]; then
      echo "No VCF match found for $VARIANT_ID in $FILE" >> "$TIME_LOG"
      rm "$TMP_ROW_FILE"
      continue
    fi

    row=$(cat "$TMP_ROW_FILE")
    IFS=$'\t' read -r -a header_cols <<< "$header"
    IFS=$'\t' read -r -a row_vals <<< "$row"
    echo -e "ID\tStudy\tbeta\tse\tpval\tea_freq" > "$OUT_FILE"
    FORMAT_FIELD="${row_vals[8]}"
    IFS=':' read -r -a format_keys <<< "$FORMAT_FIELD"
    for i in $(seq 9 $((${#row_vals[@]} - 1))); do
      study="${header_cols[$i]}"
      IFS=':' read -r -a values <<< "${row_vals[$i]}"
      beta=""; se=""; pval=""; af=""
      for j in "${!format_keys[@]}"; do
        key="${format_keys[$j]}"
        [[ "$key" == "ES" ]] && beta="${values[$j]}"
        [[ "$key" == "SE" ]] && se="${values[$j]}"
        [[ "$key" == "LP" ]] && pval=$(awk -v lp="${values[$j]}" 'BEGIN{ printf "%.6g", 10^(-lp) }')
        [[ "$key" == "AF" ]] && af="${values[$j]}"
      done
      echo -e "${VARIANT_ID}\t${study}\t${beta}\t${se}\t${pval}\t${af}" >> "$OUT_FILE"
    done
    rm "$TMP_ROW_FILE"

  # ==== Indexed VCF ====
    FILE="${BASE_NAME}.vcf.gz"
    EXT_NAME="indexed_vcf"
    OUT_FILE="${OUTPUT_DIR}/${CLEAN_ID}_${BASE_NAME}_${EXT_NAME}.tsv"
    TMP_ROW_FILE=$(mktemp)

    header_line=$(bcftools view -h "$FILE" | grep -v '^##' | tail -n 1)
    IFS=$'\t' read -r -a header_cols <<< "$header_line"

    echo -e "\n=== File: $FILE | Query: $VARIANT_ID ===" >> "$TIME_LOG"

    START=$(date +%s%3N)
    bcftools view -H -r "$REGION" "$FILE" > "$TMP_ROW_FILE"
    END=$(date +%s%3N)
    ELAPSED_MS=$((END - START))
    echo "Elapsed (ms): $ELAPSED_MS" >> "$TIME_LOG"

    if [[ ! -s "$TMP_ROW_FILE" ]]; then
      echo "No variant data found for $VARIANT_ID in $FILE" >> "$TIME_LOG"
      rm "$TMP_ROW_FILE"
      continue
    fi

    row=$(cat "$TMP_ROW_FILE")
    IFS=$'\t' read -r -a row_vals <<< "$row"

    echo -e "ID\tStudy\tbeta\tse\tpval\tea_freq" > "$OUT_FILE"

    FORMAT_FIELD="${row_vals[8]}"
    IFS=':' read -r -a format_keys <<< "$FORMAT_FIELD"

    for i in $(seq 9 $((${#row_vals[@]} - 1))); do
      study="${header_cols[$i]}"
      IFS=':' read -r -a values <<< "${row_vals[$i]}"

      beta=""; se=""; pval=""; af=""
      for j in "${!format_keys[@]}"; do
        key="${format_keys[$j]}"
        [[ "$key" == "ES" ]] && beta="${values[$j]}"
        [[ "$key" == "SE" ]] && se="${values[$j]}"
        [[ "$key" == "LP" ]] && pval=$(awk -v lp="${values[$j]}" 'BEGIN{ printf "%.6g", 10^(-lp) }')
        [[ "$key" == "AF" ]] && af="${values[$j]}"
      done
      echo -e "${VARIANT_ID}\t${study}\t${beta}\t${se}\t${pval}\t${af}" >> "$OUT_FILE"
    done

    rm "$TMP_ROW_FILE"
  done
done
