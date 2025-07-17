#!/bin/bash

DATA_DIR="/scratch/prj/mmg_data/gwas_catalog/downloads/"
cd "$DATA_DIR" || { echo "❌ ERROR: Directory not found!"; exit 1; }

echo "📂 Processing sumstat files in: $DATA_DIR"

# Step 1: Unzip all .gz files
echo "🔍 Checking for .gz compressed files..."
for file in *.gz; do
    [[ -f "$file" ]] && echo "📦 Unzipping: $file" && gunzip "$file"
done

# Step 2: Convert .txt files to .tsv (keep originals)
echo "🔍 Checking for .txt files to convert..."
for file in *.txt; do
    [[ -f "$file" ]] || continue
    new_name="${file%.txt}.tsv"
    echo "🔄 Converting $file → $new_name (keeping original)"

    # Convert to tab-delimited and write to new .tsv file
    awk '{$1=$1}1' OFS='\t' "$file" > "$new_name"
done

echo "✅ All files are now ready for processing by the Python script!"
