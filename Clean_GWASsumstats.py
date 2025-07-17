import os
import pandas as pd
import numpy as np

data_dir = "/scratch/prj/mmg_data/gwas_catalog/downloads"
output_dir = "/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/cleaned_final"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# List all .tsv files, excluding already processed ones
tsv_files = [f for f in os.listdir(data_dir) if f.endswith(".tsv")]

# Define column mappings
expected_mapping = {
    "SNP": ["variant_id","rs_id", "RSID", "rsids", "SNPID", "Name", "rsid", "rsID", "identifier", "SNP"],
    "#CHR": ["chromosome", "Chromosome", "CHR", "#CHROM", "CHROM"],
    "POS": ["base_pair_location", "Position", "POS", "BP", "GENPOS", "pos"],
    "EA": ["effect_allele", "A1", "ALT", "EA", "Effect_allele"],
    "NEA": ["Other_allele","other_allele", "A2", "REF", "ref_allele", "NEA"],
    "BETA": ["beta", "BETA", "Beta", "log_odds"],
    "SE": ["standard_error", "SE", "se", "sebeta"],
    "PVAL": ["pval","P-value", "PVAL", "P", "p.value", "P_VALUE", "p_value", "PVALUE"],
    "EAF": ["effect_allele_frequency", "af", "AF_ALT", "eaf", "MAF", "EAF", "A1.AF", "af_alt", "AF"],
    "INFO": ["INFO"],
    "N": ["sample_size", "N", "n", "European_N", "AllEthnicities_N"]
}



required_columns = ["#CHR", "POS", "ID", "EA", "NEA", "BETA", "SE", "PVAL", "SNP", "EAF"]

def process_chunk(chunk, file_name):
    try:
        chunk.columns = chunk.columns.str.strip()

        # Rename P-value columns
        if "p_value" in chunk.columns:
            chunk.rename(columns={"p_value": "P"}, inplace=True)
        if "P" in chunk.columns:
            chunk.rename(columns={"P": "PVAL"}, inplace=True)

        # Rename based on expected mapping
        columns_to_rename = {old: new for new, olds in expected_mapping.items() for old in olds if old in chunk.columns}
        chunk.rename(columns=columns_to_rename, inplace=True)

        # Convert EA and NEA to uppercase
        if "EA" in chunk.columns:
            chunk["EA"] = chunk["EA"].str.upper()
        if "NEA" in chunk.columns:
            chunk["NEA"] = chunk["NEA"].str.upper()

        # Check if SE might actually be PVAL
        if "PVAL" not in chunk.columns and "SE" in chunk.columns:
            se_values = pd.to_numeric(chunk["SE"], errors="coerce")
            if se_values.min() < 1e-5:
                print(f"⚠️ Detected possible P-values in SE column for {file_name}. Moving SE to PVAL.")
                chunk.rename(columns={"SE": "PVAL"}, inplace=True)
                chunk["SE"] = "NA"

        # Ensure all required columns exist
        for col in required_columns:
            if col not in chunk.columns:
                chunk[col] = "NA"

        # Generate unique ID if missing
        if all(col in chunk.columns for col in ["#CHR", "POS", "NEA", "EA"]):
            chunk["ID"] = chunk.apply(lambda row: f"{row['#CHR']}:{row['POS']}_{row['NEA']}_{row['EA']}", axis=1)

        # Step 5: Move gene-like values to GENE column
        if "SNP" in chunk.columns:
            gene_like = chunk["SNP"].astype(str).str.contains(r"ENSG\d+|GENE", regex=True, na=False)
            if gene_like.any():
                print("⚠️ Detected gene annotations in 'SNP' column. Moving to 'GENE' column.")
                if "GENE" not in chunk.columns:
                    chunk["GENE"] = "NA"
                chunk.loc[gene_like, "GENE"] = chunk.loc[gene_like, "SNP"]
                chunk.loc[gene_like, "SNP"] = "NA"

        # ✅ Step 6: Add LP column = -log10(PVAL)
        if "PVAL" in chunk.columns:
            chunk["LP"] = pd.to_numeric(chunk["PVAL"], errors="coerce").apply(
                lambda p: -1 * np.log10(p) if pd.notnull(p) and p > 0 else "."
            )
        else:
            chunk["LP"] = "."

        return chunk

    except Exception as e:
        print(f"❌ Error processing {file_name}: {str(e)}")
        return None

# Process files
for file_name in tsv_files:
    file_path = os.path.join(data_dir, file_name)
    cleaned_file_path = os.path.join(output_dir, file_name.replace(".tsv", "_cleaned.tsv"))

    print(f"📝 Processing file: {file_path}")

    try:
        # Read first few lines to check if headers exist
        sample_df = pd.read_csv(file_path, sep="\t", engine="c", skipinitialspace=True, nrows=5, header=None)
        has_headers = all(isinstance(col, str) for col in sample_df.iloc[0])

        # If no headers, assign default column names
        if not has_headers:
            print(f"⚠️ No headers detected in {file_name}. Assigning default column names once.")
            num_cols = sample_df.shape[1]
            default_headers = ["SNP", "#CHR", "POS", "EA", "NEA", "BETA", "SE", "PVAL", "INFO", "EAF"][:num_cols]
            chunk_iter = pd.read_csv(file_path, sep="\t", engine="c", skipinitialspace=True, low_memory=False, chunksize=100000, header=None, names=default_headers)
        else:
            chunk_iter = pd.read_csv(file_path, sep="\t", engine="c", skipinitialspace=True, low_memory=False, chunksize=100000, header=0)

        first_chunk = True
        for chunk in chunk_iter:
            cleaned_chunk = process_chunk(chunk, file_name)
            if cleaned_chunk is not None:
                mode = "w" if first_chunk else "a"
                header = first_chunk  # Write header only for the first chunk
                cleaned_chunk.to_csv(cleaned_file_path, sep="\t", index=False, mode=mode, header=header)
                first_chunk = False  # Set to False after first write

        print(f"✅ Cleaned file saved: {cleaned_file_path}")

    except Exception as e:
        print(f"❌ Error processing {file_name}: {str(e)}")

print("🎉 Processing complete for all .tsv files!")

