import os
import pandas as pd

# Directories
input_dir = "/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/cleaned_final"
output_dir = "/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Required output columns in exact order
required_columns = ["#CHR", "POS", "ID", "EA", "NEA", "BETA", "SE", "PVAL", "EAF", "LP"]

# List all _cleaned.tsv files
cleaned_files = [f for f in os.listdir(input_dir) if f.endswith("_cleaned.tsv")]

# Process files in chunks to reduce memory usage
for file_name in cleaned_files:
    input_path = os.path.join(input_dir, file_name)
    output_path = os.path.join(output_dir, file_name.replace("_cleaned.tsv", "_ready.tsv"))
    
    print(f"📝 Processing file: {input_path}")
    
    try:
        chunk_size = 100000  # Process 100,000 rows at a time
        first_chunk = True

        chunk_iter = pd.read_csv(input_path, sep="\t", dtype=str, chunksize=chunk_size)
        
        for chunk in chunk_iter:
            # Ensure all required columns exist; if missing, create and fill with 'NA'
            for col in required_columns:
                if col not in chunk.columns:
                    chunk[col] = "NA"
            
            # Ensure columns are in the correct order
            chunk = chunk[required_columns]
            
            # Fill missing values with 'NA' explicitly
            chunk = chunk.fillna("NA")
            
            # Save the extracted data
            mode = "w" if first_chunk else "a"
            header = first_chunk  # Write header only for the first chunk
            chunk.to_csv(output_path, sep="\t", index=False, mode=mode, header=header)
            first_chunk = False  # Set to False after first write
        
        print(f"✅ Extracted file saved: {output_path}")
    
    except Exception as e:
        print(f"❌ Error processing {file_name}: {str(e)}")

print("🎉 Extraction complete for all cleaned files!")
