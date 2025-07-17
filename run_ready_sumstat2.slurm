#!/bin/bash -l
#SBATCH --job-name=harmonize_gwas
#SBATCH --partition=cpu
#SBATCH --output=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/pipeline_%j.out
#SBATCH --error=/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/ready_backup/pipeline_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --ntasks=1
#SBATCH --nodes=1



echo "Starting job at $(date)"
echo "Running on host: $(hostname)"

# ====== Load the correct Python interpreter directly ======
PYTHON_EXEC="/users/k24056829/.conda/envs/pythonday/bin/python"
SCRIPT="/scratch/prj/mmg_data/gwas_catalog/downloads/HarmonizingExampleGwas/cleaned_final/ready.py"

if [ ! -x "$PYTHON_EXEC" ]; then
    echo "Python executable not found at: $PYTHON_EXEC"
    exit 1
fi

if [ ! -f "$SCRIPT" ]; then
    echo "Python script not found at: $SCRIPT"
    exit 1
fi

echo "Using Python from: $PYTHON_EXEC"
echo "Running script: $SCRIPT"

# ====== Run the script ======
$PYTHON_EXEC "$SCRIPT"

echo "Job completed at $(date)"
