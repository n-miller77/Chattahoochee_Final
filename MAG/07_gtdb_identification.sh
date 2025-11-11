#!/bin/bash
#SBATCH -Jgtdbtk
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 10
#SBATCH --mem=200G
#SBATCH -t24:00:00
#SBATCH -q inferno
#SBATCH --array=0-45%20
#SBATCH -o logs/gtdbtk_nohuman_%A_%a.out


eval "$(conda shell.bash hook)"
conda activate gtdb



# Set the parent directory where your samples live
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"
GTDB_OUTPUT_ROOT="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/gtdb_concat_output_final_final"
GTDB_DATA_PATH="/storage/project/r-ktk3-0/shared/rich_project_bio-konstantinidis/shared3/release226"
EXT="fasta"  # File extension for input genomes

# Set GTDB data path
export GTDBTK_DATA_PATH="$GTDB_DATA_PATH"

# Get sample from file
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" sample_ids2.txt)
SAMPLE_NAME="$SAMPLE"
SAMPLE_DIR="${PARENT_DIR}/${SAMPLE_NAME}"
GENOME_DIR="${SAMPLE_DIR}/analysis_nohuman_final/Dereplicated/dereplicated_genomes"

# Check if genome dir exists
if [[ ! -d "$GENOME_DIR" ]]; then
    echo "$SAMPLE_NAME: genome dir not found — skipping"
    exit 0
fi

# Check if genome files with correct extension exist
GENOME_COUNT=$(find "$GENOME_DIR" -maxdepth 1 -type f -name "*.$EXT" | wc -l)
if [[ "$GENOME_COUNT" -eq 0 ]]; then
    echo "$SAMPLE_NAME: no .$EXT files in genome dir — skipping"
    exit 0
fi

# Define and create output directory
OUT_DIR="${GTDB_OUTPUT_ROOT}/${SAMPLE_NAME}"
mkdir -p "$OUT_DIR"

echo "Running GTDB-Tk for $SAMPLE_NAME with $GENOME_COUNT genome(s)..."

gtdbtk classify_wf \
    --genome_dir "$GENOME_DIR" \
    -x "$EXT" \
    --out_dir "$OUT_DIR" \
    --cpus 10 \
    --pplacer_cpus 10 \
    --mash_db "$GTDBTK_DATA_PATH"

echo "✅ Done: $SAMPLE_NAME → $OUT_DIR"
