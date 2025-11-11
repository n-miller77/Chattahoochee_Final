#!/bin/bash
#SBATCH -J microbecensus
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 4
#SBATCH --mem=100G
#SBATCH -t2:00:00
#SBATCH -q inferno
#SBATCH --array=0-16%20
#SBATCH -o logs/microbecensus_%A_%a.out

# Load environment
eval "$(conda shell.bash hook)"
conda activate microbecensus_final



# Set parent directory
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"


SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" sample_ids4.txt)
SAMPLE_DIR="${PARENT_DIR}/${SAMPLE}"

# Skip directories with 'ds' in the name
if [[ "$SAMPLE_DIR" == *ds* ]]; then
    echo "Skipping directory with 'ds' in name: $SAMPLE_DIR"
    exit 0
fi

echo "Processing: $SAMPLE_DIR"

# Find the R1 fastq.gz file
R1_FILE=$(ls "$SAMPLE_DIR"/*R1*clean*.fastq.gz 2>/dev/null | head -n 1)

if [[ -z "$R1_FILE" ]]; then
    echo "  No *_R1_clean.fastq.gz found in $SAMPLE_DIR"
    exit 0
fi

# Define decompressed fastq
FASTQ_FILE="${R1_FILE%.gz}"

# Decompress only if not already present
if [[ -f "$FASTQ_FILE" ]]; then
    echo "  $FASTQ_FILE already exists, skipping decompression."
else
    echo "  Decompressing $R1_FILE → $FASTQ_FILE"
    gunzip -c "$R1_FILE" > "$FASTQ_FILE"
fi

# Define output file
DIRNAME=$(basename "$SAMPLE_DIR")
OUTPUT_FILE="${SAMPLE_DIR}/${DIRNAME}microbecensus_output_final.txt"

# Skip if output already exists
if [[ -f "$OUTPUT_FILE" ]]; then
    echo "  $OUTPUT_FILE already exists, skipping census."
    exit 0
fi

# Run microbe census
echo "  Running microbe census on $FASTQ_FILE..."
run_microbe_census.py -n 100000000000000 -t 8 "$FASTQ_FILE" "$OUTPUT_FILE"
