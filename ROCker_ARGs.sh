#!/bin/bash
#SBATCH -J rocker_array
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 10
#SBATCH --mem=100G
#SBATCH -t 24:00:00
#SBATCH -q inferno
#SBATCH --array=0-45%20
#SBATCH -o logs/rocker_concat_MBL_addtl_%A_%a.out

# Load environment
eval "$(conda shell.bash hook)"
conda activate gem

# Set parent directory
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/scratch/Chattahooche-samples-working/"

# Generate sample_ids.txt if it doesn't exist
if [[ ! -f sample_ids2.txt ]]; then
    echo "Generating sample_ids.txt from $PARENT_DIR"
    find "$PARENT_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort > sample_ids2.txt
fi

# Get sample ID for this array index
SAMPLE_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" sample_ids2.txt)
SAMPLE_DIR="${PARENT_DIR}/${SAMPLE_ID}"

echo "Processing: $SAMPLE_ID"


OUTPUT_FILE="${SAMPLE_DIR}/rocker_concat_BLAA.blast"
if [[ -s "$OUTPUT_FILE" ]]; then
    echo "Output file already exists for $SAMPLE_ID, skipping..."
    exit 0
fi


# Define input FASTA path
INPUT_FASTA="${SAMPLE_DIR}/R1_for_rocker.fasta"

# Check for existing and non-empty FASTA
if [[ ! -s "$INPUT_FASTA" ]]; then
    echo "FASTA missing or empty, regenerating from FASTQ..."

    # Find the corresponding FASTQ file
    FASTQ_FILE=$(find "$SAMPLE_DIR" -name "*_R1*clean.fastq.gz" | head -n 1)

    if [[ -f "$FASTQ_FILE" ]]; then
        echo "Found FASTQ: $FASTQ_FILE"
        seqtk seq -A "$FASTQ_FILE" > "$INPUT_FASTA"

        # Double-check output
        if [[ ! -s "$INPUT_FASTA" ]]; then
            echo "Error: FASTA file is still empty after conversion!"
            exit 1
        fi
    else
        echo "No FASTQ file found in $SAMPLE_DIR to regenerate FASTA."
        exit 1
    fi
fi

# Run ROCker
ruby /storage/home/hcoda1/9/nmiller304/rocker/bin/ROCker search \
    -i "$INPUT_FASTA" \
    -k /storage/home/hcoda1/9/nmiller304/rocker/models/BLA/BlaA.150.rocker \
    -o "${SAMPLE_DIR}/rocker_concat_BLAA.blast"

echo "Finished: $SAMPLE_ID"
