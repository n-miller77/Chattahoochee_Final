#!/bin/bash
#SBATCH --job-name=metaspades_array
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH -t 72:00:00
#SBATCH -q inferno
#SBATCH --array=0-45%20
#SBATCH -o logs/metaspades_nohuman_%A_%a.out

# Load conda
eval "$(conda shell.bash hook)"
conda activate spades

# Base directory
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"

# Get sample name from array ID
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" sample_ids2.txt)
SAMPLE_ID="${SAMPLE}"

# Navigate to subdirectory
SAMPLE_DIR="${PARENT_DIR}/${SAMPLE}"
cd "$SAMPLE_DIR" || { echo "Directory $SAMPLE_DIR not found"; exit 1; }

# Find paired reads
R1=$(ls *_R1*clean_removed.fastq.gz 2>/dev/null)
R2=$(ls *_R2*clean_removed.fastq.gz 2>/dev/null)

if [[ -z "$R1" || -z "$R2" ]]; then
  echo "Missing reads in $SAMPLE_DIR"
  exit 1
fi

# Find existing spades_output* directory if it exists
SPADES_DIR=$(find . -maxdepth 1 -type d -name "spades_output*" | head -n 1)

# If contigs.fasta already exists in a metaspades_3_5973873_8.outspades_output* directory, skip
if [[ -n "$SPADES_DIR" && -f "$SPADES_DIR/contigs.fasta" ]]; then
  echo "$SAMPLE_ID already has contigs in $SPADES_DIR. Skipping."
  exit 0
fi

# Otherwise, define a new output directory name for metaSPAdes
OUTDIR="${SAMPLE_DIR}/nohuman_spades_output_${SAMPLE_ID}"

echo "Running metaSPAdes for $SAMPLE_ID"
metaspades.py \
  -1 "$R1" -2 "$R2" \
  -o "$OUTDIR" \
  -m 400 -t 16 -k 11,33,55,77,127

echo "Done with $SAMPLE_ID"
