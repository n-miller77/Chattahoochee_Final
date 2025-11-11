#!/bin/bash
#SBATCH -J sourceapp2_loop
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t12:00:00
#SBATCH -q inferno
#SBATCH --array=1-45%20
#SBATCH -o logs/sourceapp_%A_%a.out


#### Download metaG ####
eval "$(conda shell.bash hook)"
conda activate spades

set -euo pipefail

# Paths
READS_DIR="/scicomp/home/reads"
KRAKEN2_BIN="/scicomp/home/kraken2/kraken2"
DB="/scicomp/databases/kraken/k2_standard_full"
OUTDIR="/scicomp/home/data/kraken2"
SAMPLE_LIST="/scicomp/home/sample_list.txt"  # full path to your sample list

# Get the sample ID for this array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

# File paths
R1="${READS_DIR}/${SAMPLE}/${SAMPLE}_R1_001.fastq.gz"
R2="${READS_DIR}/${SAMPLE}/${SAMPLE}_R2_001.fastq.gz"

# Output paths
mkdir -p "$OUTDIR"

echo "Running Kraken2 for sample: $SAMPLE"
echo "R1: $R1"
echo "R2: $R2"

"$KRAKEN2_BIN" \
  --db "$DB" \
  --threads "$SLURM_CPUS_PER_TASK" \
  --gzip-compressed \
  --paired \
  --unclassified-out "$OUTDIR/${SAMPLE}.unclassified.#.fq" \
  --classified-out "$OUTDIR/${SAMPLE}.classified.#.fq" \
  --report "$OUTDIR/${SAMPLE}.kreport" \
  --output "$OUTDIR/${SAMPLE}.kout" \
  "$R1" "$R2"










#!/usr/bin/env bash
#SBATCH -Jindex
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 1
#SBATCH --mem=32G
#SBATCH -t6:00:00
#SBATCH -q inferno


# Load environment
eval "$(conda shell.bash hook)"
conda activate bracken


set -euo pipefail

# Paths
KREPORT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/kraken_output"        # where your *.kreport files live
DB="/scicomp/home/databases/kraken/k2_standard_full"   # path to kraken2/bracken database
OUTDIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/kraken_output/bracken_output"

# Make sure output directory exists
mkdir -p "$OUTDIR"

echo "Running Bracken on all kreport files in: $KREPORT_DIR"

for KREPORT in "$KREPORT_DIR"/*.kreport; do
    BASENAME=$(basename "$KREPORT")
    SAMPLE="${BASENAME%.kreport}"

    echo "Processing sample: $SAMPLE"

    # Species level
    bracken \
        -t 32 \
        -d "$DB" \
        -i "$KREPORT" \
        -o "$OUTDIR/${SAMPLE}.bracken.S.out" \
        -w "$OUTDIR/${SAMPLE}.bracken.kraken.S.out"

    # Genus level
    bracken \
        -t 32 \
        -d "$DB" \
        -i "$KREPORT" \
        -l G \
        -o "$OUTDIR/${SAMPLE}.bracken.G.out" \
        -w "$OUTDIR/${SAMPLE}.bracken.kraken.G.out"

done

echo "All samples processed."
