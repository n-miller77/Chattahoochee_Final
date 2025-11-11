#!/bin/bash
#SBATCH -J abundance
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 4
#SBATCH --mem=64GB
#SBATCH -t 24:00:00
#SBATCH -q inferno
#SBATCH --array=0-45%20
#SBATCH -o logs/maxbin_nohuman_%A_%a.out


# Load conda
eval "$(conda shell.bash hook)"
conda activate Maxbin2



# Base directory
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"

# Get sample directory from array file
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" sample_ids2.txt)
SAMPLE_DIR="${PARENT_DIR}/${SAMPLE}"

# Skip directories with 'ds' in the name
if [[ "$SAMPLE_DIR" == *ds* ]]; then
    echo "Skipping directory with 'ds' in name: $SAMPLE_DIR"
    exit 0
fi

echo "Processing: $SAMPLE_DIR"

# Define analysis directory and output prefix
ANALYSIS_DIR="${SAMPLE_DIR}/analysis_nohuman_again"
OUT_PREFIX="${ANALYSIS_DIR}/BinningMax"

# Skip if bins already exist
if compgen -G "${OUT_PREFIX}*.fasta" > /dev/null; then
    echo "  MaxBin2 bins already exist, skipping."
    exit 0
fi

# Find new_spades_output* directory
SPADES_DIR=$(find "$SAMPLE_DIR" -maxdepth 1 -type d -name "nohuman_spades_output_*" | head -n 1)
if [[ -z "$SPADES_DIR" ]]; then
    echo "  No new_spades_output* directory found, skipping."
    exit 0
fi

CONTIGS="${SPADES_DIR}/contigs.fasta"
if [[ ! -f "$CONTIGS" ]]; then
    echo "  contigs.fasta not found in $SPADES_DIR, skipping."
    exit 0
fi

COVERAGE_FILE="${SAMPLE_DIR}/analysis_nohuman_again/coverage_for_maxbin2.txt"
if [[ ! -f "$COVERAGE_FILE" ]]; then
    echo "  coverage_for_maxbin.txt not found, skipping."
    exit 0
fi

mkdir -p "$ANALYSIS_DIR"

echo "  Running MaxBin2..."
run_MaxBin.pl \
    -contig "$CONTIGS" \
    -abund "$COVERAGE_FILE" \
    -out "$OUT_PREFIX" \
    -thread 6

echo "Done with $SAMPLE_DIR"
echo "-------------------------------------------"







###############iteration
#!/bin/bash
#SBATCH -J Maxbin
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 6
#SBATCH --mem=50G
#SBATCH -t2:00:00
#SBATCH -q inferno


# Load environment
eval "$(conda shell.bash hook)"
conda activate Maxbin2

# Set parent directory
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"


# Loop through each sample directory
for SAMPLE_DIR in "$PARENT_DIR"/*/; do
    # Skip directories with 'ds' in the name
    if [[ "$SAMPLE_DIR" == *ds* ]]; then
        echo "Skipping directory with 'ds' in name: $SAMPLE_DIR"
        continue
    fi

    echo "Processing: $SAMPLE_DIR"

    # Define analysis directory and output prefix
    ANALYSIS_DIR="${SAMPLE_DIR}/analysis_nohuman_again"
    OUT_PREFIX="${ANALYSIS_DIR}/BinningMax"

    # Skip if bins already exist
    if compgen -G "${OUT_PREFIX}*.fasta" > /dev/null; then
        echo "  MaxBin2 bins already exist, skipping."
        continue
    fi

    # Find new_spades_output* directory
    SPADES_DIR=$(find "$SAMPLE_DIR" -maxdepth 1 -type d -name "nohuman_spades_output_*" | head -n 1)
    if [[ -z "$SPADES_DIR" ]]; then
        echo "  No new_spades_output* directory found, skipping."
        continue
    fi

    CONTIGS="${SPADES_DIR}/contigs.fasta"
    if [[ ! -f "$CONTIGS" ]]; then
        echo "  contigs.fasta not found in $SPADES_DIR, skipping."
        continue
    fi

    COVERAGE_FILE="${ANALYSIS_DIR}/coverage_for_maxbin2.txt"
    if [[ ! -f "$COVERAGE_FILE" ]]; then
        echo "  coverage_for_maxbin.txt not found, skipping."
        continue
    fi

    mkdir -p "$ANALYSIS_DIR"

    echo "  Running MaxBin2..."
    run_MaxBin.pl \
        -contig "$CONTIGS" \
        -abund "$COVERAGE_FILE" \
        -out "$OUT_PREFIX" \
        -thread 6

    echo "Done with $SAMPLE_DIR"
    echo "-------------------------------------------"
done
