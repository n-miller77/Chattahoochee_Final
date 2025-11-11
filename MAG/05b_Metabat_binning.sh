#!/bin/bash
#SBATCH -J abundance
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 4
#SBATCH --mem=64GB
#SBATCH -t 24:00:00
#SBATCH -q inferno
#SBATCH --array=0-45%20
#SBATCH -o logs/metabat2_nohuman_%A_%a.out


# Load conda
eval "$(conda shell.bash hook)"
conda activate MetaBat



# Base directory
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"


SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" sample_ids2.txt)
SAMPLE_DIR="${PARENT_DIR}/${SAMPLE}"

# Skip directories with 'ds' in the name
if [[ "$SAMPLE_DIR" == *ds* ]]; then
    echo "Skipping directory with 'ds' in name: $SAMPLE_DIR"
    exit 0
fi

echo "Processing: $SAMPLE_DIR"

# Check if analysis output already exists
METABAT_DIR="${SAMPLE_DIR}/analysis_nohuman_again"
OUT_PREFIX="${METABAT_DIR}/BinningMeta"
if compgen -G "${OUT_PREFIX}*.fa" > /dev/null; then
    echo "  MetaBat2 bins already exist, skipping."
    exit 0
fi

# Find normal2* directory
SPADES_DIR=$(find "$SAMPLE_DIR" -maxdepth 1 -type d -name "nohuman_spades_output_*" | head -n 1)
if [[ -z "$SPADES_DIR" ]]; then
    echo "  No normal2* directory found, skipping."
    exit 0
fi

CONTIGS="${SPADES_DIR}/contigs.fasta"
if [[ ! -f "$CONTIGS" ]]; then
    echo "  contigs.fasta not found in $SPADES_DIR, skipping."
    exit 0
fi

COVERAGE_FILE="${SAMPLE_DIR}/analysis_nohuman_again/coverage_for_metabat2.txt"
if [[ ! -f "$COVERAGE_FILE" ]]; then
    echo "  coverage_for_metabat.txt not found, skipping."
    exit 0
fi

mkdir -p "$METABAT_DIR"

echo "  Running MetaBat2..."
metabat2 \
    -i "$CONTIGS" \
    -o "$OUT_PREFIX" \
    -a "$COVERAGE_FILE" \
    -t 6

echo "Done with $SAMPLE_DIR"
echo "-------------------------------------------"









###############interation
#!/bin/bash
#SBATCH -J Metabatnohuman
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 6
#SBATCH --mem=50G
#SBATCH -t2:00:00
#SBATCH -q inferno


# Load environment
eval "$(conda shell.bash hook)"
conda activate MetaBat

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

    # Check if analysis output already exists
    METABAT_DIR="${SAMPLE_DIR}/analysis_nohuman_again"
    OUT_PREFIX="${METABAT_DIR}/BinningMeta"
    if compgen -G "${OUT_PREFIX}*.fa" > /dev/null; then
        echo "  MetaBat2 bins already exist, skipping."
        continue
    fi

    # Find normal2* directory
    SPADES_DIR=$(find "$SAMPLE_DIR" -maxdepth 1 -type d -name "nohuman_spades_output_*" | head -n 1)
    if [[ -z "$SPADES_DIR" ]]; then
        echo "  No normal2* directory found, skipping."
        continue
    fi

    CONTIGS="${SPADES_DIR}/contigs.fasta"
    if [[ ! -f "$CONTIGS" ]]; then
        echo "  contigs.fasta not found in $SPADES_DIR, skipping."
        continue
    fi

    COVERAGE_FILE="${SAMPLE_DIR}/analysis_nohuman_again/coverage_for_metabat2.txt"
    if [[ ! -f "$COVERAGE_FILE" ]]; then
        echo "  coverage_for_metabat.txt not found, skipping."
        continue
    fi

    mkdir -p "$METABAT_DIR"

    echo "  Running MetaBat2..."
    metabat2 \
        -i "$CONTIGS" \
        -o "$OUT_PREFIX" \
        -a "$COVERAGE_FILE" \
        -t 6

    echo "Done with $SAMPLE_DIR"
    echo "-------------------------------------------"
done
