#!/bin/bash
#SBATCH -J abundance
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 4
#SBATCH --mem=64GB
#SBATCH -t 24:00:00
#SBATCH -q inferno
#SBATCH --array=0-45%20
#SBATCH -o logs/drep_nohuman_%A_%a.out


# Load conda
eval "$(conda shell.bash hook)"
conda activate dRep2


# Base directory
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"



SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" sample_ids2.txt)
SAMPLE_DIR="${PARENT_DIR}/${SAMPLE}"

# Skip ds directories
if [[ "$SAMPLE_DIR" == *ds* ]]; then
    echo "Skipping directory with 'ds' in name: $SAMPLE_DIR"
    exit 0
fi

echo "📂 Processing sample: $SAMPLE_DIR"

###############################################
###############################################
# Step 1: Rename .fa → .fasta if present
###############################################
#mkdir -p "$SAMPLE_DIR/analysis_nohuman_again/Binning_files"

MAXBIN_DIR="$SAMPLE_DIR/analysis_nohuman_again"
BINNING_DIR="$MAXBIN_DIR/Binning_files_backup"

# ✅ Skip if .fasta files already in Binning_files
if compgen -G "$BINNING_DIR/*.fasta" > /dev/null; then
    echo "✅ .fasta files already exist in Binning_files — skipping rename and move."
else
    if [[ -d "$MAXBIN_DIR" ]]; then
        echo "🔍 Checking for .fa files in $MAXBIN_DIR"
        for FILE in "$MAXBIN_DIR"/*.fa; do
            [ -e "$FILE" ] || continue
            NEW_FILE="${FILE%.fa}.fasta"
            mv "$FILE" "$NEW_FILE"
            echo "✅ Renamed: $(basename "$FILE") -> $(basename "$NEW_FILE")"
        done

        echo "📦 Moving .fasta files to Binning_files..."
        mv "$MAXBIN_DIR"/*.fasta "$BINNING_DIR"/
    else
        echo "⚠️ No MaxBin directory in $SAMPLE_DIR"
    fi
fi




###############################################
# Step 2: Run CheckM2
###############################################


export CHECKM2DB="/storage/home/hcoda1/9/nmiller304/databases/CheckM2_database/uniref100.KO.1"



INPUT_DIR="$SAMPLE_DIR/analysis_nohuman_again"
CHECKM_OUT="$SAMPLE_DIR/checkm2_output_again_againnn"

if [[ -d "$INPUT_DIR" ]]; then
    echo "🧪 Running CheckM2 on $INPUT_DIR"
    mkdir -p "$CHECKM_OUT"
    checkm2 predict \
        --input "$INPUT_DIR/Binning_files_backup" \
        --output-directory "$CHECKM_OUT" \
        -x .fasta \
        --force \
        --threads 10
else
    echo "⚠️ Skipping CheckM2 — no analysis directory in $SAMPLE_DIR"
fi

###############################################
# Step 3: Run dRep dereplication
###############################################


QUALITY_FILE="$CHECKM_OUT/quality_report.tsv"
CSV_FILE="$CHECKM_OUT/CheckM2_table.csv"
FIXED_CSV="$CHECKM_OUT/CheckM2_table_fix.csv"
OUTPUT_DIR="$INPUT_DIR/Dereplicated_take3"

if [[ ! -f "$QUALITY_FILE" ]]; then
    echo "❌ Missing quality_report.tsv in $CHECKM_OUT, skipping dRep"
    exit 0
fi

echo "📄 Converting CheckM2 quality report for dRep"
awk 'BEGIN {OFS=","; print "genome,completeness,contamination"} NR>1 {print $1,$2,$3}' "$QUALITY_FILE" > "$CSV_FILE"
awk -F',' 'NR==1 {print; next} { $1 = $1 ".fasta"; print }' OFS=',' "$CSV_FILE" > "$FIXED_CSV"

echo "🚀 Running dRep dereplication"
dRep dereplicate "$OUTPUT_DIR" \
    --genomes "$INPUT_DIR"/Binning_files_backup/*.fasta \
    --genomeInfo "$FIXED_CSV" \
    --S_algorithm fastANI \
    -sa 0.95 \
    -comp 50 \
    -p 10

echo "✅ Finished $SAMPLE_DIR"
echo "-------------------------------------------"
