#!/bin/bash
#SBATCH --job-name=bbmap
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH -t 72:00:00
#SBATCH -q inferno
#SBATCH --array=0-185%20
#SBATCH -o logs/bbmap_step1_%A_%a.out



# Load environment
eval "$(conda shell.bash hook)"
conda activate bbmap



set -euo pipefail
shopt -s nullglob


# ==== EDIT THESE ====
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work"



SAMPLE_NAME="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" sample_ids.txt || true)"
if [[ -z "${SAMPLE_NAME:-}" ]]; then
  echo "No entry for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
  exit 0
fi

SAMPLE_DIR="${PARENT_DIR}/${SAMPLE_NAME}"
if [[ ! -d "$SAMPLE_DIR" ]]; then
  echo "❌ Sample dir not found: $SAMPLE_DIR"
  exit 1
fi

echo "▶ Processing: $SAMPLE_NAME"

# --- Find inputs (first matching pair) ---
R1=$(find "$SAMPLE_DIR" -maxdepth 1 -type f -name "*_R1*.fastq.gz" | head -n 1)
R2=$(find "$SAMPLE_DIR" -maxdepth 1 -type f -name "*_R2*.fastq.gz" | head -n 1)

if [[ -z "$R1" || -z "$R2" ]]; then
  echo "❌ Missing R1/R2 in $SAMPLE_DIR"
  exit 1
fi

# --- Outputs ---
OUT1="${SAMPLE_DIR}/${SAMPLE_NAME}_R1_clean.fastq.gz"
OUT2="${SAMPLE_DIR}/${SAMPLE_NAME}_R2_clean.fastq.gz"


# Skip if already done (non-empty files)
if [[ -s "$OUT1" && -s "$OUT2" ]]; then
  echo "⏭️  Skipping $SAMPLE_NAME: outputs already exist"
  exit 0
fi

# Use temp names to avoid partial files appearing complete
TMP1="${OUT1}.tmp"
TMP2="${OUT2}.tmp"
rm -f "$TMP1" "$TMP2"

# --- Run BBDuk ---
echo "🚀 Running BBDuk on $SAMPLE_NAME"
bbduk.sh \
  in="$R1" in2="$R2" \
  out="$TMP1" out2="$TMP2" \
  qtrim=w trimq=17 minlength=70 tbo tossjunk=t cardinalityout=t

# Promote to final names on success
mv -f "$TMP1" "$OUT1"
mv -f "$TMP2" "$OUT2"

echo "✅ Done: $SAMPLE_NAME -> $OUT1, $OUT2"
