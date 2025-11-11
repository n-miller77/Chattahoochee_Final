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









######concatnate cleaned files
#!/bin/bash

# Set the parent directory (you can change this to an argument if needed)
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work"
OUTPUT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Step 1: Find all directories matching the pattern *_L00*
find "$PARENT_DIR" -type d -name "*_L00*" | while read -r dir; do
    # Get the base name of the directory
    base=$(basename "$dir")

    # Extract the prefix before _L00
    prefix=${base%%_L00*}

    # Add to an associative array (simulate using temp files for portability)
    echo "$dir" >> "/tmp/${prefix}_dirs.txt"
done

# Step 2: Loop through each prefix group
for file in /tmp/*_dirs.txt; do
    # Get the prefix from the filename
    prefix=$(basename "$file" | sed 's/_dirs.txt//')

    # Make the output subdirectory
    out_subdir="${OUTPUT_DIR}/${prefix}"
    mkdir -p "$out_subdir"

    # Output files
    out_R1="${out_subdir}/${prefix}_R1_clean.fastq.gz"
    out_R2="${out_subdir}/${prefix}_R2_clean.fastq.gz"

    echo "Processing group: $prefix"

    # Empty the output files if they exist
    : > "$out_R1"
    : > "$out_R2"

    # Loop through each directory in the group
    while read -r subdir; do
        echo "  Looking in: $subdir"

        # Find R1 files
        for r1file in "$subdir"/*R1_clean.fastq.gz; do
            [ -f "$r1file" ] || continue
            echo "    Adding R1: $r1file"
            cat "$r1file" >> "$out_R1"
        done

        # Find R2 files
        for r2file in "$subdir"/*R2_clean.fastq.gz; do
            [ -f "$r2file" ] || continue
            echo "    Adding R2: $r2file"
            cat "$r2file" >> "$out_R2"
        done
    done < "$file"

    # Optionally compress if not already compressed
    # gzip "$out_R1"
    # gzip "$out_R2"

    # Cleanup temp file
    rm -f "$file"
done

echo "✅ Done. Concatenated files are in: $OUTPUT_DIR
