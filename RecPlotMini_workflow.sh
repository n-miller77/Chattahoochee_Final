    echo "Running BLAST for ${SAMPLE} against combined rMAG database..."
    blastn \
        -query "${FASTA_FILE}" \
        -db "${COMBINED_DB_BASE}" \
        -evalue 0.01 \
        -num_threads 20 \
        -max_target_seqs 10 \
        -perc_identity 70 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
        -out "${COMP_BLAST_RAW}"




        
# Loop through all subdirectories
for dir in */ ; do
    # Look for a .blast.tsv file in the subdirectory
    blast_file=$(find "$dir" -maxdepth 1 -type f -name "*.blast.tsv")

    # If a blast file exists, process it
    if [[ -n "$blast_file" ]]; then
        # Create output filename: insert _filtered before .blast.tsv
        filtered_file="${blast_file%.blast.tsv}_filtered.blast.tsv"

        echo "Processing: $blast_file"
        echo "Output:     $filtered_file"

        # Filter based on alignment length ≥100 and coverage ≥90%
        awk '
        {
            len = $4; qlen = $13;
            if (len >= 100 && len/qlen >= 0.9) print $0
        }' "$blast_file" > "$filtered_file"
    fi
done

echo "Done."






#!/usr/bin/env bash

# Loop through all subdirectories
for dir in */ ; do
    # Look for the filtered .blast.tsv file in the subdirectory
    filtered_file=$(find "$dir" -maxdepth 1 -type f -name "*_filtered.blast.tsv")

    # If such a file exists, process it
    if [[ -n "$filtered_file" ]]; then
        # Create output best-hit filename
        besthit_file="${filtered_file%.blast.tsv}_besthit.blast.tsv"

        echo "Selecting best hits for: $filtered_file"
        echo "Output:                  $besthit_file"

        # Select best hit per query (highest bitscore, column 12)
        awk '
        {
            qid = $1; bits = $12
            if (!(qid in best) || bits > best[qid]) {
                best[qid] = bits
                line[qid] = $0
            }
        }
        END {
            for (i in line) print line[i]
        }' "$filtered_file" > "$besthit_file"
    fi
done

echo "Done."






#!/usr/bin/env bash

# Loop through all subdirectories
for dir in */ ; do
    # Look for the best-hit .blast.tsv file in the subdirectory
    besthit_file=$(find "$dir" -maxdepth 1 -type f -name "*_filtered_besthit.blast.tsv")

    # If present, deconcatenate
    if [[ -n "$besthit_file" ]]; then
        echo "Deconcatenating rMAGs for: $besthit_file"

        # Remove any old r???.blast.tsv files to avoid mixing runs
        find "$dir" -maxdepth 1 -type f -name "r???.blast.tsv" -delete

        # Split by first 3 characters of sseqid (column 2)
        awk -F"\t" -v outdir="$dir" '
        {
            mag_id = "r" substr($2,1,3);     # rXYZ
            outfile = outdir "/" mag_id ".blast.tsv";
            print $0 >> outfile;
            close(outfile);
        }' "$besthit_file"
    fi
done

echo "Done."







####then, move files of the same MAG into its own directory and run:
#!/bin/bash
#SBATCH -J tad80
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 10
#SBATCH --mem=100G
#SBATCH -t8:00:00
#SBATCH -q inferno
#SBATCH --array=1-190%45
#SBATCH -o logs/tad80_%A_%a.out


#### Download metaG ####
eval "$(conda shell.bash hook)"
conda activate roth_env



PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/tad80_attempt/output/Dereplicated/filter1/dereplicated_genomes/blast_out/Chattahoochee"

FASTA_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/tad80_attempt/output/Dereplicated/filter1/dereplicated_genomes"

# Get the sample ID for this array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_ids2.txt)

# Construct paths
B_DIR="${PARENT_DIR}/${SAMPLE}"
FASTA="${FASTA_DIR}/${SAMPLE}.fasta"
OUT_DIR="${PARENT_DIR}/${SAMPLE}_results_attempt2"

# Ensure output dir exists
mkdir -p "$OUT_DIR"

# Check input BLAST directory
if [[ ! -d "$B_DIR" ]]; then
  echo "❌ Error: Input BLAST directory $B_DIR does not exist."
  exit 1
fi

# Check input FASTA file
if [[ ! -f "$FASTA" ]]; then
  echo "❌ Error: Input FASTA file $FASTA not found."
  exit 1
fi

# Log info
echo "✅ Running sample: $SAMPLE"
echo "📂 Input BLAST dir: $B_DIR"
echo "📄 FASTA file: $FASTA"
echo "📁 Output dir: $OUT_DIR"

# Run the script
python3 "/storage/home/hcoda1/9/nmiller304/shared_project/Metagenomic_Population_Tracking/06f_TabBlast_RecPlot_Mini_Auto_v4_edit.py" \
  -f "$FASTA" \
  -b "$B_DIR" \
  -o "$OUT_DIR" \
  -t 0 \
  -d
