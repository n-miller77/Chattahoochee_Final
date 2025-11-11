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
