#!/bin/bash

# Set your base and output directories
base_dir="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"
output_dir="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/all_1030_mags"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each subdirectory in base_dir
for subdir in "$base_dir"/KT*/; do
    # Get the name of the subdirectory (no path)
    subdir_name=$(basename "$subdir")

    # Construct the target path
    derep_dir="${subdir}analysis_nohuman_again/Dereplicated_take3/dereplicated_genomes"

    # Check if the target path exists and is a directory
    if [ -d "$derep_dir" ]; then
        # Loop through all fasta files in that directory
        for fasta_file in "$derep_dir"/*.fasta; do
            # Check if glob matched any file
            [ -e "$fasta_file" ] || continue
            # Get the base name of the fasta file
            fasta_name=$(basename "$fasta_file")
            # Copy with subdir name prefixed
            cp "$fasta_file" "$output_dir/${subdir_name}_${fasta_name}"
        done
    fi
done





#!/bin/bash

# Hardcoded paths
input_dir="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/all_1030_mags"
output_dir="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/all_1030_mags/output"
prefix="Chattahoochee_MAG_"
counter=1

# Make output directory if it doesn't exist
mkdir -p "$output_dir"

# Create mapping file
mapping_file="$output_dir/file_mapping.tsv"
echo -e "Original_File\tNew_File" > "$mapping_file"

# Loop through fasta files (sorted for consistency)
for file in $(ls "$input_dir" | grep -E '\.fa(sta)?$' | sort); do
    new_name=$(printf "%s%03d.fasta" "$prefix" "$counter")
    cp "$input_dir/$file" "$output_dir/$new_name"
    echo -e "$file\t$new_name" >> "$mapping_file"
    counter=$((counter + 1))
done

echo "Copied $(($counter - 1)) files into $output_dir"
echo "Mapping written to $mapping_file"







#!/bin/bash
#SBATCH -J drep_withhuman
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 10
#SBATCH --mem=100G
#SBATCH -t12:00:00
#SBATCH -q inferno


# Load environment
eval "$(conda shell.bash hook)"
conda activate dRep2


INPUT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/tad80_attempt"
CHECKM2DB="/storage/home/hcoda1/9/nmiller304/databases/CheckM2_database/uniref100.KO.1"
OUTPUT_ROOT="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/tad80_attempt/output"

CHECKM_OUT="$OUTPUT_ROOT/checkm2_output"

mkdir -p "$CHECKM_OUT"

echo "📂 Using fasta files from: $INPUT_DIR"


###############################################
# Step 1: Run CheckM2
###############################################
echo "🧪 Running CheckM2..."
checkm2 predict \
    --input "$INPUT_DIR" \
    --output-directory "$CHECKM_OUT" \
    -x .fasta \
    --force \
    --threads 10







   #!/bin/bash
#SBATCH -J drep_withhuman
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 10
#SBATCH --mem=100G
#SBATCH -t12:00:00
#SBATCH -q inferno


# Load environment
eval "$(conda shell.bash hook)"
conda activate dRep2


INPUT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/tad80_attempt/output"
#CHECKM2DB="/storage/home/hcoda1/9/nmiller304/databases/CheckM2_database/uniref100.KO.1"
OUTPUT_ROOT="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/tad80_attempt/output"

CHECKM_OUT="$OUTPUT_ROOT/checkm2_output/filter1_files"
DREP_OUT="$OUTPUT_ROOT/Dereplicated"



###############################################
# Step 2: Run dRep dereplication
###############################################
FIXED_CSV="$CHECKM_OUT/filtered_checkm_table.csv"


echo "🚀 Running dRep dereplication"
dRep dereplicate "$DREP_OUT/filter1" \
    --genomes "/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/tad80_attempt/output/checkm2_output/filter1_files/filtered_MAG_list.txt" \
    --genomeInfo "$FIXED_CSV" \
    --S_algorithm fastANI \
    -sa 0.95 \
    -comp 50 \
    -p 10

echo "✅ Finished processing all fasta files in $INPUT_DIR" 
