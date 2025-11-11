#!/bin/bash
#SBATCH -J abundance
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 4
#SBATCH --mem=64GB
#SBATCH -t 24:00:00
#SBATCH -q inferno
#SBATCH --array=0-45%20
#SBATCH -o logs/abundance_nohuman_%A_%a.out


#!/bin/bash

# Load environment
eval "$(conda shell.bash hook)"
conda activate MetaBat

# Set parent directory
PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"

# Generate sample_ids2.txt if it doesn't exist
if [[ ! -f sample_ids2.txt ]]; then
    echo "Generating sample_ids2.txt from $PARENT_DIR"
    find "$PARENT_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort > sample_ids2.txt
fi

# Get sample ID and directory for this array index
SAMPLE_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" sample_ids2.txt)
SAMPLE_DIR="${PARENT_DIR}/${SAMPLE_ID}"

echo "Processing: $SAMPLE_DIR"

if [[ ! -d "$SAMPLE_DIR" ]]; then
    echo "ERROR: Sample directory not found: $SAMPLE_DIR"
    exit 1
fi


# Find the normal2_spades_output* directory
SPADES_DIR=$(find "$SAMPLE_DIR" -maxdepth 1 -type d -name "nohuman_spades_output_*" | head -n 1)

if [[ -z "$SPADES_DIR" ]]; then
    echo "  No normal2_spades_output_* directory found in $SAMPLE_DIR, skipping."
    exit 0
fi

# Check for contigs.fasta
CONTIGS="${SPADES_DIR}/contigs.fasta"
if [[ ! -f "$CONTIGS" ]]; then
    echo "  contigs.fasta not found in $SPADES_DIR, skipping."
    exit 0
fi

# Find paired-end reads
R1=$(find "$SAMPLE_DIR" -type f -name "*_R1*clean_removed.fastq.gz" | head -n 1)
R2=$(find "$SAMPLE_DIR" -type f -name "*_R2*clean_removed.fastq.gz" | head -n 1)

if [[ -z "$R1" || -z "$R2" ]]; then
    echo "  Missing read files in $SAMPLE_DIR, skipping."
    exit 0
fi

# Prepare analysis directory
ANALYSIS_DIR="${SAMPLE_DIR}/analysis_nohuman_again"
mkdir -p "$ANALYSIS_DIR"
cd "$ANALYSIS_DIR" || exit 0

# Run analysis tools
echo "  Building Bowtie2 index..."
bowtie2-build "$CONTIGS" bowtie_contigs2.fasta

echo "  Aligning reads to contigs..."
bowtie2 -p 16 -x bowtie_contigs2.fasta -1 "$R1" -2 "$R2" -S contigs_vs_reads2.sam

echo "  Converting and sorting SAM/BAM..."
samtools view -S --threads 10 -b contigs_vs_reads2.sam > contigs_vs_reads2.bam
samtools sort contigs_vs_reads2.bam -o contigs_vs_reads2.sorted.bam

echo "  Generating depth file..."
jgi_summarize_bam_contig_depths --outputDepth coverage_for_metabat2.txt contigs_vs_reads2.sorted.bam
cut -f1,3 coverage_for_metabat2.txt | sed '1d' > coverage_for_maxbin2.txt

echo "Done with $SAMPLE_DIR"
