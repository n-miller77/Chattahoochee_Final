#!/bin/bash
#SBATCH -Jgtdbtk
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 4
#SBATCH --mem=36G
#SBATCH -t6:00:00
#SBATCH -q inferno



eval "$(conda shell.bash hook)"
conda activate bowtiesam



PARENT_DIR="/storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples"
INDEX="/storage/home/hcoda1/9/nmiller304/scratch/T2T_bowtie/T2T_genome"

for SUBDIR in "$PARENT_DIR"/*/; do
  cd "$SUBDIR"


  if [[ -d "human_read_files" ]]; then
    echo "Skipping $SUBDIR: human_read_files directory already exists"
    cd ..
    continue
  fi



  R1=$(ls *_R1*clean.fastq.gz)
  R2=$(ls *_R2*clean.fastq.gz)

  R1_OUT="${R1/_clean.fastq.gz/_clean_removed.fastq.gz}"
  R2_OUT="${R2/_clean.fastq.gz/_clean_removed.fastq.gz}"

  # Skip if output files already exist
  if [[ -f "$R1_OUT" && -f "$R2_OUT" ]]; then
    echo "Skipping $SUBDIR: output files already exist"
    cd ..
    continue
  fi

  bowtie2 -p 8 -x "$INDEX" -1 "$R1" -2 "$R2" -S mapped_and_unmapped.sam
  samtools view -bS mapped_and_unmapped.sam > mapped_and_unmapped.bam
  samtools view -b -f 12 -F 256 mapped_and_unmapped.bam > bothReadsUnmapped.bam
  samtools sort -n -m 5G -@ 2 bothReadsUnmapped.bam -o bothReadsUnmapped_sorted.bam
  samtools fastq -@ 8 bothReadsUnmapped_sorted.bam \
    -1 "$R1_OUT" \
    -2 "$R2_OUT" \
    -0 /dev/null -s /dev/null -n

mkdir human_read_files

mv *bam human_read_files

mv *sam human_read_files

  cd ..
done
