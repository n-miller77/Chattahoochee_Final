###############running sourceapp
#!/bin/bash
#SBATCH -J sourceapp2_loop
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t12:00:00
#SBATCH -q inferno
#SBATCH --array=1-45%20
#SBATCH -o logs/sourceapp_%A_%a.out


#### Download metaG ####
eval "$(conda shell.bash hook)"
conda activate sourceapp



RAW_DIR=/storage/home/hcoda1/9/nmiller304/scratch/Chattahooche-samples-working
OUT_DIR=/storage/home/hcoda1/9/nmiller304/scratch/Sourceapp_Output_Combined
FILELIST=/storage/home/hcoda1/9/nmiller304/scratch/Chatt_scripts/sample_ids2.txt

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILELIST)


READ1=$(find "${RAW_DIR}/${SAMPLE}" -name "*R1*clean.fastq.gz")
READ2=$(find "${RAW_DIR}/${SAMPLE}" -name "*R2*clean.fastq.gz")

INPUT="${READ1},${READ2}"
OUTPUT="${OUT_DIR}/${SAMPLE}"

python /storage/home/hcoda1/9/nmiller304/SourceApp/pipelines/sourceapp.py \
  --use-geq \
  --drop-env \
  --aggregate-human \
  -i "$INPUT" \
  -o "$OUTPUT" \
  -d /storage/home/hcoda1/9/nmiller304/p-ktk3-0/SourceApp_database \
  -t 8

echo "SourceApp finished for $SAMPLE"
exit 0
