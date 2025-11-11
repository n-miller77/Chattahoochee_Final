
#Competative Read Mapping 

for f in *.{fa,fna,fas,fasta}; do
  [ -e "$f" ] || continue
  tmp="${f}.tmp"
  sed -E 's/^>([^[:space:]]+).*/>\1/' "$f" > "$tmp" && mv "$tmp" "$f"
done


for f in *.fasta; do
  [ -e "$f" ] || continue
  tag=$(basename "$f" .fasta)
  tag="${tag: -3}"
  tmp="${f}.tmp"
  awk -v t="$tag" '/^>/{print ">" t substr($0, 2); next}1' "$f" > "$tmp" && mv "$tmp" "$f"
done



bowtie2-build all_mags.fasta all_mags 

bowtie2 -x all_mags \
  -1 /storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/KTK04_2_DAM/KTK04_2_DAM_R1_clean.fastq \
  -2 /storage/home/hcoda1/9/nmiller304/shared_project/Chattahooche-samples-work/concat_samples/KTK04_2_DAM/KTK04_2_DAM_R2_clean.fastq \
  -S mapped_mags_DAM.sam \
  -k 1 \
  --no-unal \
  --threads 8 


samtools view -bS mapped_mags_DAM.sam > mapped_mags_DAM.bam
samtools sort -o mapped_mags_DAM_sorted.bam mapped_mags_DAM.bam
samtools index mapped__mags_DAM_sorted.bam
samtools flagstat mapped_sorted.bam

samtools view -c -F 4 -f 64 mapped__mags_DAM_sorted.bam




rpe build -d cc_2_down.db -gf all_chatt_mags.txt --mag -r /storage/home/hcoda1/9/nmiller304/shared_project/tad80_calc/KTK04_2_CC_Down2/KTK04_2_CC_Down2.sorted.sam

rpe plot -d cc_2_down.db --mag_file mags_for_plot.txt --width 3000 -o cc_2_down_plots/













# Total reads in the metagenome (from BAM)
total_reads=$(samtools view -c mapped_mags_DAM.bam)

# Number of mapped reads
mapped_reads=$(samtools view -c -F 4 mapped_mags_DAM.bam)

# Fraction mapped
fraction=$(echo "scale=4; $mapped_reads/$total_reads" | bc)
echo "Fraction of metagenome reads mapping to MAG: $fraction"
