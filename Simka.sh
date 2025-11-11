##first need to create a txt file where you have Sample_name: /path/to/sample

#!/bin/bash
#SBATCH -Jsimka
#SBATCH -A gts-ktk3-coda20
#SBATCH -N1 --ntasks-per-node=24
#SBATCH --mem=100G
#SBATCH -t2:00:00
#SBATCH -q inferno


#### Download metaG ####
#eval "$(conda shell.bash hook)"
#conda activate simka

/storage/home/hcoda1/9/nmiller304/scratch/simka/build/bin/simka -in ./simka_files.txt -out ./simka_out -max-reads 0 -complex-dist -max-memory 360000 -nb-cores 24 -out-tmp /storage/home/hcoda1/9/nmiller304/scratch/simka_tmpdir2
