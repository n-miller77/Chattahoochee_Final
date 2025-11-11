#!/bin/bash
#SBATCH -J nonpareil
#SBATCH -A gts-ktk3-coda20
#SBATCH -n 8
#SBATCH --mem=10G
#SBATCH -t24:00:00
#SBATCH -q inferno

eval "$(conda shell.bash hook)"
conda activate base

nonpareil -s ${input}_clean_1.fastq.gz -T kmer -f fastq -b ${input%%.fastq.gz} -X 50000 -t 

# For pair-end reads use only one sister pair. That is, use only the forward or reverse reads.



















library(Nonpareil)

samples <- read.table('npo_files6.tsv', sep='\t', header=TRUE, as.is=TRUE)
nps <- Nonpareil.set(samples$File, col=samples$Col, labels=samples$Name)
summary(nps)  # Displays completeness and Nd

library(Nonpareil)

samples <- read.table('npo_files3.tsv', sep='\t', header=TRUE, as.is=TRUE)
nps <- Nonpareil.set(samples$File, labels=samples$Name)
summary(nps)  # Displays completeness and Nd
