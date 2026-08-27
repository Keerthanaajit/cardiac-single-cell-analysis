#!/bin/bash
#SBATCH -J mkref
#SBATCH -n 16
#SBATCH -t 10:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -o mkref.out

source ~/.bashrc
source activate spipe

# Set this to the location of the project directory on your system
PROJECT_DIR="/path/to/my_project"

split-pipe \
  --mode mkref \
  --genome_name GRCm38_mCherry \
  --nthreads 16 \
  --fasta "${PROJECT_DIR}/genomes/GRCm38_mCherry.fa" \
  --genes "${PROJECT_DIR}/genomes/Mus_musculus.GRCm38.99.gtf" \
  --output_dir GRCm38_mCherry_ref
