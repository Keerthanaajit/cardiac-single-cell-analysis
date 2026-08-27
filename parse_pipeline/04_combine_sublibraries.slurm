#!/bin/bash
#SBATCH -J splitpipe_combine
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH -o splitpipe_combine.out
#SBATCH -e splitpipe_combine.err

module load miniforge
conda activate spipe

# Set this to the location of the project directory on your system
PROJECT_DIR="/path/to/my_project"

cd "${PROJECT_DIR}"

split-pipe \
  --mode comb \
  --sublibraries \
  "${PROJECT_DIR}/parse/sub1" \
  "${PROJECT_DIR}/parse/sub2" \
  --output_dir "${PROJECT_DIR}/parse/combined_final"
