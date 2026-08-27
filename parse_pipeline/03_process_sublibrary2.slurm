#!/bin/bash
#SBATCH -J splitpipe_sub2
#SBATCH --cpus-per-task=24
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH -o splitpipe_sub2.out
#SBATCH -e splitpipe_sub2.err

module load miniforge
conda activate spipe

# Set this to the location of the project directory on your system
PROJECT_DIR="/path/to/my_project"

cd "${PROJECT_DIR}"

split-pipe \
  --mode all \
  --kit WT_mini \
  --chemistry v3 \
  --genome_dir "${PROJECT_DIR}/genomes/GRCm38_mCherry/GRCm38_mCherry_ref" \
  --fq1 \
  "${PROJECT_DIR}/expdata/DFA_161224_2/DFA_161224_2_1.fq.gz" \
  "${PROJECT_DIR}/expdata/resequencing/X204SC24103974-Z01-F001/01.RawData/DFA_180325_2_R1.fq.gz" \
  --fq2 \
  "${PROJECT_DIR}/expdata/DFA_161224_2/DFA_161224_2_2.fq.gz" \
  "${PROJECT_DIR}/expdata/resequencing/X204SC24103974-Z01-F001/01.RawData/DFA_180325_2_R2.fq.gz" \
  --output_dir "${PROJECT_DIR}/parse/sub2" \
  --sample 288WT_280WT A1 \
  --sample 295WT_296WT A2 \
  --sample 290WT_289WT A3 \
  --sample 286WT_291WT A4 \
  --sample 295BT A5 \
  --sample 294BT A6 \
  --sample 278BT A7 \
  --sample 296BT A8 \
  --sample 286BT A9 \
  --sample 286BT2 A10 \
  --sample 289BT A11 \
  --sample 291BT A12 \
  --nthreads 24
