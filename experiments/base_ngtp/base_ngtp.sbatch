#!/usr/bin/bash
#SBATCH --job-name=basen
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --output=/scratch/users/minwoos/nestedcv-confint/experiments/base_ngtp/base_ngtp.log
#SBATCH -c 20
#SBATCH --mail-type=END,FAIL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=minwoos@stanford.edu

pwd; hostname; date

# load the module
module load R/4.0.2

# R script to run
EXPERIMENT=/scratch/users/minwoos/nestedcv-confint/experiments/base_ngtp/base_ngtp.R

# R script arguments
INPUT_PATH=None
OUTPUT_PATH=/scratch/users/minwoos/ncv_output/
JOB_NAME=base_ngtp

Rscript $EXPERIMENT $INPUT_PATH $OUTPUT_PATH $JOB_NAME
