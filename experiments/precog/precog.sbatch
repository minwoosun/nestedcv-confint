#!/usr/bin/bash
#SBATCH --job-name=precog
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=/scratch/users/minwoos/nestedcv-confint/experiments/precog/precog.log
#SBATCH -c 32
#SBATCH --mail-type=BEGIN,END,FAIL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=minwoos@stanford.edu

pwd; hostname; date

# load the module
module load R/4.0.2

# R script to run
EXPERIMENT=/scratch/users/minwoos/nestedcv-confint/experiments/precog/precog.R

# R script arguments
INPUT_PATH=None
OUTPUT_PATH=/scratch/users/minwoos/ncv_output/
JOB_NAME=precog

Rscript $EXPERIMENT $INPUT_PATH $OUTPUT_PATH $JOB_NAME
