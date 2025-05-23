#!/bin/bash
#SBATCH --job-name=parse_anno
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=ac05869@uga.edu
#SBATCH --mail-type=ALL

#load the prerequisite modules

eval "$(conda shell.bash hook)"

conda activate genespace

# run genespace
Rscript parse_annotation.R
