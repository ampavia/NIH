#!/bin/bash
#SBATCH --job-name=CR_ARC_count
#SBATCH --partition=batch		# Partition name (batch, heighten_p, or gpu_p), _required_
#SBATCH --ntasks=1 		# Run job in single task or in paralelle, _required_
#SBATCH --cpus-per-task=16		# CPU cores per task
#SBATCH --mem=128G			# How much memory per node, _required_
#SBATCH --time=168:00:00		# Time Limit hrs:min:sec or day-hrs:min:sec 2-12:00:00 is 2.5 d, _required_
#SBATCH --export=NONE		# Don't export submit node variables to compute node
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Standard output log
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu # Send an email when job is done or dead
#SBATCH --mail-type=ALL	# Mail events (BEGIN, END, FAIL, ALL)

WD='/scratch/ac05869/10X_Multiome/KRT_leaf'
cd $WD
ml CellRanger-ARC/2.0.2

cellranger-arc count --id=KRT_AC_AE \
 --reference=mitr_v1/ \
 --libraries=./KRT1.csv \
 --localcores=16 \
 --localmem=128G
 
 cellranger-arc count --id=KRT_AD_AF \
 --reference=mitr_v1/ \
 --libraries=./KRT2.csv \
 --localcores=16 \
 --localmem=128G
 
