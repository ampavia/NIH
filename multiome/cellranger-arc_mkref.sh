#!/bin/bash
#SBATCH --job-name=CR_ARC_mkref_2
#SBATCH --partition=batch		# Partition name (batch, heighten_p, or gpu_p), _required_
#SBATCH --ntasks=1 		# Run job in single task or in paralelle, _required_
#SBATCH --cpus-per-task=16		# CPU cores per task
#SBATCH --mem=128G			# How much memory per node, _required_
#SBATCH --time=12:00:00		# Time Limit hrs:min:sec or day-hrs:min:sec 2-12:00:00 is 2.5 d, _required_
#SBATCH --export=NONE		# Don't export submit node variables to compute node
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Standard output log
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu # Send an email when job is done or dead
#SBATCH --mail-type=ALL	# Mail events (BEGIN, END, FAIL, ALL)

WD='/scratch/ac05869/10X_Multiome/KRT_leaf/err_out'
GTF='/scratch/ac05869/10X_Multiome/KRT_leaf/genome/mitr_v1.working_models.agat.gtf3'
CONFIG='/scratch/ac05869/10X_Multiome/KRT_leaf/genome/config.txt'
cd $WD
ml CellRanger-ARC/2.0.2

cellranger-arc mkref --config=$CONFIG --nthreads=16 --memgb=128 --ref-version=1