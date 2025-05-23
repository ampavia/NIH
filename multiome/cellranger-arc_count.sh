#!/bin/bash
#SBATCH --job-name=CR_ARC_count
#SBATCH --partition=highmem_p		# Partition name (batch, heighten_p, or gpu_p), _required_
#SBATCH --ntasks=1 		# Run job in single task or in paralelle, _required_
#SBATCH --cpus-per-task=16		# CPU cores per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=300G			# How much memory per node, _required_
#SBATCH --time=7-00:00:00		# Time Limit hrs:min:sec or day-hrs:min:sec 2-12:00:00 is 2.5 d, _required_
#SBATCH --export=NONE		# Don't export submit node variables to compute node
#SBATCH --output=%x_%j.out	# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu # Send an email when job is done or dead
#SBATCH --mail-type=ALL	# Mail events (BEGIN, END, FAIL, ALL)

LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`
WD='/scratch/ac05869/10X_Multiome/KRT_leaf'
cd $WD

ml CellRanger-ARC/2.0.2-QVupdate

cellranger-arc count --id=${LIB}_CR-arc_count \
 --reference=./mitr_v1 \
 --libraries=./${LIB}/${LIB}.csv \
 --localcores=16 \
 --localmem=300
 
 #sbatch --array 1-2 --export=INFILE=/scratch/ac05869/10X_Multiome/KRT_leaf/multiome_ID.txt ~/NIH/multiome/cellranger-arc_count.sh