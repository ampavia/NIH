#!/bin/bash
#SBATCH --job-name=merge_barcoded_fastqs		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=1	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=50G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Kratom - merge pipseeker barcoded reads for leaf
#       Script function: concatenate barcoded files together into one to run STAR
#       Input: barcoded_R*.fastq
#       Output 1: all_barcoded_R*.fastq.gz
################################################################################

#INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt

OUTDIR="/scratch/ac05869/KRT_AA_AB_sc/wd"
cd ${OUTDIR}
LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`

cat ../${LIB}/barcoded_fastqs/barcoded_*_R1.fastq.gz > ${LIB}_all_barcoded_R1.fastq.gz
cat ../${LIB}/barcoded_fastqs/barcoded_*_R2.fastq.gz > ${LIB}_all_barcoded_R2.fastq.gz

#Parameters 
#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt /home/ac05869/NIH/Single_cell/merge_barcoded_fastqs_pipseeker.sh
