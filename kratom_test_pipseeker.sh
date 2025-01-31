#!/bin/bash
#SBATCH --job-name=kratom_initial_pipseeker		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=20	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=50G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Kratom leaf pipseeker MIT_AH and MIT_AI
#       Script function: scRNAseq pre-processing
#       Input: raw_reads_R1.fastq
#       Output 1: barcode_stats.csv
#       Output 2: barcode_whitelist.txt
#       Output 3: generated_barcode_read_info_table.csv
#		Output 4: barcoded_fastqs/*
################################################################################

#INFILE=/scratch/ac05869/nih/kratom/rawdata/Kratom_test.txt

LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`

/home/ac05869/software/pipseeker-v3.3.0-linux/pipseeker barcode --threads 20 --verbosity 1 --fastq /scratch/ac05869/nih/kratom/rawdata/${LIB} \
--chemistry v4 --output-path /scratch/ac05869/nih/kratom/test_pipseeker/MIT_${LIB}

#Parameters 

#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/nih/kratom/rawdata/Kratom_test.txt /home/ac05869/nih/kratom/kratom_test_pipseeker.sh

