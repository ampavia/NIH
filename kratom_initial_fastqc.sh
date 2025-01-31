#!/bin/bash
#SBATCH --job-name=kratom_initial_fastqc		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=10	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=50G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=4:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Kratom - QC Reads
#       Script function: Run fastqc on raw reads
#       Input: raw_reads.fastq
#       Output: fastqc_report.html
################################################################################

#INFILE=/scratch/ac05869/nih/kratom/rawdata/Kratom_test.txt

LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 2 | tail -n 1`

ml purge
ml FastQC/0.11.9-Java-11

fastqc -t 10 -o /scratch/ac05869/nih/kratom/test_fastqc \
/scratch/ac05869/nih/kratom/rawdata/${LIB}*.fastq.gz

#Parameters 

#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/nih/kratom/rawdata/Kratom_test.txt /home/ac05869/nih/kratom/rawdata/kratom_initial_fastqc.sh
