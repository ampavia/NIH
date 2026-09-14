#!/bin/bash
#SBATCH --job-name=caac_v4_gex_pair		# Job name
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=16	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=128G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=48:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=BEGIN,END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Single Cell - Pair the Cleaned R2 and uncleaned R1 reads back together, so they match
#       Script function: Pair reads
#       Input: raw_R1.fastq + cleaned_R2.fastq
#       Output: r1_paired.fastq + r2_paired.fastq
################################################################################
LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`
OUT="/scratch/ac05869/camptotheca_sc/LEAF/paired"
if [ ! -d ${OUT} ]
then
    mkdir -p ${OUT}
fi
cd ${OUT}

ml SeqKit/2.9.0

seqkit pair -1 ../raw/${LIB}_R1.fastq.gz -2 ../cleaned/${LIB}_R2.trim.fastq.gz \
-O /scratch/ac05869/camptotheca_sc/LEAF/paired

#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/camptotheca_sc/10x_libs.txt /scratch/ac05869/camptotheca_sc/scripts/caa_multi_10x_pair_cleaned_reads.sh

#Parameters
#Usage:
#  seqkit pair [flags]
#
#Flags:
#  -f, --force            overwrite output directory
#  -h, --help             help for pair
#  -O, --out-dir string   output directory
#  -1, --read1 string     (gzipped) read1 file
#  -2, --read2 string     (gzipped) read2 file
#  -u, --save-unpaired    save unpaired reads if there are

