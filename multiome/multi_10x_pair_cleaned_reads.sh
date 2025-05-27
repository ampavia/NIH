#!/bin/bash
#SBATCH --job-name=pair_reads		# Job name
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=5	 	# CPU core count per task, by default 1 CPU core per task
#SBATCH --mem=125G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=24:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=ALL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Single Cell - Pair the Cleaned R2 and uncleaned R1 reads back together, so they match
#       Script function: Pair reads
#       Input: raw_R1.fastq + cleaned_R2.fastq
#       Output: r1_paired.fastq + r2_paired.fastq
################################################################################

R2='/scratch/ac05869/10X_Multiome/KRT_leaf/raw_fastq/cutadapt_out'
WD='/scratch/ac05869/10X_Multiome/KRT_leaf/raw_fastq/err_out'
OUT1='/scratch/ac05869/10X_Multiome/KRT_leaf/raw_fastq/AE_gex_paired'
OUT2='/scratch/ac05869/10X_Multiome/KRT_leaf/raw_fastq/AF_gex_paired'
[ -d $OUT1 ] || mkdir -p $OUT1
[ -d $OUT2 ] || mkdir -p $OUT2
ml SeqKit/0.16.1

cd $WD

#AE
seqkit pair -1 ./UMGC.220425.Aviti4_Run_123/KRT_AE_S2_R1_001.fastq.gz -2 $R2/KRT_AE_S2_R2_001.trim.fastq.gz \
-O $OUT1

#AF
seqkit pair -1 ./UMGC.220425.Aviti4_Run_123/KRT_AF_S3_R1_001.fastq.gz -2 $R2/KRT_AF_S3_R2_001.trim.fastq.gz \
-O $OUT2



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

#sbatch ~/NIH/multiome/multi_10x_pair_cleaned_reads.sh