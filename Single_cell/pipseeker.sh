#!/bin/bash
#SBATCH --job-name=pipseeker		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=20	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=50G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/KRT_AA_AB_sc/err_out/%x_%j.out		# Standard output log
#SBATCH --error=/scratch/ac05869/KRT_AA_AB_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: UGA Kratom leaf sc analysis
#       Script function: scRNAseq pre-processing and merging at end
#       Input: raw_reads_R1.fastq and raw_read_R2.fastq
#       Output 1: barcode_stats.csv
#       Output 2: barcode_whitelist.txt
#       Output 3: generated_barcode_read_info_table.csv
#		Output 4: barcoded_fastqs/
################################################################################
#INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt #list lib ID's

OUTDIR="/scratch/ac05869/KRT_AA_AB_sc/wd"
#if output directory doesn't exist, create it
if [ ! -d ${OUTDIR} ]
then
    mkdir -p ${OUTDIR}
fi

#set working directory
cd ${OUTDIR}

#INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt #list lib ID's

LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`

mkdir ../${LIB}

/home/ac05869/software/pipseeker-v3.3.0-linux/pipseeker barcode --threads 20 --verbosity 1 --fastq ../rawdata/${LIB} \
--chemistry V --output-path ../${LIB}

cat ../${LIB}/barcoded_*_R1.fastq.gz > ${LIB}_all_barcoded_R1.fastq.gz
cat ../${LIB}/barcoded_*_R2.fastq.gz > ${LIB}_all_barcoded_R2.fastq.gz

#Parameters 
#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt /home/ac05869/NIH/Single_cell/pipseeker.sh

