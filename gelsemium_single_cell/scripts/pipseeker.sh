#!/bin/bash
#SBATCH --job-name=pipseeker		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=16	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-4				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=100G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=24:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.out		# Standard output log
#SBATCH --error=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: plumeria sc analysis
#       Script function: scRNAseq pre-processing and merging at end
#       Input: raw_reads_R1.fastq and raw_read_R2.fastq (can do multiple flowcells, each prefixed by ${LIB})
#       Output 1: barcode_stats.csv
#       Output 2: barcode_whitelist.txt
#       Output 3: generated_barcode_read_info_table.csv
#		Output 4: barcoded_fastqs/
################################################################################
#INFILE is list lib ID's
LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`
pipseeker='/home/ac05869/software/pipseeker-v3.3.0-linux/pipseeker'
IN='/scratch/ac05869/gelsemium_sc/JES/raw'
PIP='/scratch/ac05869/gelsemium_sc/JES/pipseeker'
if [ ! -d ${PIP} ]
then
    mkdir -p ${PIP}
fi

#set working directory
cd ${PIP}
mkdir ${LIB}

${pipseeker} barcode \
--threads 16 \
--verbosity 1 \
--fastq ${IN}/${LIB} \
--chemistry V \
--output-path ${LIB}

cat ${LIB}/barcoded_fastqs/barcoded_*_R1.fastq.gz > ${LIB}/barcoded_fastqs/${LIB}_all_barcoded_R1.fastq.gz
cat ${LIB}/barcoded_fastqs/barcoded_*_R2.fastq.gz > ${LIB}/barcoded_fastqs/${LIB}_all_barcoded_R2.fastq.gz

#Parameters 
#sbatch --array 1-4 --export=INFILE=/scratch/ac05869/gelsemium_sc/JES_libraries.txt /scratch/ac05869/gelsemium_sc/scripts/pipseeker.sh

