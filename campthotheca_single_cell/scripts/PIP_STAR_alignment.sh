#!/bin/bash
#SBATCH --job-name=STAR_alignment		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1		# aka threads. Each task by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=24	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=100GB			# Memory per node (30GB); by default using M as unit
#SBATCH --time=48:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.out		# Standard output log
#SBATCH --error=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Single Cell processing for pipseq cDNA libraries
#       Script function: mapping R2 and R1 to genome index
#       Input: genome index and fastq reads 2 and 1 and barcodes from Pipseeker
#       Output: BAM
################################################################################
LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`
PIP='/scratch/ac05869/camptotheca_sc/pipseeker'
OUT='/scratch/ac05869/camptotheca_sc/STEM/STAR'

if [ ! -d ${OUT} ]
then
    mkdir -p ${OUT}
fi

mkdir ${OUT}/${LIB}
cd ${OUT}/${LIB}

ml STAR/2.7.11b-GCC-13.3.0

STAR --genomeDir /scratch/ac05869/camptotheca_sc/genomes \
--readFilesIn ${PIP}/${LIB}/barcoded_fastqs/${LIB}_all_barcoded_R2.fastq.gz \
${PIP}/${LIB}/barcoded_fastqs/${LIB}_all_barcoded_R1.fastq.gz \
--runThreadN 24 \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--alignIntronMax 5000 \
--soloUMIlen 12 \
--soloBarcodeReadLength 0 \
--soloCBwhitelist ${PIP}/${LIB}/barcodes/barcode_whitelist.txt \
--soloCellFilter EmptyDrops_CR \
--soloFeatures GeneFull \
--soloMultiMappers EM \
--soloType CB_UMI_Simple \
--soloOutFileNames ${LIB}

#Parameters 
#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/camptotheca_sc/pip_libs.txt /scratch/ac05869/camptotheca_sc/scripts/PIP_STAR_alignment.sh
