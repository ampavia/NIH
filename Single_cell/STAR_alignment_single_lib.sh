#!/bin/bash
#SBATCH --job-name=STAR_KRT_AA		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1		# aka threads. Each task by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=24	 	# CPU core count per task, by default 1 CPU core per task
#SBATCH --mem=100GB			# Memory per node (30GB); by default using M as unit
#SBATCH --time=48:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/KRT_AA_AB_sc/err_out/%x_%j.out		# Standard output log
#SBATCH --error=/scratch/ac05869/KRT_AA_AB_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Single Cell pipeline for failed array job library KRT_AA
#       Script function: mapping R2 and R1 to genome index
#       Input: genome index and barcoded fastq reads 2 and 1
#       Output: BAM
################################################################################

LIB="KRT_AA"

OUTDIR="/scratch/ac05869/KRT_AA_AB_sc/${LIB}/STARsolo"
if [ ! -d ${OUTDIR} ]
then
    mkdir -p ${OUTDIR}
fi

#set working directory
cd ${OUTDIR}

ml STAR/2.7.10b-GCC-11.3.0

STAR --runThreadN 24 \
--genomeDir /scratch/ac05869/nih/kratom/star_index_mitr_v1/ \
--readFilesIn ../barcoded_fastqs/${LIB}_all_barcoded_R2.fastq.gz \
../barcoded_fastqs/${LIB}_all_barcoded_R1.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix ${LIB}_ \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 17296580376 \
--alignIntronMax 5000 \
--soloType CB_UMI_Simple \
--soloUMIlen 12 \
--soloBarcodeReadLength 0 \
--soloCBwhitelist ../barcodes/barcode_whitelist.txt \
--soloFeatures GeneFull \
--soloCellFilter EmptyDrops_CR \
--soloMultiMappers EM \

#Parameters 
#sbatch ~/NIH/Single_cell/STAR_alignment.sh
