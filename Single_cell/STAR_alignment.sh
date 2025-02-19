#!/bin/bash
#SBATCH --job-name=STAR_alignment		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1		# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=24	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=175G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=48:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/KRT_AA_AB_sc/err_out/%x_%j.out		# Standard output log
#SBATCH --error=/scratch/ac05869/KRT_AA_AB_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Star alignment R2 and R1
#       Script function: map reads
#       Input: genome index and barcoded fastq reads 2 and 1
#       Output: BAM
################################################################################
#INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt
LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`

OUTDIR="/scratch/ac05869/KRT_AA_AB_sc/${LIB}/STARsolo"
if [ ! -d ${OUTDIR} ]
then
    mkdir -p ${OUTDIR}
fi

#set working directory
cd ${OUTDIR}

#merge pipseeker barcoded reads for input
cat ../barcoded_fastqs/barcoded_*_R1.fastq.gz > ${LIB}_all_barcoded_R1.fastq.gz
cat ../barcoded_fastqs/barcoded_*_R2.fastq.gz > ${LIB}_all_barcoded_R2.fastq.gz

ml STAR/2.7.10b-GCC-11.3.0

STAR --runThreadN 24 \
--genomeDir /scratch/ac05869/nih/kratom/star_index \
--readFilesIn ../barcoded_fastqs/${LIB}_all_barcoded_R2.fastq.gz \
../barcoded_fastqs/${LIB}_all_barcoded_R1.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix ${LIB}_ \
--outSAMtype BAM SortedByCoordinate \
--alignIntronMax 5000 \
--soloType CB_UMI_Simple \
--soloUMIlen 12 \
--soloBarcodeReadLength 0 \
--soloCBwhitelist ../barcodes/barcode_whitelist.txt \
--soloFeatures GeneFull \
--soloCellFilter EmptyDrops_CR \
--soloMultiMappers EM \

#index bam file
ml purge
ml SAMtools/1.14-GCC-11.2.0
cd ${OUTDIR}
samtools index -@ 24 ${LIB}_Aligned.sortedByCoord.out.bam

#Parameters 
#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt ~/NIH/Single_cell/STAR_alignment.sh
