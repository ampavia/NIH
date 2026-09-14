#!/bin/bash
#SBATCH --job-name=gex-STAR		# Job name 
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
#Project: Single Cell - Align RNAseq Reads to Genome
#       Script function: Align RNAseq Reads to Genome
#       Input: trimmed_reads.fastq
#       Output: alignments.sam --> alignmnets.sorted.bam
################################################################################

LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`
IN='/scratch/ac05869/camptotheca_sc/LEAF/paired' 
OUT='/scratch/ac05869/camptotheca_sc/LEAF/STAR'
if [ ! -d ${OUT} ]
then
    mkdir -p ${OUT}
fi

mkdir ${OUT}/${LIB}
cd ${OUT}/${LIB} #very important to start run in unique directories so that files do not overwrite other files in the array

ml STAR/2.7.11b-GCC-13.3.0
ml CellRanger-ARC/2.0.2-QVupdate
#cp /apps/eb/CellRanger-ARC/2.0.2-QVupdate/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz ${OUT}
#gunzip ${OUT}/737K-arc-v1.txt.gz

STAR --genomeDir /scratch/ac05869/camptotheca_sc/genomes \
--readFilesIn ${IN}/${LIB}_R2.trim.fastq.gz ${IN}/${LIB}_R1.fastq.gz \
--readFilesCommand zcat \
--runThreadN 24 \
--alignIntronMax 5000 \
--soloBarcodeReadLength 0 \
--soloUMIlen 12 \
--soloCellFilter EmptyDrops_CR \
--soloFeatures GeneFull \
--soloMultiMappers EM \
--soloType CB_UMI_Simple \
--soloCBwhitelist ${OUT}/737K-arc-v1.txt \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 26843545600 \
--soloOutFileNames $LIB


#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/camptotheca_sc/10x_libs.txt /scratch/ac05869/camptotheca_sc/scripts/caa_multi_10x_starsolo.sh
