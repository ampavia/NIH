#!/bin/bash
#SBATCH --job-name=star_MIT_AE		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=24	 	# CPU core count per task, by default 1 CPU core per task
#SBATCH --mem=175G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=48:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Kratom star alignment MIT_AE R1 and R2
#       Script function: map reads
#       Input: genome index and barcoded fastq reads 2 and 1
#       Output: BAM
################################################################################

cd /scratch/ac05869/nih/kratom/star_results/
mkdir MIT_AE_redo
cd MIT_AE_redo
ml purge
ml STAR/2.7.10b-GCC-11.3.0

STAR --runThreadN 24 \
--genomeDir /scratch/ac05869/nih/kratom/star_index \
--readFilesCommand zcat \
--readFilesIn /scratch/ac05869/nih/kratom/test_pipseeker/MIT_AE_redo/barcoded_fastqs/MIT_AE_all_barcoded_R2.fastq.gz \
/scratch/ac05869/nih/kratom/test_pipseeker/MIT_AE_redo/barcoded_fastqs/MIT_AE_all_barcoded_R1.fastq.gz \
--outFileNamePrefix MIT_AE_ \
--soloBarcodeReadLength 0 \
--soloUMIlen 12 \
--alignIntronMax 5000 \
--soloCellFilter EmptyDrops_CR \
--soloMultiMappers EM \
--soloFeatures GeneFull \
--soloType CB_UMI_Simple \
--soloCBwhitelist /scratch/ac05869/nih/kratom/test_pipseeker/MIT_AE_redo/barcodes/barcode_whitelist.txt \
--outSAMtype BAM SortedByCoordinate


#Parameters 

#sbatch /home/ac05869/nih/kratom/star_alignment_single_library.sh
