#!/bin/bash
#SBATCH --job-name=star_index		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=16	 	# CPU core count per task, by default 1 CPU core per task
#SBATCH --mem=50G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: single cell processing
#       Script function: create genome index of camptotheca genome assemblies
#       Input: genome.fa and annotation.gtf
#       Output 1: genome index
################################################################################
ml AGAT/1.4.2-GCC-13.3.0
ml STAR/2.7.11b-GCC-13.3.0
IN='/work/crblab/ac05869/genomes/caac_v4_anno/export/'
#mkdir /scratch/ac05869/camptotheca_sc/genomes
cd /scratch/ac05869/camptotheca_sc/genomes

cp ${IN}/caac_v4.asm.fa .
#cp ${IN}/caac_v4.working_models.gff3 . #Need to filter for the longest isoform

#agat_sp_keep_longest_isoform.pl -gff caac_v4.working_models.gff3 -o caac_v4.repr.working_models.gff3
agat_convert_sp_gff2gtf.pl --gff caac_v4.repr.working_models.gff3 -o caac_v4.repr.working_models.agat.gtf

STAR --runThreadN 16 \
--runMode genomeGenerate mode \
--genomeDir . \
--genomeFastaFiles ./caac_v4.asm.fa \
--sjdbGTFfile ./caac_v4.repr.working_models.agat.gtf


#sbatch /scratch/ac05869/camptotheca_sc/scripts/STAR_genome_index.sh