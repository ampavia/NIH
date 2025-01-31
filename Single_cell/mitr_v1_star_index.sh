#!/bin/bash
#SBATCH --job-name=mitr_star_index		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=24	 	# CPU core count per task, by default 1 CPU core per task
#SBATCH --mem=50G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Kratom new genome index
#       Script function: create genome index
#       Input: genome.fa and annotation.gtf
#       Output 1: genome index
################################################################################

#mkdir /scratch/ac05869/nih/kratom/star_index_mitr_v1
#interact -c 4 --mem=10G
#ml AGAT/1.1.0
#cd /scratch/ac05869/nih/kratom/mitr_v1_genome
#agat_convert_sp_gff2gtf.pl --gff mitr_v1.working_models.gff3 -o mitr_v1.working_models.gtf
#exit

cd /scratch/ac05869/nih/kratom

ml purge

ml STAR/2.7.10b-GCC-11.3.0

STAR --runThreadN 24 \
--runMode genomeGenerate mode \
--genomeDir /scratch/ac05869/nih/kratom/star_index_mitr_v1 \
--genomeFastaFiles /scratch/ac05869/nih/kratom/mitr_v1_genome/mitr_v1.asm.fa \
--sjdbGTFfile /scratch/ac05869/nih/kratom/mitr_v1_genome/mitr_v1.working_models.gtf

#Parameters 

#sbatch /home/ac05869/nih/kratom/mitr_v1_star_index.sh
