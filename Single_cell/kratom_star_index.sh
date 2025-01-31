#!/bin/bash
#SBATCH --job-name=kratom_star_index		# Job name 
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
#Project: Kratom genome index Julia's genome
#       Script function: create genome index
#       Input: genome.fa and annotation.gtf
#       Output 1: genome index
################################################################################

#mkdir /scratch/ac05869/nih/kratom/star_index

#run interactive job: interact
#ml AGAT/1.1.0
#cd /scratch/ac05869/nih/kratom/kratom_genome
#agat_convert_sp_gff2gtf.pl --gff kratom.working_models.gff3 -o kratom.working_models.agat.gtf

cd /scratch/ac05869/nih/kratom

ml purge

ml STAR/2.7.10b-GCC-11.3.0

STAR --runThreadN 24 \
--runMode genomeGenerate mode \
--genomeDir /scratch/ac05869/nih/kratom/star_index \
--genomeFastaFiles /scratch/ac05869/nih/kratom/kratom_genome/kratom_final_assembly.fasta \
--sjdbGTFfile /scratch/ac05869/nih/kratom/kratom_genome/kratom.working_models.agat.gtf

#Parameters 


#sbatch /home/ac05869/nih/kratom/kratom_star_index.sh
