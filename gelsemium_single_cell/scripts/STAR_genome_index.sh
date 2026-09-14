#!/bin/bash
#SBATCH --job-name=star_index		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=16	 	# CPU core count per task, by default 1 CPU core per task
#SBATCH --mem=50G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: single cell processing
#       Script function: create genome index of gelsemium genome assemblies
#       Input: genome.fa and annotation.gtf
#       Output 1: genome index
################################################################################
ml AGAT/1.4.2-GCC-13.3.0
ml STAR/2.7.11b-GCC-13.3.0
IN_UGA='/work/crblab/ac05869/genomes/gese_uga_v1/'
IN_MPI='/work/crblab/ac05869/genomes/gese_mpi_v2/'
cd /scratch/ac05869/gelsemium_sc/genomes

## Start with UGA genome assembly

cp /work/crblab/ac05869/genomes/gese_uga_v1/gese_uga_v1.asm.fa ./gese_uga_v1

agat_convert_sp_gff2gtf.pl --gff ${IN_UGA}/gese_uga_v1.working_models.repr.gff3 -o ./gese_uga_v1.working_models.repr.agat.gtf

STAR --runThreadN 16 \
--runMode genomeGenerate mode \
--genomeDir ./gese_uga_v1 \
--genomeFastaFiles gese_uga_v1/gese_uga_v1.asm.fa \
--sjdbGTFfile ./gese_uga_v1.working_models.repr.agat.gtf

## MPI genome assembly
mkdir gese_mpi_v2
cp /work/crblab/ac05869/genomes/gese_mpi_v2/gese_mpi_v2.asm.fa ./gese_mpi_v2

agat_convert_sp_gff2gtf.pl --gff ${IN_MPI}/gese_mpi_v2.working_models.repr.gff3 -o gese_mpi_v2.working_models.repr.agat.gtf

STAR --runThreadN 16 \
--runMode genomeGenerate mode \
--genomeDir ./gese_mpi_v2 \
--genomeFastaFiles gese_mpi_v2/gese_mpi_v2.asm.fa \
--sjdbGTFfile ./gese_mpi_v2.working_models.repr.agat.gtf

#sbatch /scratch/ac05869/gelsemium_sc/scripts/STAR_genome_index.sh