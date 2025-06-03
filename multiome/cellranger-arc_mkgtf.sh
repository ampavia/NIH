#!/bin/bash
#SBATCH --job-name=CR_ARC_mkgtf
#SBATCH --partition=batch		# Partition name (batch, heighten_p, or gpu_p), _required_
#SBATCH --ntasks=1 		# Run job in single task or in paralelle, _required_
#SBATCH --cpus-per-task=16		# CPU cores per task
#SBATCH --mem=128G			# How much memory per node, _required_
#SBATCH --time=168:00:00		# Time Limit hrs:min:sec or day-hrs:min:sec 2-12:00:00 is 2.5 d, _required_
#SBATCH --export=NONE		# Don't export submit node variables to compute node
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Standard output log
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu # Send an email when job is done or dead
#SBATCH --mail-type=ALL	# Mail events (BEGIN, END, FAIL, ALL)

WD='/scratch/ac05869/10X_Multiome/KRT_leaf/genome'
GTF='/scratch/ac05869/10X_Multiome/KRT_leaf/genome/hap1_mitr_v1.working_models.agat.gtf3'
GFF='/scratch/ac05869/10X_Multiome/KRT_leaf/genome/hap1_mitr_v1.working_models.gff3'
cd $WD
ml CellRanger-ARC/2.0.2
ml AGAT/1.1.0

agat_convert_sp_gff2gtf.pl --gff $GFF -o $GTF

cellranger-arc mkgtf $GTF hap1_mitr_v1.working_models.agat.filtered.gtf3 --attribute=gene_biotype:protein_coding