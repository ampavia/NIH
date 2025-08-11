#!/bin/bash
#SBATCH --job-name=MACS2_Cebola
#SBATCH --partition=batch		# Partition name (batch, heighten_p, or gpu_p), _required_
#SBATCH --ntasks=1 		# Run job in single task or in paralelle, _required_
#SBATCH --cpus-per-task=6		# CPU cores per task
#SBATCH --mem=128G			# How much memory per node, _required_
#SBATCH --time=168:00:00		# Time Limit hrs:min:sec or day-hrs:min:sec 2-12:00:00 is 2.5 d, _required_
#SBATCH --export=NONE		# Don't export submit node variables to compute node
#SBATCH --output=/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis/err_out/%x_%j.out	# Standard output log
#SBATCH --error=/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu # Send an email when job is done or dead
#SBATCH --mail-type=ALL	# Mail events (BEGIN, END, FAIL, ALL)

###After mapping paired end sc-atac reads to genome, call peaks using CebolaLab github code

cd /scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis/data
OUT="/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis/MACS2"
[ -d $OUT ] || mkdir -p $OUT 

ml MACS2/2.2.7.1-foss-2021b

#Call peaks
for bam in *_sorted.marked.filtered.shifted.bam
do
macs2 callpeak -f BEDPE \
--nomodel \
--shift -37 \
--extsize 73 \
-g 1433922938 \
-B --broad \
--keep-dup all \
--cutoff-analysis \
-n $bam \
-t ${bam/.bam/}.bed \
--outdir ${OUT}/
done

#sbatch ~/NIH/multiome/call_peaks_MACS2.sh
