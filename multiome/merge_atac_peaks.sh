#!/bin/bash
#SBATCH --job-name=atac_peaks_merge	                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=150gb			                            # Total memory for job
#SBATCH --time=48:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


WD='/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_data/'
OUT='/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_data/bedtools'
BED1='../KRT_AC_AE_atac_CR/outs/atac_peaks.bed'
BED2='../KRT_AD_AF_atac_CR/outs/atac_peaks.bed'
cd $WD
ml BEDTools/2.30.0-GCC-12.2.0

cat $BED1 $BED2 | sort -k1,1 -k2,2n > $OUT/AC_AD.atac_peaks.sorted.bed

bedtools merge -i $OUT/AC_AD.atac_peaks.sorted.bed > $OUT/atac_peaks.sorted.merged.bed

#sbatch ~/NIH/multiome/merge_atac_peaks.sh


