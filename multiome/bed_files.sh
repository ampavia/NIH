#!/bin/bash
#SBATCH --job-name=gene_atac_bed                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=150gb			                            # Total memory for job
#SBATCH --time=48:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


WD='/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis'
OUT='/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis/bedtools'
BED1='../KRT_AC_AE_atac_CR/outs/atac_peaks.bed'
BED2='../KRT_AD_AF_atac_CR/outs/atac_peaks.bed'
GFF='/scratch/ac05869/10X_Multiome/KRT_leaf/genome/mitr_v1.working_models.repr.gff3'
cd $WD
ml BEDTools/2.30.0-GCC-12.2.0

#merge bed files for atac peaks
cat $BED1 $BED2 | sort -k1,1 -k2,2n > $OUT/AC_AD.atac_peaks.sorted.bed
bedtools merge -i $OUT/AC_AD.atac_peaks.sorted.bed > $OUT/atac_peaks.sorted.merged.bed

#make bed file from gene gff3
ml AGAT
agat_convert_sp_gff2bed.pl --gff $GFF -o 1_test_agat.bed
awk '$1 ~ /_1$/' 1_test_agat.bed | sort -k1,1 -k2,2n > hap1_mitr_v1.sorted.agat.bed




#sbatch ~/NIH/multiome/bed_files.sh