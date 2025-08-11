#!/bin/bash
#SBATCH --job-name=sc-atac_bed                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=150gb			                            # Total memory for job
#SBATCH --time=48:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_qc/err_out/%x_%j.out	
#SBATCH --error=/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_qc/err_out/%x_%j.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

###
###After calling peaks with MACS2, take the bed file of peaks and merge reps
###

WD='/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_qc'
OUT='/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_qc/merged_peaks'
BED1='./MACS2/KRT_AC_all/KRT_AC_all.sorted.bam_peaks.broadPeak'
BED2='./MACS2/KRT_AD_all/KRT_AD_all.sorted.bam_peaks.broadPeak'
GFF='/scratch/ac05869/10X_Multiome/KRT_leaf/genome/mitr_v1.working_models.repr.gff3'
[ -d $OUT ] || mkdir -p $OUT 
cd $WD
ml BEDTools/2.30.0-GCC-12.2.0

#merge bed files for atac peaks
cat $BED1 $BED2 | sort -k1,1 -k2,2n > $OUT/AC_AD.atac_peaks.sorted.bed
bedtools merge -i $OUT/AC_AD.atac_peaks.sorted.bed > $OUT/atac_peaks.sorted.merged.bed

#make bed file from gene gff3
ml AGAT
agat_convert_sp_gff2bed.pl --gff $GFF -o ./data/mitr_v1.agat.bed
sort -k1,1 -k2,2n ./data/mitr_v1.agat.bed > ./data/mitr_v1.sorted.agat.bed

#sbatch ~/NIH/multiome/bed_files.sh