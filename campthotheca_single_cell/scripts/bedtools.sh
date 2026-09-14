#!/bin/bash
#SBATCH --job-name=merge_peaks
#SBATCH --partition=batch		# Partition name (batch, heighten_p, or gpu_p), _required_
#SBATCH --ntasks=1 		# Run job in single task or in paralelle, _required_
#SBATCH --cpus-per-task=4		# CPU cores per task
#SBATCH --mem=50G			# How much memory per node, _required_
#SBATCH --time=168:00:00		# Time Limit hrs:min:sec or day-hrs:min:sec 2-12:00:00 is 2.5 d, _required_
#SBATCH --export=NONE		# Don't export submit node variables to compute node
#SBATCH --output=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.out	# Standard output log
#SBATCH --error=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu # Send an email when job is done or dead
#SBATCH --mail-type=ALL	# Mail events (BEGIN, END, FAIL, ALL)

###After mapping paired end sc-atac reads to genome, merge replicate atac peaks

WD='/scratch/ac05869/camptotheca_sc/CR-ARC/'
OUT='/scratch/ac05869/camptotheca_sc/CR-ARC/merged_peaks'
GFF='/scratch/ac05869/camptotheca_sc/genomes/caac_v4.working_models.agat.repr.gff3'

L1BED='CAA_BG_BF/outs/atac_peaks.bed'
L2BED='CAA_BG_BF/outs/atac_peaks.bed'

[ -d $OUT ] || mkdir -p $OUT 
cd $WD

ml BEDTools/2.31.1-GCC-13.3.0
ml AGAT/1.4.2-GCC-13.3.0

cat $L1BED $L2BED | sort -k1,1 -k2,2n > $OUT/CAA_leaf.atac_peaks.sorted.bed

bedtools merge -i $OUT/CAA_leaf.atac_peaks.sorted.bed > $OUT/CAA_leaf.atac_peaks.sorted.merged.bed


#make gene bed file from gff3
agat_convert_sp_gff2bed.pl --gff $GFF -o $OUT/caac_v4.working_models.agat.repr.bed
cat $OUT/caac_v4.working_models.agat.repr.bed | sort -k1,1 -k2,2n > $OUT/caac_v4.working_models.agat.repr.sorted.bed

#sbatch /scratch/ac05869/camptotheca_sc/scripts/bedtools.sh

#these atac bed files get fed into Seurat
