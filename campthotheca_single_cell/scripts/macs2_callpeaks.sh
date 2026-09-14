#!/bin/bash
#SBATCH --job-name=MACS2                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=150gb			                            # Total memory for job
#SBATCH --time=48:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.out	
#SBATCH --error=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

IN='/scratch/ac05869/camptotheca_sc/CR-ARC/'
OUT="/scratch/ac05869/camptotheca_sc/MACS2"
[ -d $OUT ] || mkdir -p $OUT 
cd $IN

ml MACS2/2.2.9.1-foss-2023a


#Call peaks
for celltype in *_CAA_leaf_fragments.bed
do
macs2 callpeak -f BED \
--nomodel \
--shift -37 \
--extsize 73 \
-g 429188713 \
--keep-dup all \
--broad \
--cutoff-analysis \
-n ${celltype/.bed/} \
-t ${celltype} \
--outdir ${OUT}
done

module purge
ml BEDTools/2.31.1-GCC-13.3.0

cat $OUT/*_CAA_leaf_fragments_peaks.broadPeak | sort -k1,1 -k2,2n > $OUT/CAA_leaf.celltype_broadPeak.sorted.bed
bedtools merge -i $OUT/CAA_leaf.celltype_broadPeak.sorted.bed > $OUT/CAA_leaf.celltype_broadPeak.sorted.merged.bed

#sbatch /scratch/ac05869/camptotheca_sc/scripts/macs2_callpeaks.sh


#notes on --shift and --extsize options
#https://www.biostars.org/p/209592/
#For certain nucleosome-seq data, we need to pileup the centers of nucleosomes using a half-nucleosome size
#for wavelet analysis (e.g. NPS algorithm). Since the DNA wrapped on nucleosome is about 147bps, this option
#can be used: '--nomodel --shift 37 --extsize 73'.
