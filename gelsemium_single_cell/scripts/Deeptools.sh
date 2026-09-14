#!/bin/bash
#SBATCH --job-name=metaplot_matrix
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=48gb
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.out	# Standard output log
#SBATCH --error=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu
#SBATCH --mail-type=ALL

LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`

BED='/scratch/ac05869/gelsemium_sc/CR-ARC/merged_peaks'
WD='/scratch/ac05869/gelsemium_sc/CR-ARC/deeptools_out'
[ -d $WD ] || mkdir -p $WD 
cd $WD

# deepTools bam to bigwig
ml deepTools/3.5.5-gfbf-2023a

# bamCoverage --bam ../${LIB}/outs/atac_possorted_bam.bam \
# --outFileName ${LIB}.bigwig \
# --outFileFormat bigwig \
# --ignoreDuplicates \
# --normalizeUsing CPM

# Make matrix at genes
# computeMatrix scale-regions -S ${LIB}.bigwig \
#                             -R ${BED}/gese_mpi_v2.working_models.agat.repr.sorted.bed \
#                             --beforeRegionStartLength 3000 \
#                             --regionBodyLength 5000 \
#                             --afterRegionStartLength 3000 \
#                             -out ${LIB}.gene.tab.gz \
#                             --skipZeros
#                             
# # Matrix at peaks
# computeMatrix reference-point -S ${LIB}.bigwig \
#                             -R ${BED}/GEL_leaf.atac_peaks.sorted.merged.bed \
#                             --referencePoint center \
#                             -b 2000 \
#                             -a 2000 \
#                             -out ${LIB}.peak.tab.gz \
#                             --skipZeros

# Heatmaps
plotHeatmap -m ${LIB}.peak.tab.gz  \
      --colorMap YlGnBu \
      --heatmapHeight 4 \
      --heatmapWidth 6 \
      --missingDataColor "white" \
      --legendLocation none \
      --xAxisLabel "" \
      --refPointLabel Peak \
      --regionsLabel Peaks \
      -out ${LIB}_peaks.svg

plotHeatmap -m ${LIB}.gene.tab.gz \
            --colorMap YlGnBu \
            --heatmapHeight 4 \
            --heatmapWidth 6 \
            --missingDataColor "white" \
            --legendLocation none \
            -out ${LIB}_genes.svg

#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/gelsemium_sc/LEAF_multi_libraries.txt /scratch/ac05869/gelsemium_sc/scripts/Deeptools.sh