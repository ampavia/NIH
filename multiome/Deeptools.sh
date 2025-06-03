#!/bin/bash
#SBATCH --job-name=atac_deeptools
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ac05869/10X_Multiome/KRT_leaf/err_out/%x.%j.out
#SBATCH --error=/scratch/ac05869/10X_Multiome/KRT_leaf/err_out/%x_%j.err
#SBATCH --mail-user=ac05869@uga.edu
#SBATCH --mail-type=ALL

WD='/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis/deeptools_out'
BAM1='/scratch/ac05869/10X_Multiome/KRT_leaf/KRT_AC_AE_atac_CR/outs/atac_possorted_bam.bam'
BAM2='/scratch/ac05869/10X_Multiome/KRT_leaf/KRT_AD_AF_atac_CR/outs/atac_possorted_bam.bam'
# deepTools bam to bigwig
ml deepTools/3.5.2-foss-2022a

bamCoverage --bam $BAM1 \
--outFileName KRT_AC.bigwig \
--outFileFormat bigwig \
--ignoreDuplicates \
--normalizeUsing CPM

bamCoverage --bam $BAM2 \
--outFileName KRT_AD.bigwig \
--outFileFormat bigwig \
--ignoreDuplicates \
--normalizeUsing CPM

# Make matrix at genes
computeMatrix scale-regions -S KRT_AC.bigwig \
                            -R ../../ATAC_analysis/bedtools/hap1_mitr_v1.sorted.agat.bed \
                            --beforeRegionStartLength 3000 \
                            --regionBodyLength 5000 \
                            --afterRegionStartLength 3000 \
                            -out KRT_AC.gene.tab.gz \
                            --skipZeros


computeMatrix scale-regions -S KRT_AD.bigwig \
                            -R ../../ATAC_analysis/bedtools/hap1_mitr_v1.sorted.agat.bed \
                            --beforeRegionStartLength 3000 \
                            --regionBodyLength 5000 \
                            --afterRegionStartLength 3000 \
                            -out KRT_AD.gene.tab.gz \
                            --skipZeros

# Matrix at peaks
computeMatrix reference-point -S KRT_AC.bigwig \
                            -R ../../ATAC_analysis/bedtools/atac_peaks.sorted.merged.bed \
                            --referencePoint center \
                            -b 2000 \
                            -a 2000 \
                            -out KRT_AC.tab.gz \
                            --skipZeros


computeMatrix reference-point -S KRT_AD.bigwig \
                            -R ../../ATAC_analysis/bedtools/atac_peaks.sorted.merged.bed \
                            --referencePoint center \
                            -b 2000 \
                            -a 2000 \
                            -out KRT_AD.tab.gz \
                            --skipZeros
                            
                            
#Heatmaps KRT_AC
#atac peaks
plotHeatmap -m KRT_AC.tab.gz  \
      --colorMap YlGnBu \
      --heatmapHeight 4 \
      --missingDataColor "white" \
      --legendLocation none \
      --xAxisLabel "" \
      --refPointLabel Peak \
      --regionsLabel Peaks \
      -out KRT_AC_peaks.png
#gene intervals     
plotHeatmap -m KRT_AC.gene.tab.gz \
            --colorMap YlGnBu \
            --heatmapHeight 4 \
            --missingDataColor "white" \
            --legendLocation none \
            -out KRT_AC_genes.png
#Heatmaps KRT_AD
#atac peaks     
plotHeatmap -m KRT_AD.tab.gz  \
      --colorMap YlGnBu \
      --heatmapHeight 4 \
      --missingDataColor "white" \
      --legendLocation none \
      --xAxisLabel "" \
      --refPointLabel Peak \
      --regionsLabel Peaks \
      -out KRT_AD_peaks.png
#gene intervals 
plotHeatmap -m KRT_AD.gene.tab.gz \
            --colorMap YlGnBu \
            --heatmapHeight 4 \
            --missingDataColor "white" \
            --legendLocation none \
            -out KRT_AD_genes.png
            
#sbatch ~/NIH/multiome/Deeptools.sh