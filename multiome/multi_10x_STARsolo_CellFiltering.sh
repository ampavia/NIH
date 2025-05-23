#!/bin/bash
#SBATCH --job-name=CellFilter		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=12	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-4				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=50G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=ALL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Kratom multiome
#       Script function: STARsolo CellFiltering (EmptyDrops_CR) algorithm
#       Input: matrix.mtx (renamed from UniqueAndMult-EM.mtx)
#       Output: EM_EmptyDrops_Combined/
################################################################################

LIB1='KRT_AE'
LIB2='KRT_AF'
WD1='/scratch/ac05869/10X_Multiome/KRT_leaf/KRT_AC_AE'
WD2='/scratch/ac05869/10X_Multiome/KRT_leaf/KRT_AD_AF'
ml STAR/2.7.10b-GCC-11.3.0

#KRT_AE
cd ${WD1}/${LIB1}GeneFull/raw/
mkdir data_for_R
cp features.tsv data_for_R
cp barcodes.tsv data_for_R
cp UniqueAndMult-EM.mtx data_for_R
cd data_for_R
mv UniqueAndMult-EM.mtx matrix.mtx

STAR --runThreadN 12 --runMode soloCellFiltering ${WD1}/${LIB1}GeneFull/raw/data_for_R \
${WD1}/${LIB1}GeneFull/raw/data_for_R/EM_EmptyDrops_Combined/ --soloCellFilter EmptyDrops_CR
#The “EM_EmptyDrops_Combined” directory we specify will be created and output will be put there
#Go into the EM_EmptyDrops_Combined directory and gzip everything to be ready for input into Seurat
cd ./EM_EmptyDrops_Combined/
gzip *

#KRT_AF
cd ${WD2}/${LIB2}GeneFull/raw/
mkdir data_for_R
cp features.tsv data_for_R
cp barcodes.tsv data_for_R
cp UniqueAndMult-EM.mtx data_for_R
cd data_for_R
mv UniqueAndMult-EM.mtx matrix.mtx

STAR --runThreadN 12 --runMode soloCellFiltering ${WD2}/${LIB2}GeneFull/raw/data_for_R \
${WD2}/${LIB2}GeneFull/raw/data_for_R/EM_EmptyDrops_Combined/ --soloCellFilter EmptyDrops_CR
#The “EM_EmptyDrops_Combined” directory we specify will be created and output will be put there
#Go into the EM_EmptyDrops_Combined directory and gzip everything to be ready for input into Seurat
cd ./EM_EmptyDrops_Combined/
gzip *
