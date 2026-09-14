#!/bin/bash
#SBATCH --job-name=JES_EM		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=12	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-4				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=20G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: single cell processing
#       Script function: STARsolo CellFiltering (EmptyDrops_CR) algorithm
#       Input: matrix.mtx (renamed from UniqueAndMult-EM.mtx)
#       Output: EM_EmptyDrops_Combined/
################################################################################

LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`

IN='/scratch/ac05869/gelsemium_sc/JES/STAR'

cd ${IN}/${LIB}/${LIB}GeneFull/raw
mkdir data_for_R
cp features.tsv data_for_R/
cp barcodes.tsv data_for_R/
cp UniqueAndMult-EM.mtx data_for_R/
cd data_for_R/
mv UniqueAndMult-EM.mtx matrix.mtx

ml STAR/2.7.11b-GCC-13.3.0

STAR --runThreadN 12 \
--runMode soloCellFiltering \
../data_for_R \
./EM_EmptyDrops_Combined/ \
--soloCellFilter EmptyDrops_CR \

cd ./EM_EmptyDrops_Combined/
gzip *

#sbatch --array 1-4 --export=INFILE=/scratch/ac05869/gelsemium_sc/JES_libraries.txt /scratch/ac05869/gelsemium_sc/scripts/JES_STAR_EM.sh
