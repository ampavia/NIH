#!/bin/bash
#SBATCH --job-name=EmptyDrops		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=12	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=20G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/KRT_AA_AB_sc/err_out/%x_%j.out		# Standard output log
#SBATCH --error=/scratch/ac05869/KRT_AA_AB_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Single Cell pipeline
#       Script function: STARsolo EM CellFiltering and matrix for R
#       Input: UniqueAndMult-EM.mtx
#       Output: EM_EmptyDrops_Combined
################################################################################
#INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt
LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`

OUTDIR="/scratch/ac05869/KRT_AA_AB_sc/${LIB}/data_for_R"
#if output directory doesn't exist, create it
if [ ! -d ${OUTDIR} ]
then
    mkdir -p ${OUTDIR}
fi

cd ${OUTDIR}
cp ../STARsolo/${LIB}_Solo.out/GeneFull/raw/* ${OUTDIR}
mv UniqueAndMult-EM.mtx matrix.mtx #rename so that it is recognized by STARsolo

ml STAR/2.7.10b-GCC-11.3.0
STAR --runMode soloCellFiltering ${OUTDIR} ${OUTDIR}/EM_EmptyDrops_Combined/ --soloCellFilter EmptyDrops_CR

gzip ${OUTDIR}/EM_EmptyDrops_Combined/*

#Parameters 
#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt ~/NIH/Single_cell/STAR_CellFiltering.sh
