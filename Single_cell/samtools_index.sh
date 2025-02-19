#!/bin/bash
#SBATCH --job-name=samtools		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=12			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=1	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-4				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=20G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/KRT_AA_AB_sc/err_out/%x_%j.out		# Standard output log
#SBATCH --error=/scratch/ac05869/KRT_AA_AB_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Single cell pipeline
#       Script function: Index bam files from STARsolo
#       Input: ${LIB}_Aligned.sortedByCoord.out.bam
#       Output: bam.bai
# interact --mem 20G --ntasks 12 --time 2:00:00
################################################################################
#INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt
LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`

#set working directory
cd ${OUTDIR}


#index bam file
ml purge
ml SAMtools/1.14-GCC-11.2.0
samtools index -@ 12 ${LIB}_Aligned.sortedByCoord.out.bam

#sbatch --array 1-4 --export=INFILE=/scratch/ac05869/KRT_AA_AB_sc/leaf_libs.txt ~/NIH/Single_cell/samtools_index.sh