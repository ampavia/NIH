#!/bin/bash
#SBATCH --job-name=samtools		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=12			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=1	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-4				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=20G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=12:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Kratom MIT_AD, MIT_AE, MIT_AH, and MIT_AI
#       Script function: Index BAM output from star
#       Input: {LIB}_Aligned.sortedByCoord.out.bam
#       Output: bam.bai
# interact --mem 20G --ntasks 12 --time 2:00:00
## srun --pty  --cpus-per-task=1 --ntasks=12 --nodes=1 --partition=inter_p --time=2:00:00 --mem=20G /bin/bash -l
################################################################################

#INFILE=/scratch/ac05869/nih/kratom/rawdata/kratom_leaf_and_stem.txt

LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`

ml SAMtools

cd /scratch/ac05869/nih/kratom/star_results/${LIB}/
samtools index -@ 12 ${LIB}_Aligned.sortedByCoord.out.bam

#sbatch --array 1-4 --export=INFILE=/scratch/ac05869/nih/kratom/rawdata/kratom_leaf_and_stem.txt /home/ac05869/nih/kratom/samtools_bam_to_bam.bai.sh