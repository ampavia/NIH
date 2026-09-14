#!/bin/bash
#SBATCH --job-name=rename_fastq
#SBATCH --partition=batch		# Partition name (batch, heighten_p, or gpu_p), _required_
#SBATCH --ntasks=1 		# Run job in single task or in paralelle, _required_
#SBATCH --cpus-per-task=16		# CPU cores per task
#SBTACH --array=1-4				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=50G			# How much memory per node, _required_
#SBATCH --time=7-00:00:00		# Time Limit hrs:min:sec or day-hrs:min:sec 2-12:00:00 is 2.5 d, _required_
#SBATCH --export=NONE		# Don't export submit node variables to compute node
#SBATCH --output=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.out	# Standard output log
#SBATCH --error=/scratch/ac05869/camptotheca_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu # Send an email when job is done or dead
#SBATCH --mail-type=ALL	# Mail events (BEGIN, END, FAIL, ALL)

LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`
WD='/scratch/ac05869/camptotheca_sc/LEAF/raw'
cd $WD
mkdir ../renamed/

cp ${LIB}*I1.fastq.gz ../renamed/${LIB}_S1_I1_001.fastq.gz
cp ${LIB}*I2.fastq.gz ../renamed/${LIB}_S1_I2_001.fastq.gz
cp ${LIB}*R1.fastq.gz ../renamed/${LIB}_S1_R1_001.fastq.gz
cp ${LIB}*R2.fastq.gz ../renamed/${LIB}_S1_R2_001.fastq.gz

 #sbatch --array 1-4 --export=INFILE=/scratch/ac05869/camptotheca_sc/leaf_atac_gex_libs.txt /scratch/ac05869/camptotheca_sc/scripts/rename_multi_fastq.sh
