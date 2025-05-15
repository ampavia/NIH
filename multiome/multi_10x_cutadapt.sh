#!/bin/bash
#SBATCH --job-name=KRT_R2_cutadapt		# Job name
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=15	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=50G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=48:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=ALL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Single Cell - Clean Reads by Quality Trimming, Removing Adapters, and PolyA Tails
#       Script function: Run Cutadapt on Read 2 gex libs
#       Input: raw_reads.fastq
#       Output: trimmed_reads.fastq
################################################################################
LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`

WD='/scratch/ac05869/10X_Multiome/KRT_leaf/raw_fastq'
cd $WD
OUT='/scratch/ac05869/10X_Multiome/KRT_leaf/raw_fastq/cutadapt_out'
[ -d $OUT ] || mkdir -p $OUT


ml cutadapt/3.5-GCCcore-11.2.0
cutadapt UMGC.220425.Aviti4_Run_123/${LIB}_R2_001.fastq.gz \
--cores=15 -q 30 -m 30 --trim-n -n 2 \
-g AAGCAGTGGTATCAACGCAGAGTACATGGG \
-a "A{20}" \
-o $OUT/${LIB}_R2_001.trim.fastq.gz

#Parameters
#file type fastq (autodetects)
#-q: trim bases with a quality score less than 30 from the beginning of the read (3' end)
#-m: minimum length of 100nt for each read
#--trim-n: trim all N bases (bases with no call) (done after adapter trimming)
#-n 2: remove up to two adapters
#-g: remove 5' TSO (flagged as Clontech SMARTer) adapter from read 2 AAGCAGTGGTATCAACGCAGAGTACATGGG
#NOTE: I am trimming this way as I am trimming read 2! for 10X data
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview
#-a: "A{100}" trims polyA tails from reads, used as directed by cutadapt manual for PolyA trimming. 3' end
#-o: output file
#-p: second output file for paired end cleaning

#sbatch --array 1-2 --export=INFILE=/scratch/ac05869/10X_Multiome/KRT_leaf/multiome_ID.txt ~/NIH/multiome/multi_10x_cutadapt.sh