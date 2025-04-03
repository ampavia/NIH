#!/bin/bash
#SBATCH --job-name=bam_file		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=6		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set output directory and input file variables
OUTDIR="/work/crblab/ac05869/yahs_gelsemium"                  #
#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd ${OUTDIR}

#module
module load BWA/0.7.17-GCCcore-11.3.0 #will create indeces for reference and map reads against reference 
module load SAMtools/1.16.1-GCC-11.3.0 #will sort the bam file that is being created, then index the bam file

##- map paired end reads to reference genome, output as sort BAM with index file
bwa index ${OUTDIR}/gese_v1.asm.fa
bwa mem -t 6 ${OUTDIR}/gese_v1.asm.fa ${OUTDIR}/gel-an_1438200/gel-an_1438201_S3HiC_R1.fastq.gz ${OUTDIR}/gel-an_1438200/gel-an_1438201_S3HiC_R2.fastq.gz | samtools sort -O BAM --threads 6 - > ${OUTDIR}/gel-an_1438201_S3HiC.sorted.bam
samtools index -@ 6 ${OUTDIR}/gel-an_1438201_S3HiC.sorted.bam
