#!/bin/bash
#SBATCH --job-name=shredHIFI            				 # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=128		                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

HIFI='/scratch/ac05869/gese_final_yahs/hifi_reads/GEL_AO_hifi_reads.fastq'
OUT='/scratch/ac05869/gese_final_yahs/hifi_reads/wgsim_illumina'
OUT2='/scratch/ac05869/gese_final_yahs/hifi_reads/seqkit_sliding_out'
[ -d $OUT ] || mkdir -p $OUT 
[ -d $OUT2 ] || mkdir -p $OUT2 


#shred hifi reads into illumina length reads
~/wgsim/wgsim -e 0 -N 5000000 -1 150 -2 150 -r 0 -R 0 -X 0 \
$HIFI \
$OUT/5M_reads_1.fq $OUT/5M_reads_2.fq

#seqkit may work better
module load SeqKit/2.9.0
seqkit sliding $HIFI -j 16 -W 150 -s 150 -g -o $OUT2/GEL_AO_hifi_reads_150bp.fastq

#sbatch ~/NIH/scaffolding/shredHIFIreads_asm_QC.sh
