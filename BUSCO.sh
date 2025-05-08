#!/bin/bash
#SBATCH --job-name=BUSCO                    # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=32		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=250gb			                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

STATS='/scratch/ac05869/gese_final_yahs/yahs/final_stats'
BUSCO='/scratch/ac05869/gese_final_yahs/yahs/busco_downloads'
SCAF='/scratch/ac05869/gese_final_yahs/juicebox/gese_v1.org_filter.asm.yahs_JBAT.FINAL.fa'
CPU=32
LIN='embryophyta_odb10'
cd $STATS

module load SeqKit/2.9.0
module load BUSCO/5.4.7-foss-2022a

#>&2 echo"### Step 1: check stats"
#seqkit stats -a -o $STATS/FINAL_stats.txt -t dna -j $CPU $SCAF #finding statistics about your assembly, checking before/after filtering 
>&2 echo"### Step 2: BUSCO embryophyta_odb10"
busco -i $SCAF -m genome -l $LIN -c $CPU -o gese_v1_scaf --out_path $STATS# --download_path $BUSCO
