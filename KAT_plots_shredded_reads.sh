#!/bin/bash
#SBATCH --job-name=shredded_reads_plot                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=128G		                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error/scratch/ac05869/err_out/=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

WD='/scratch/ac05869/gese_final_yahs/kat_31mer_illumina/seqkit_reads'
ASM='/scratch/ac05869/gese_final_yahs/gese_v2.asm/100kb/gese_v2.asm.fa'
OUT='/scratch/ac05869/gese_final_yahs/kat_31mer_illumina/seqkit_reads/plots'
IN='/scratch/ac05869/gese_final_yahs/hifi_reads/seqkit_sliding_out'
[ -d $WD ] || mkdir -p $WD 
[ -d $OUT ] || mkdir -p $OUT 

ml KAT/2.4.2

kat comp -t 16 -o $WD/gese_v2.asm -m 31 -h -v ${IN}/GEL_AO_hifi_reads_150bp.fastq ${ASM}

kat plot spectra-cn -o ${OUT}/gese_v2.asm_cn-plot -x 600 -y 8000000 $WD/gese_v2.asm-main.mx

