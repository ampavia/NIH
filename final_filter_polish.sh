#!/bin/bash
#SBATCH --job-name=gese_v2_polish	                    # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=32		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=300G		                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/gese_final_yahs/gese_v2.asm/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/gese_final_yahs/gese_v2.asm/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

CPU=32
WD='/scratch/ac05869/gese_final_yahs/gese_v2.asm'
CONTIG='/scratch/ac05869/gese_final_yahs/gese_v2.asm/gese_v2.asm_JBAT.FINAL.fa'
OUT='/scratch/ac05869/gese_final_yahs/gese_v2.asm/100kb'
M100='/scratch/ac05869/gelsemium_yahs/gese_v2.asm/100kb/gese_v2.asm.fa'
[ -d $OUT ] || mkdir -p $OUT 

module load SeqKit/2.9.0
module load BUSCO/5.8.3-foss-2023a

export NUMEXPR_MAX_THREADS=$CPU

#seqkit stats -a -o $OUT/prefilter_stats.txt -t dna -j $CPU $CONTIG #finding statistics about your assembly, checking before/after filtering 
#busco -i $CONTIG -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_prefilter --out_path $WD #check for completeness of your assembly, before and after filtering 

# use for filtering 100kb
#seqkit seq -m 100000 -t dna -j $CPU $CONTIG > $OUT/gese_v2.asm.fa
seqkit stats -a -o $OUT/100kb_filter_stats.txt -t dna -j $CPU $M100
busco -i $M100 -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_100kb --out_path $WD

#sbatch ~/NIH/final_filter_polish.sh