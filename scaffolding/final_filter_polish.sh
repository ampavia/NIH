#!/bin/bash
#SBATCH --job-name=final_polish_stats	                    # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=32		                        # Number of cores per task - match this to the num_threads used
#SBATCH --mem=300G		                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/gese_final_yahs/err_out/%x_%j.out	# Location of standard output and error log files
#SBATCH --error=/scratch/ac05869/gese_final_yahs/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail 
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

CPU=32
WD='/scratch/ac05869/gese_final_yahs/gese_v2.asm'
CONTIG='/scratch/ac05869/gese_final_yahs/juicer_post/gese_v1.org_filter.50kb.yahs_scaf.JBAT_review.FINAL.fa'
OUT='/scratch/ac05869/gese_final_yahs/gese_v2.asm/no_filter'
OUT2='/scratch/ac05869/gese_final_yahs/gese_v2.asm/50kb'
OUT3='/scratch/ac05869/gese_final_yahs/gese_v2.asm/100kb'
[ -d $WD ] || mkdir -p $WD 
[ -d $OUT ] || mkdir -p $OUT 
[ -d $OUT2 ] || mkdir -p $OUT2
[ -d $OUT3 ] || mkdir -p $OUT3

cd $WD

module load SeqKit/2.9.0
module load BUSCO/5.5.0-foss-2022a

export NUMEXPR_MAX_THREADS=$CPU
cp $CONTIG $OUT/gese_v2.asm.prefilter.fa
seqkit stats -a -o $OUT/prefilter_stats.txt -t dna -j $CPU $CONTIG #finding statistics about your assembly, checking before/after filtering 
busco -i $CONTIG -m genome -l embryophyta_odb10 -c $CPU -o BUSCO_prefilter --out_path $WD #check for completeness of your assembly, before and after filtering 

# use for filtering 100kb
seqkit seq -m 100000 -t dna -j $CPU $CONTIG > $OUT3/gese_v2.asm.100kb_filter.fa
seqkit stats -a -o $OUT3/100kb_filter_stats.txt -t dna -j $CPU $OUT3/gese_v2.asm.100kb_filter.fa
busco -i $OUT3/gese_v2.asm.100kb_filter.fa -m genome -l embryophyta_odb10 -c $CPU -o BUSCO_100kb --out_path $WD

# use for filtering 50kb
seqkit seq -m 50000 -t dna -j $CPU $CONTIG > $OUT2/gese_v2.asm.50kb_filter.fa
seqkit stats -a -o $OUT2/50kb_filter_stats.txt -t dna -j $CPU $OUT2/gese_v2.asm.50kb_filter.fa
busco -i $OUT2/gese_v2.asm.50kb_filter.fa -m genome -l embryophyta_odb10 -c $CPU -o BUSCO_50kb --out_path $WD

#sbatch ~/NIH/scaffolding/final_filter_polish.sh