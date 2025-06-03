#!/bin/bash
#SBATCH --job-name=final_polish_stats	                    # Job name
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
OUT2='/scratch/ac05869/gese_final_yahs/gese_v2.asm/50kb'
[ -d $OUT ] || mkdir -p $OUT 
[ -d $OUT2 ] || mkdir -p $OUT2
cd $WD

module load SeqKit/2.9.0
module load BUSCO/5.8.3-foss-2023a

export NUMEXPR_MAX_THREADS=$CPU

seqkit stats -a -o $OUT/prefilter_stats.txt -t dna -j $CPU $CONTIG #finding statistics about your assembly, checking before/after filtering 
busco -i $CONTIG -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_prefilter --out_path $WD #check for completeness of your assembly, before and after filtering 

# use for filtering 100kb
seqkit seq -m 100000 -t dna -j $CPU $CONTIG > $OUT/gese_v2.asm.fa
seqkit stats -a -o $OUT/100kb_filter_stats.txt -t dna -j $CPU $OUT/gese_v2.asm.fa
busco -i $OUT/gese_v2.asm.fa -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_100kb --out_path $WD

# use for filtering 50kb
seqkit seq -m 50000 -t dna -j $CPU $CONTIG > $OUT2/gese_v2.asm.50kb.fa
seqkit stats -a -o $OUT2/50kb_filter_stats.txt -t dna -j $CPU $OUT2/gese_v2.asm.50kb.fa
busco -i $OUT2/gese_v2.asm.50kb.fa -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_50kb --out_path $WD

#sbatch ~/NIH/scaffolding/final_filter_polish.sh