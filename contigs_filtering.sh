#!/bin/bash
#SBATCH --job-name=asm_qc	                    # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=32		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=300G		                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

NUMEXPR_MAX_THREADS=32
REF='/scratch/ac05869/gelsemium_yahs/gese_v1.asm.fa'
CPU=32
CONTIG='/scratch/ac05869/gelsemium_yahs/filter_contigs'
M10='/scratch/ac05869/gelsemium_yahs/filter_contigs/min10_gese_v1.asm.filtered.fa'
M50='/scratch/ac05869/gelsemium_yahs/filter_contigs/min50_gese_v1.asm.filtered.fa'
[ -d $CONTIG ] || mkdir -p $CONTIG 

module load SeqKit/2.9.0
module load BUSCO/5.8.3-foss-2023a

seqkit stats -a -o $CONTIG/prefilter_stats.txt -t dna -j $CPU $REF #finding statistics about your assembly, checking before/after filtering 
busco -i $REF -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_$(basename "REF") --out_path $CONTIG #check for completeness of your assembly, before and after filtering 

# use for filtering 10kb
seqkit seq -m 10000 -o $M10 -t dna -j $CPU $REF
seqkit stats -a -o $CONTIG/10kb_filter_stats.txt -t dna -j $CPU $M10
busco -i $M10 -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_$(basename "$M10") --out_path $CONTIG

# use for filtering 50kb
seqkit seq -m 50000 -o $M50 -t dna -j $CPU $REF
seqkit stats -a -o $CONTIG/50kb_filter_stats.txt -t dna -j $CPU $M50
busco -i $M50 -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_$(basename "$M50") --out_path $CONTIG
