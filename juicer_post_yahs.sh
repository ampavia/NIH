#!/bin/bash
#SBATCH --job-name=juicer_post	                    # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=32		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=250gb			                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/gese_final_yahs/gese_v2.asm/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/gese_final_yahs/gese_v2.asm/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

JBAT='/scratch/ac05869/gese_final_yahs/juicebox'
OUT='/scratch/ac05869/gese_final_yahs/gese_v2.asm'
PREF='gese_v1.org_filter.asm.yahs' #prefix of the files
FINAL='gese_v2.asm'
REF='/scratch/ac05869/gese_final_yahs/assembly/gese_v1_organellar_filter.asm.fa' #the contig file

module load YaHS/1.2.2-GCC-11.3.0
module load Juicebox/1.9.9

echo "### Final step: generate final genome assembly file after manual curation with JuiceBox (JBAT)"
## I do not want to keep the manual curation of the assembly, so I am using "${JBAT}/${PREF}_JBAT.assembly" instead of "${JBAT}/${PREF}_JBAT.review.assembly"
## the final output is ${JBAT}/${PREF}_JBAT.FINAL.agp and ${JBAT}/${PREF}_JBAT.FINAL.fa
juicer post -o ${OUT}/${FINAL}_JBAT ${JBAT}/${PREF}_JBAT.assembly ${JBAT}/${PREF}_JBAT.liftover.agp ${REF}


#sbatch ~/NIH/juicer_post_yahs.sh