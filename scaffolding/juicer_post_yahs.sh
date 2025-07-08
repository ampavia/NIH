#!/bin/bash
#SBATCH --job-name=juicer_post_final	                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=12		                        # Number of cores per task - match this to the num_threads used
#SBATCH --mem=100gb			                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/gese_final_yahs/err_out/%x_%j.out	# Location of standard output and error log files
#SBATCH --error=/scratch/ac05869/gese_final_yahs/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

JBAT='/scratch/ac05869/gese_final_yahs/juicebox' #Input directory. already exists
OUT='/scratch/ac05869/gese_final_yahs/juicer_post' #Final assembly output directory
REF='/scratch/ac05869/gese_final_yahs/assembly/gese_v1_organellar_filter.asm.fa' #the ref contig file
PREF='gese_v1.org_filter.asm.yahs' #prefix of the files output by yahs earlier in pipeline
FINAL='gese_v1.org_filter.50kb.yahs_scaf.JBAT_review' #Final assembly ID
[ -d $OUT ] || mkdir -p $OUT


module load YaHS/1.2.2-GCC-11.3.0
module load Juicebox/1.9.9

echo "### Final step: generate final genome assembly file after manual curation with JuiceBox (JBAT)"
## Manual curation consisted of flipping superscaffolds 1,2,5,6,7, and 8 to make short arm appear first. Juicebox output '.review.assembly' was uploaded to "${JBAT}/flipped_chr/${PREF}_JBAT.review.assembly"
juicer post -o ${OUT}/${FINAL} ${JBAT}/flipped_chr/${PREF}_JBAT.review.assembly ${JBAT}/${PREF}_JBAT.liftover.agp ${REF}

## the final output is ${JBAT}/flipped_chr/${PREF}_JBAT.FINAL.agp and ${JBAT}/flipped_chr/${PREF}_JBAT.FINAL.fa
#sbatch ~/NIH/scaffolding/juicer_post_yahs.sh