#!/bin/bash
#SBATCH --job-name=minimap	                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=12		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=50G			                            # Total memory for job
#SBATCH --time=48:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

WD="/scratch/ac05869/gese_final_yahs/minimap"                  #
#if output directory doesn't exist, create it
if [ ! -d $WD ]
then
    mkdir -p $WD
fi
cd ${WD}

module load minimap2/2.28-GCCcore-12.3.0

minimap2 -x asm5 /scratch/ac05869/nih/Cro_v3_genome/cro_v3.asm.fa /scratch/ac05869/gese_final_yahs/juicebox/gese_v1.org_filter.asm.yahs_JBAT.FINAL.fa \
-secondary=no > ${WD}/cro_v3_vs_gese_v1_scaffolded.paf

#sbatch ~/NIH/minimap.sh