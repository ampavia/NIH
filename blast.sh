#!/bin/bash
#SBATCH --job-name=BLAST-leaf		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/nih/leaf_markers/log.%j			# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)
################################################################################
#Project: blast C.roseus leaf marker genes against mitr_v1
#       Script function: blastn
#       Query: c.ros_v3 marker gene fa
#       Output: mitr_v1 marker genes tsv
################################################################################
#set output directory and input file variables
OUTDIR="/scratch/ac05869/nih/leaf_markers/blast"                  # replace cbergman in the following line with your myid
QUERY="/home/ac05869/nih/leaf_markers_query.fa"                # replace cbergman in the following line with your myid

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load blast module
module load BLAST+/2.14.1-gompi-2023a
#move into db
cd /scratch/ac05869/nih/leaf_markers/db
makeblastdb -in mitr_v1.working_models.cdna.fa -parse_seqids -dbtype nucl

# run blast against local copy of NCBI nucleotide database
# note this is one command split over multiple lines ending in \ for improved readability
blastn -num_threads 4 \
       -query $QUERY \
       -db /scratch/ac05869/nih/leaf_markers/db/mitr_v1.working_models.cdna.fa \
       -max_target_seqs 2 \
       -outfmt 6 \
       -out $OUTDIR/query.fa.blastn.${SLURM_JOB_ID}.tsv
       
#sbatch 