#!/bin/bash
#SBATCH --job-name=blastn_kratom		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files 
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=ALL                         # Mail events (BEGIN, END, FAIL, ALL)
################################################################################
#Project: blast kratom MIA genes
################################################################################
#set output directory and input file variables
QUERY="/scratch/ac05869/nih/kratom_pathway/kratom_pathway_query.na.fa"                # replace cbergman in the following line with your myid
CDNA="/scratch/ac05869/nih/kratom/mitr_v1_anno/mitr_v1.working_models.cdna.fa"
OUTDIR="/scratch/ac05869/nih/kratom_pathway/blast" 
#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

mkdir ${OUTDIR}/na_db
cd ${OUTDIR}/na_db
cp ${CDNA} . 

module load BLAST+/2.14.1-gompi-2023a
makeblastdb -in mitr_v1.working_models.cdna.fa -parse_seqids -dbtype nucl

blastn -num_threads 4 \
       -query ${QUERY} \
       -db mitr_v1.working_models.cdna.fa \
       -max_target_seqs 10 \
       -outfmt 6 \
       -out ../blastn.${SLURM_JOB_ID}.tsv
       
#sbatch /scratch/ac05869/nih/kratom_pathway/mitr_v1_query.blastn.sh
