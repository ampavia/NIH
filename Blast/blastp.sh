#!/bin/bash
#SBATCH --job-name=blastn_RED03		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/nih/kratom_pathway/%x_%j.out			# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)
################################################################################
#Project: blast new allwin gene
#       Script function: blastn
#       Query: RED03 marker gene fa
#       Output: top 4 hits
################################################################################
#set output directory and input file variables
QUERY="/scratch/ac05869/nih/kratom_pathway/RED03.fasta"                # replace cbergman in the following line with your myid
CDNA="/scratch/ac05869/nih/kratom/mitr_v1_anno/mitr_v1.working_models.cdna.fa"
FILE=$(basename $CDNA)
OUTDIR="/scratch/ac05869/nih/kratom_pathway/blast/db" 
#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

cd ${OUTDIR}
cp ${CDNA} . 

module load BLAST+/2.14.1-gompi-2023a
makeblastdb -in ${FILE} -parse_seqids -dbtype nucl

blastn -num_threads 4 \
       -query $QUERY \
       -db ${CDNA} \
       -max_target_seqs 10 \
       -outfmt 6 \
       -out ../blastn.${SLURM_JOB_ID}.tsv
       
#sort hits to get top 4
cd ..
export LANG=C; export LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr blastn.${SLURM_JOB_ID}.tsv | sort -u -k1,4 --merge > top_hit.${SLURM_JOB_ID}.tsv

#sbatch blast.sh
