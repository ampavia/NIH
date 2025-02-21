#!/bin/bash
#SBATCH --job-name=BLAST-pop_markers		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/nih/pop_markers/blast/log.%j			# Location of standard output and error log files
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)
################################################################################
#Project: blast poplar stem and leaf marker peptide against mitr_v1
#       Script function: blastp
#       Query: poplar_leaf_stem_ath_pep_query.fa
#       Output: blastp_hits.
################################################################################

#set output directory and input file variables
OUTDIR="/scratch/ac05869/nih/pop_markers/blast"
QUERY="/scratch/ac05869/nih/pop_markers/poplar_leaf_stem_ath_pep_query.fa"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load blast module and make database
ml BLAST+/2.11.0-gompi-2019b
mkdir /scratch/ac05869/nih/pop_markers/blast/db
cd /scratch/ac05869/nih/pop_markers/blast/db
cp /scratch/ac05869/nih/kratom/mitr_v1_anno/mitr_v1.hc_gene_models.repr.pep.fa .
makeblastdb -in mitr_v1.hc_gene_models.repr.pep.fa -parse_seqids -dbtype prot

# run blast against local copy of NCBI nucleotide database
blastp -num_threads 4 \
       -query $QUERY \
       -db mitr_v1.hc_gene_models.repr.pep.fa \
       -max_target_seqs 10 \
       -outfmt 6 \
       -out $OUTDIR/blastp_hits.${SLURM_JOB_ID}.tsv

#sort hits
module purge
cd $OUTDIR
export LANG=C; export LC_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr blastp_hits.${SLURM_JOB_ID}.tsv | sort -u -k1,1 --merge > bestHits_pop_mitsp_v1_markers.${SLURM_JOB_ID}.tsv
       
#sbatch /scratch/ac05869/nih/pop_markers/blast.sh
