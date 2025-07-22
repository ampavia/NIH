#!/bin/bash
#SBATCH --job-name=Cro_blastp		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files 
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)
################################################################################
#Project: blast catharanthus MIA peptide sequences with mitr_v1
################################################################################

#set output directory and input file variables
OUTDIR="/scratch/ac05869/nih/kratom_pathway/blast"
QUERY="/scratch/ac05869/nih/kratom_pathway/Cro_v3_query.pep.fa"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

# load blast module and make database
ml BLAST+/2.11.0-gompi-2019b
mkdir /scratch/ac05869/nih/kratom_pathway/blast/pep_db
cd /scratch/ac05869/nih/kratom_pathway/blast/pep_db
cp /scratch/ac05869/nih/kratom/mitr_v1_anno/mitr_v1.hc_gene_models.repr.pep.fa .
makeblastdb -in mitr_v1.hc_gene_models.repr.pep.fa -parse_seqids -dbtype prot

# run blast against local copy of NCBI nucleotide database
blastp -num_threads 4 \
       -query $QUERY \
       -db mitr_v1.hc_gene_models.repr.pep.fa \
       -max_target_seqs 10 \
       -outfmt 6 \
       -out $OUTDIR/Mitsp.v1_MIA_genes_${SLURM_JOB_ID}.tsv

       
#sbatch /scratch/ac05869/nih/kratom_pathway/Cro_v3_query.blastp.sh
