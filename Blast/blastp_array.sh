#!/bin/bash
#SBATCH --job-name=MEP_blast		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --array=1-13
#SBATCH --time=12:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/MEP_evolution/err_out/%x_%j.out	# Location of standard output and error log files 
#SBATCH --error=/scratch/ac05869/MEP_evolution/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)
################################################################################
#Project: blast catharanthus MEP peptide sequences against all other species, then blast Arabidopsis MEP pathway
################################################################################
SPECIES=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`
#set output directory and input file variables
OUTDIR="/scratch/ac05869/MEP_evolution/Blast/out"
OUTDIR2="/scratch/ac05869/MEP_evolution/Blast/top_hit"
Q1="/scratch/ac05869/MEP_evolution/Blast/Catharanthus_query.pep.fa"
Q2="/scratch/ac05869/MEP_evolution/Blast/Arabidopsis_query.pep.fa"
Query1=$(basename "$Q1" .pep.fa)
Query2=$(basename "$Q2" .pep.fa)


#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

if [ ! -d $OUTDIR2 ]
then
    mkdir -p $OUTDIR2
fi

# load blast module and make database
ml BLAST+/2.14.1-gompi-2023a

cd /scratch/ac05869/MEP_evolution/Blast
makeblastdb -in ${SPECIES}.fa -parse_seqids -dbtype prot

# run blast against local copy of NCBI nucleotide database
blastp -num_threads 4 \
       -query $Q1 \
       -db ${SPECIES}.fa \
       -max_target_seqs 50 \
       -outfmt 6 \
       -out $OUTDIR/${SPECIES}_${Query1}.tsv
       
blastp -num_threads 4 \
       -query $Q2 \
       -db ${SPECIES}.fa \
       -max_target_seqs 50 \
       -outfmt 6 \
       -out $OUTDIR/${SPECIES}_${Query2}.tsv

#sort post blast: Make sure the file is sorted based on query and best hits (here bitscore > evalue > perc identity):
export LC_ALL=C LC_LANG=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr ${OUTDIR}/${SPECIES}_${Query1}.tsv > ${OUTDIR}/${SPECIES}_${Query1}_sorted.tsv
export LC_ALL=C LC_LANG=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr ${OUTDIR}/${SPECIES}_${Query2}.tsv > ${OUTDIR}/${SPECIES}_${Query2}_sorted.tsv


#get top queries - top 4 according to bitscore. Change -m value if needing different number
for next in $(cut -f1 ${OUTDIR}/${SPECIES}_${Query1}_sorted.tsv | sort -u); do grep -w -m 4 "$next" ${OUTDIR}/${SPECIES}_${Query1}_sorted.tsv; done > ${OUTDIR2}/${SPECIES}_${Query1}_sorted_top.tsv
for next in $(cut -f1 ${OUTDIR}/${SPECIES}_${Query2}_sorted.tsv | sort -u); do grep -w -m 4 "$next" ${OUTDIR}/${SPECIES}_${Query2}_sorted.tsv; done > ${OUTDIR2}/${SPECIES}_${Query2}_sorted_top.tsv

       
# sbatch --export=INFILE=/scratch/ac05869/MEP_evolution/species.txt /scratch/ac05869/MEP_evolution/run_BLASTp.sh

