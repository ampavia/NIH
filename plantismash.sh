#!/bin/bash
#SBATCH --job-name=plantismash_coex	                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=12		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=200gb			                            # Total memory for job
#SBATCH --time=48:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out			# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err			# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

WD='/scratch/ac05869'
GFF="/scratch/ac05869/plantismash_mitr.v1_data/mitr_v1.working_models.repr.gff3"
FA='/scratch/ac05869/plantismash_mitr.v1_data/mitr_v1.asm.fa'
COEX='/scratch/ac05869/plantismash_mitr.v1_data/KRT_RNA_integrated_cluster_av_expression_repr_transcript.csv'
OUT='/scratch/ac05869/plantismash_kratom'

[ -d $OUT ] || mkdir -p $OUT 

conda activate plantismash
#python ~/plantismash/run_antismash.py -h

python ~/plantismash/run_antismash.py --verbose --debug --statusfile $OUT/status.txt --limit -1 --taxon plants --outputfolder $OUT --use_phase --coexpress $COEX --gff $GFF $FA 
# Please check the error message. The genome names in the gff3 file may differ from those in the fasta file, causing an error.

#sbatch ~/NIH/plantismash.sh