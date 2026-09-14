#!/bin/bash
#SBATCH --job-name=mkgtf_mkref
#SBATCH --partition=batch		# Partition name (batch, heighten_p, or gpu_p), _required_
#SBATCH --ntasks=1 		# Run job in single task or in paralelle, _required_
#SBATCH --cpus-per-task=16		# CPU cores per task
#SBATCH --mem=64G			# How much memory per node, _required_
#SBATCH --time=12:00:00		# Time Limit hrs:min:sec or day-hrs:min:sec 2-12:00:00 is 2.5 d, _required_
#SBATCH --export=NONE		# Don't export submit node variables to compute node
#SBATCH --output=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.out	# Standard output log
#SBATCH --error=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu # Send an email when job is done or dead
#SBATCH --mail-type=ALL	# Mail events (BEGIN, END, FAIL, ALL)

WD='/scratch/ac05869/gelsemium_sc/CR-ARC'
GFF='/work/crblab/ac05869/genomes/gese_mpi_v2/gese_mpi_v2.working_models.gff3' #non repr models
GFF2='/scratch/ac05869/gelsemium_sc/genomes/gese_mpi_v2.working_models.agat.repr.gff3'
GTF='/scratch/ac05869/gelsemium_sc/genomes/gese_mpi_v2.working_models.agat.repr.gtf' #non repr
GTF2='gese_mpi_v2.working_models.agat.repr.CR_ARC.gtf'
CONFIG='/scratch/ac05869/gelsemium_sc/config.txt'
[ -d $WD ] || mkdir -p $WD 

cd $WD

ml CellRanger-ARC/2.0.2
ml AGAT/1.4.2-GCC-13.3.0

agat_sp_keep_longest_isoform.pl -f $GFF -o $GFF2 #keep only one transcript per gene before making reference.

module purge
ml gffread/0.12.7-GCCcore-12.3.0 #cellranger only likes this version of gtf, not agat's version.
gffread $GFF2 -T -o $GTF

module purge
ml CellRanger-ARC/2.0.2
cellranger-arc mkgtf $GTF $GTF2 --attribute=gene_biotype:protein_coding #
cellranger-arc mkref --config=$CONFIG --nthreads=16 --memgb=64 --ref-version=1 #make sure config file has correct paths and filenames for gtf and the genome fasta


#sbatch /scratch/ac05869/gelsemium_sc/scripts/ATAC_1.sh
