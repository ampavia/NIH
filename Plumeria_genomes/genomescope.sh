#!/bin/bash
#SBATCH --job-name=genomescope                   # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=300G		                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#conda init
#conda activate GenomeScope

source ~/miniconda/bin/activate GenomeScope
OUT='/scratch/ac05869/plumeria_datanote/GScope'
[ -d $OUT ] || mkdir -p $OUT 
cd $OUT

##The input histogram_file (from KMC or jellyfish), output_dir, and k-mer_length are required parameters. 
##The optional parameter -p ploidy sets the ploidy of the model for GenomeScope to use. 
##The optional parameter -l lambda sets the initial guess for the average k-mer coverage of the sequencing. 
##The optional parameter -n 'name_prefix' sets the prefix for the output files. 
##The optional parameter -m max_kmercov specifies the cutoff for excluding high frequence k-mers from the analysis. 
##The output plots and a text file of the inferred genome characteristics will be output to the specified output_dir directory.

/PATH/TO/genomescope.R -i ara_F1_21.hist -o output -k 21
~/genomescope2.0/genomescope.R -i $OUT/FastK/PLU_AD_k21.hist -o $OUT -n PLU_AD_21mer_GS -k 21
~/genomescope2.0/genomescope.R -i $OUT/FastK/FRA_AM_k21.hist -o $OUT -n FRA_AM_21mer_GS -k 21

#sbatch ~/NIH/genomescope.sh