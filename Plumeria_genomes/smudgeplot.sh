#!/bin/bash
#SBATCH --job-name=smudgeplot                   # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=300G		                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#conda init
#conda activate FastK

source ~/miniconda/bin/activate FastK

HIFI="/scratch/ac05869/plumeria_datanote/plumeria_hifi_fastq"
WD='/scratch/ac05869/plumeria_datanote'
ASM='/scratch/ac05869/gs_plcu_plal_cro/genomeRepo'
ASM2='/scratch/ac05869/nih/Plcu_v2/plcu_v2.asm.fa'
ASM3='/scratch/ac05869/nih/Plal_v1/plal_v1.asm.fa'
OUT='/scratch/ac05869/plumeria_datanote/FastK'
SMDG='/scratch/ac05869/plumeria_datanote/FastK/smdg'
[ -d $OUT ] || mkdir -p $OUT
[ -d $SMDG ] || mkdir -p $SMDG
cd $WD
FastK -v -t1 -k21 -M300 -T16 -N${OUT}/FastK_Table_FRA_AM $HIFI/FRA_AM_150bp.fastq
Histex -G ${OUT}/FastK_Table_FRA_AM > $OUT/FRA_AM_k21.hist

FastK -v -t1 -k21 -M300 -T16 -N${OUT}/FastK_Table_PLU_AD ${HIFI}/PLU_AD_150bp.fastq
Histex -G ${OUT}/FastK_Table_PLU_AD > $OUT/PLU_AD_k21.hist

# Find all k-mer pairs in the dataset using hetmer module
smudgeplot.py hetmers -L 12 -t 16 -o $SMDG/kmerpairs_FRA_AM --verbose ${OUT}/FastK_Table_FRA_AM
smudgeplot.py hetmers -L 12 -t 16 -o $SMDG/kmerpairs_PLU_AD --verbose ${OUT}/FastK_Table_PLU_AD

# this now generated `_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot.py all -o $SMDG/kmerpairs_FRA_AM_trial_run -t P.rubra_21mer $SMDG/kmerpairs_FRA_AM_text.smu
smudgeplot.py all -o $SMDG/kmerpairs_PLU_AD_trial_run -t P.cubensis_21mer $SMDG/kmerpairs_PLU_AD_text.smu

# check that bunch files are generated (3 pdfs; some summary tables and logs)
#ls data/Scer/trial_run_*

#sbatch ~/NIH/smudgeplot.sh