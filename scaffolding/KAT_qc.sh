#!/bin/bash
#SBATCH --job-name=kat_shredded                    # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBTACH --array=1-3			# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=500G		                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

ASM=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`
WD='/scratch/ac05869/gese_final_yahs/kat_31mer_illumina'
LAB=$(basename "${ASM}")
OUT='/scratch/ac05869/gese_final_yahs/kat_31mer_illumina/plots'
IN='/scratch/ac05869/gese_final_yahs/hifi_reads'
[ -d $WD ] || mkdir -p $WD 
[ -d $OUT ] || mkdir -p $OUT 
[ -d $IN ] || mkdir -p $IN 


ml KAT/2.4.2
echo ${SLURM_ARRAY_TASK_ID}_${ASM}

#kat comp -t 16 -o $WD/$LAB -m 31 -h -v ${IN}/wgsim_illumina/5M_reads_?.fq ${ASM}
kat comp -t 16 -o $WD/$LAB -m 31 -h -v ${IN}/seqkit_sliding_out/GEL_AO_hifi_reads_150bp.fastq ${ASM}

kat plot spectra-cn -o ${OUT}/${LAB}_cn-plot -x 600 -y 8000000 $WD/${LAB}-main.mx

#sbatch --array 1-3 --export=INFILE=/scratch/ac05869/gese_final_yahs/assembly_list.txt ~/NIH/scaffolding/KAT_qc.sh

