#!/bin/bash
#SBATCH --job-name=KAT_plumeria                   # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=16		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=300G		                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)


HIFI="/scratch/ac05869/plumeria_datanote/plumeria_hifi_fastq"
WD='/scratch/ac05869/plumeria_datanote'
ASM='/scratch/ac05869/gs_plcu_plal_cro/genomeRepo'
ASM2='/scratch/ac05869/nih/Plcu_v2/plcu_v2.asm.fa'
ASM3='/scratch/ac05869/nih/Plal_v1/plal_v1.asm.fa'
OUT='/scratch/ac05869/plumeria_datanote/KAT'

[ -d $OUT ] || mkdir -p $OUT 

module load SeqKit/2.8.2
ml KAT/2.4.2

seqkit sliding $HIFI/FRA_AM_bc1008.fastq -j 16-W 150 -s 150 -g -o $HIFI/FRA_AM_150bp.fastq #UGA plumeria alba/rubra
seqkit sliding $HIFI/PLU_AD_hifi_reads.fastq -j 16 -W 150 -s 150 -g -o $HIFI/PLU_AD_150bp.fastq #MPI plumeria cubensis/obtusa

#plcu both haps, hap 1, and hap 2
kat comp -t 16 -o $OUT/plcu_v2 -m 21 -h -v ${HIFI}/PLU_AD_150bp.fastq ${ASM2}
kat comp -t 16 -o $OUT/plcu_v2_h1 -m 21 -h -v ${HIFI}/PLU_AD_150bp.fastq ${ASM}/plcu_v2_h1/plcu_v2_h1.fasta
kat comp -t 16 -o $OUT/plcu_v2_h2 -m 21 -h -v ${HIFI}/PLU_AD_150bp.fastq ${ASM}/plcu_v2_h2/plcu_v2_h2.fasta

#plal both haps, hap1, and hap 2
kat comp -t 16 -o $OUT/plal_v1 -m 21 -h -v ${HIFI}/FRA_AM_150bp.fastq ${ASM3}
kat comp -t 16 -o $OUT/plal_v1_h1 -m 21 -h -v ${HIFI}/FRA_AM_150bp.fastq ${ASM}/plal_v1_h1/plal_v1_h1.fasta
kat comp -t 16 -o $OUT/plal_v1_h2 -m 21 -h -v ${HIFI}/FRA_AM_150bp.fastq ${ASM}/plal_v1_h2/plal_v1_h2.fasta



#kat plot spectra-cn -o ${OUT}/gese_v2.asm_cn-plot -x 600 -y 8000000 $WD/gese_v2.asm-main.mx
#sbatch ~/NIH/plumeria_KAT.sh 