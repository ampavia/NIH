#!/bin/bash
#SBATCH --job-name=BUSCO                   # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=32		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=300G		                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)
STATS='/scratch/ac05869/plu_comparative_final/genomes/BUSCO_plumeria_haplotypes/final_stats'
BUSCO='/scratch/ac05869/plu_comparative_final/genomes/BUSCO_plumeria_haplotypes/busco_downloads'
ASM='/scratch/ac05869/plu_comparative_final/genomes'
CPU=32
LIN='embryophyta_odb10'
PEP='/scratch/ac05869/plu_comparative_final/genomes/repr_pep'

wd='/scratch/ac05869/plu_comparative_final/genomes/BUSCO_plumeria_haplotypes'
cd $wd

module load SeqKit/2.9.0
module load BUSCO/5.8.3-foss-2023a

>&2 echo"### Step 1: check stats"
seqkit stats -a -o $STATS/plob_v4.txt -t dna -j $CPU $ASM/plob_v4/plob_v4.asm.fa #finding statistics about total obtusa assembly
seqkit stats -a -o $STATS/Pobtusa_A.txt -t dna -j $CPU $ASM/Pobtusa_A/Pobtusa_A.asm.fa #finding statistics about obtusa assembly A
seqkit stats -a -o $STATS/Pobtusa_B.txt -t dna -j $CPU $ASM/Pobtusa_B/Pobtusa_B.asm.fa #finding statistics about obtusa assembly B

seqkit stats -a -o $STATS/plru_v2.txt -t dna -j $CPU $ASM/plru_v2/plru_v2.asm.fa #finding statistics about total rubra assembly
seqkit stats -a -o $STATS/Prubra_A.txt -t dna -j $CPU $ASM/Prubra_A/Prubra_A.asm.fa #finding statistics about obtusa assembly A
seqkit stats -a -o $STATS/Prubra_B.txt -t dna -j $CPU $ASM/Prubra_B/Prubra_B.asm.fa #finding statistics about obtusa assembly B


>&2 echo"### Step 2: BUSCO embryophyta_odb10"
busco -i $ASM/plob_v4/plob_v4.asm.fa -m genome -l $LIN -c $CPU -o plob_v4 --out_path $wd --download_path $BUSCO
busco -i $ASM/Pobtusa_A/Pobtusa_A.asm.fa -m genome -l $LIN -c $CPU -o Pobtusa_A --out_path $wd --download_path $BUSCO
busco -i $ASM/Pobtusa_B/Pobtusa_B.asm.fa -m genome -l $LIN -c $CPU -o Pobtusa_B --out_path $wd --download_path $BUSCO

busco -i $ASM/plru_v2/plru_v2.asm.fa -m genome -l $LIN -c $CPU -o plru_v2 --out_path $wd --download_path $BUSCO
busco -i $ASM/Prubra_A/Prubra_A.asm.fa -m genome -l $LIN -c $CPU -o Prubra_A --out_path $wd --download_path $BUSCO
busco -i $ASM/Prubra_B/Prubra_B.asm.fa -m genome -l $LIN -c $CPU -o Prubra_B --out_path $wd --download_path $BUSCO

>&2 echo"### Step 2: BUSCO embryophyta_odb10 repr hc gene models"
busco -i $ASM/plob_v4/plob_v4.hc_gene_models.repr.pep.fa -m prot -l $LIN -c $CPU -o pep_plob_v4 --out_path $wd --download_path $BUSCO
busco -i $PEP/Pobtusa_A_repr_hc_pep.fa -m prot -l $LIN -c $CPU -o pep_Pobtusa_A --out_path $wd --download_path $BUSCO
busco -i $PEP/Pobtusa_B_repr_hc_pep.fa -m prot -l $LIN -c $CPU -o pep_Pobtusa_B --out_path $wd --download_path $BUSCO

busco -i $ASM/plru_v2/plru_v2.hc_gene_models.repr.pep.fa -m prot -l $LIN -c $CPU -o pep_plru_v2 --out_path $wd --download_path $BUSCO
busco -i $PEP/Prubra_A_repr_hc_pep.fa -m prot -l $LIN -c $CPU -o pep_Prubra_A --out_path $wd --download_path $BUSCO
busco -i $PEP/Prubra_B_repr_hc_pep.fa -m prot -l $LIN -c $CPU -o pep_Prubra_B --out_path $wd --download_path $BUSCO


#sbatch --export=CPU=32 ~/NIH/Plumeria_genomes/BUSCO.sh