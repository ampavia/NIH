#!/bin/bash
#SBATCH --job-name=yahs_100kb	                    # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=32		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=250gb			                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

PREF='gese.org_filter.yahs_scaf.100kb_filter.yahs2' #prefix of the output file
HIC='gel-an_1438201_S3HiC'
JBAT='/scratch/ac05869/gese_final_yahs/juicebox2'
CPU=32
REF='/scratch/ac05869/gese_final_yahs/scaf_assembly/gese_v1.org_filter.asm.yahs_scaffolds_final_100kb_filter.fa'
ASM='gese_v1.org_filter.asm.yahs_scaffolds_final_100kb_filter.fa'
YAHS='/scratch/ac05869/gese_final_yahs/yahs2'

[ -d $JBAT ] || mkdir -p $JBAT

module load YaHS/1.2.2-GCC-11.3.0
module load Juicebox/1.9.9
module load SAMtools/1.16.1-GCC-11.3.0

#>&2 echo "### Step 1: run yahs"
#yahs -o ${YAHS}/${PREF} ${REF} ${YAHS}/${HIC}.bam
#####
#The outputs include several AGP format files and a FASTA format file.
#The *_inital_break_[0-9]{2}.agp AGP files are for initial assembly error corrections.
#The *_r[0-9]{2}.agp and related *_r[0-9]{2}_break.agp AGP files are for scaffolding results in each round.
#The *_scaffolds_final.agp and *_scaffolds_final.fa files are for the final scaffolding results.
#####

#>&2 echo "### Step 1: generate HiC contact map"

#(juicer pre ${YAHS}/${PREF}.bin ${YAHS}/${PREF}_scaffolds_final.agp ${REF}.fai \
#2>${YAHS}/tmp_juicer_pre.log \
#| LC_ALL=C sort -k2,2d -k6,6d -T ${YAHS} --parallel=8 -S32G \
#| awk 'NF' > ${YAHS}/alignments_sorted.txt.part) \
#&& (mv ${YAHS}/alignments_sorted.txt.part ${YAHS}/alignments_sorted.txt)

>&2 echo "### Step 2: make scaffolds_final.chrom.sizes for juicer_tools input file"
#####
#The file for scaffold sizes should contain two columns - scaffold name and scaffold size, which can be taken from the first two columns of the FASTA index file 
#or the log file created in the previous step
#####
cat ${YAHS}/tmp_juicer_pre.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${YAHS}/scaffolds_final.chrom.sizes

>&2 echo "### Step 3: do juicer hic map"
(java -jar -Xmx32G $EBROOTJUICEBOX/juicer_tools.1.9.9.jar pre ${YAHS}/alignments_sorted.txt ${YAHS}/${PREF}.hic.part ${YAHS}/scaffolds_final.chrom.sizes) \
&& (mv ${YAHS}/${PREF}.hic.part ${YAHS}/${PREF}.hic)

>&2 echo "### Step 4: generate HiC contact map file that can be loaded by Juicebox for manual editing - assembly (JBAT) mode (-a)"

juicer pre -a -o ${JBAT}/${PREF}_JBAT ${YAHS}/${PREF}.bin ${YAHS}/${PREF}_scaffolds_final.agp ${REF}.fai 2>${JBAT}/tmp_juicer_pre_JBAT.log
##This results in 5 output files:
#out_JBAT.txt             - BED format file for hic links
#out_JBAT.liftover.agp    - coordinate file between new and old contigs
#out_JBAT.assembly        - assembly annotation file for Juicebox
#out_JBAT.assembly.agp    - AGP file contains same information as the assembly annotation file. Not a real AGP file as there are no gaps.
#tmp_juicer_pre_JBAT.log  - the output log file

>&2 echo "### Step 5: need to run juicer_tools pre with 'out_JBAT.txt' (BED file for hic links)."
cat ${JBAT}/tmp_juicer_pre_JBAT.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${JBAT}/${PREF}_JBAT.chrom.sizes
(java -jar -Xmx32G $EBROOTJUICEBOX/juicer_tools.1.9.9.jar pre ${JBAT}/${PREF}_JBAT.txt ${JBAT}/${PREF}_JBAT.hic.part ${JBAT}/${PREF}_JBAT.chrom.sizes) \
&& (mv ${JBAT}/${PREF}_JBAT.hic.part ${JBAT}/${PREF}_JBAT.hic)

# >&2 echo "### Final step: generate final genome assembly file after manual curation with JuiceBox (JBAT)"
## the output assembly file after curation is ${JBAT}/${PREF}_JBAT.review.assembly
## the final output is ${JBAT}/${PREF}_JBAT.FINAL.agp and ${JBAT}/${PREF}_JBAT.FINAL.fa
# juicer post -o ${JBAT}/${PREF}_JBAT ${JBAT}/${PREF}_JBAT.review.assembly ${JBAT}/${PREF}_JBAT.liftover.agp ${REF}
