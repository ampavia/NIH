#!/bin/bash

REF='/scratch/ac05869/gelsemium_yahs/gese_v1.asm.fa' #the contig file
HIC='gel-an_1438201_S3HiC' #same as mapping pipeline script
PREF='gese_v1.asm.yahs' #prefix of the output file
YAHS='/scratch/ac05869/gelsemium_yahs/yahs' #same as mapping pipeline output directory
JBAT='/scratch/ac05869/gelsemium_yahs/juicebox'
[ -d $JBAT ] || mkdir -p $JBAT

module load YaHS/1.2.2-GCC-11.3.0
module load Juicebox/1.9.9
module load SAMtools/1.16.1-GCC-11.3.0

echo "### Step 1: run yahs"
yahs -o ${YAHS}/${PREF} ${REF} ${YAHS}/${HIC}.bam

echo "### Step 2.A: generate HiC contact map"
(juicer pre ${YAHS}/${PREF}.bin ${YAHS}/${PREF}_scaffolds_final.agp ${REF}.fai \
2>${YAHS}/tmp_juicer_pre.log \
| LC_ALL=C sort -k2,2d -k6,6d -T ${YAHS} --parallel=8 -S32G \
| awk 'NF' > ${YAHS}/alignments_sorted.txt.part) \
&& (mv ${YAHS}/alignments_sorted.txt.part ${YAHS}/alignments_sorted.txt)

echo "### Step 2.B: make scaffolds_final.chrom.sizes for juicer_tools input file"
#####
#The file for scaffold sizes should contain two columns - scaffold name and scaffold size, which can be taken from the first two columns of the FASTA index file 
#or the log file created in the previous step
#####
cat ${YAHS}/tmp_juicer_pre.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${YAHS}/${PREF}_scaffolds_final.chrom.sizes

echo "### Step 2.C: do juicer hic map"
(java -jar -Xmx32G $EBROOTJUICEBOX/juicer_tools.1.9.9.jar pre ${YAHS}/alignments_sorted.txt ${YAHS}/${PREF}.hic.part ${YAHS}/${PREF}_scaffolds_final.chrom.sizes) && (mv ${YAHS}/${PREF}.hic.part ${YAHS}/${PREF}.hic)

echo "### Step 3: generate HiC contact map file that can be loaded by Juicebox for manual editing - assembly (JBAT) mode (-a)"

juicer pre -a -o ${JBAT}/${PREF}_JBAT ${YAHS}/${PREF}.bin ${YAHS}/${PREF}_scaffolds_final.agp ${REF}.fai 2>${JBAT}/tmp_juicer_pre_JBAT.log
##This results in 5 output files:
#out_JBAT.txt             - BED format file for hic links
#out_JBAT.liftover.agp    - coordinate file between new and old contigs
#out_JBAT.assembly        - assembly annotation file for Juicebox
#out_JBAT.assembly.agp    - AGP file contains same information as the assembly annotation file. Not a real AGP file as there are no gaps.
#tmp_juicer_pre_JBAT.log  - the output log file

echo "### Step 4: need to run juicer_tools pre with 'out_JBAT.txt' (BED file for hic links)."
cat ${JBAT}/tmp_juicer_pre_JBAT.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${JBAT}/${PREF}_JBAT.chrom.sizes
(java -jar -Xmx32G $EBROOTJUICEBOX/juicer_tools.1.9.9.jar pre ${JBAT}/${PREF}_JBAT.txt ${JBAT}/${PREF}_JBAT.hic.part ${JBAT}/${PREF}_JBAT.chrom.sizes) \
&& (mv ${JBAT}/${PREF}_JBAT.hic.part ${JBAT}/${PREF}_JBAT.hic)

# echo "### Final step: generate final genome assembly file after manual curation with JuiceBox (JBAT)"
## the output assembly file after curation is ${JBAT}/${PREF}_JBAT.review.assembly
## the final output is ${JBAT}/${PREF}_JBAT.FINAL.agp and ${JBAT}/${PREF}_JBAT.FINAL.fa
# juicer post -o ${JBAT}/${PREF}_JBAT ${JBAT}/${PREF}_JBAT.review.assembly ${JBAT}/${PREF}_JBAT.liftover.agp ${REF}
