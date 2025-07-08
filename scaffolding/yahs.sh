#!/bin/bash

REF='/scratch/ac05869/gese_final_yahs/assembly/gese_v1_organellar_filter.asm.fa' #the contig file
HIC='gel-an_1438201_S3HiC' #same as mapping pipeline script
PREF='gese_v1.org_filter.asm.yahs' #prefix of the output file
YAHS='/scratch/ac05869/gese_final_yahs/yahs' #same as mapping pipeline output directory
JBAT='/scratch/ac05869/gese_final_yahs/juicebox'
STATS='/scratch/ac05869/gese_final_yahs/yahs/filter'
CPU=32
SCAF='/scratch/ac05869/gese_final_yahs/yahs/gese_v1.org_filter.asm.yahs_scaffolds_final.fa'
Min50='/scratch/ac05869/gese_final_yahs/yahs/filter/gese_v1.org_filter.asm.yahs_scaffolds_final_50kb_filter.fa'
Min100='/scratch/ac05869/gese_final_yahs/yahs/filter/gese_v1.org_filter.asm.yahs_scaffolds_final_100kb_filter.fa'
BUSCO='/scratch/ac05869/gese_final_yahs/yahs/busco_downloads'

[ -d $JBAT ] || mkdir -p $JBAT
[ -d $STATS ] || mkdir -p $STATS 

module load YaHS/1.2.2-GCC-11.3.0
module load Juicebox/1.9.9
module load SAMtools/1.16.1-GCC-11.3.0
module load SeqKit/2.9.0
module load BUSCO/5.8.3-foss-2023a

>&2 echo "### Step 1.A: run yahs"
yahs -o ${YAHS}/${PREF} ${REF} ${YAHS}/${HIC}.bam
#####
#The outputs include several AGP format files and a FASTA format file.
#The *_inital_break_[0-9]{2}.agp AGP files are for initial assembly error corrections.
#The *_r[0-9]{2}.agp and related *_r[0-9]{2}_break.agp AGP files are for scaffolding results in each round.
#The *_scaffolds_final.agp and *_scaffolds_final.fa files are for the final scaffolding results.
#####

>&2 echo"### Step 1.B: check stats"
seqkit stats -a -o $STATS/prefilter_stats.txt -t dna -j $CPU $SCAF #finding statistics about your assembly, checking before/after filtering 
busco -i $SCAF -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_$(basename "$SCAF") --out_path $STATS# --download_path $BUSCO #check for completeness of your assembly, before and after filtering 

>&2 echo"### Step 1.C: Filter for 50kb and check stats"
seqkit seq -m 50000 -o $Min50 -t dna -j $CPU $SCAF
seqkit stats -a -o $STATS/50kb_filter_stats.txt -t dna -j $CPU $Min50
busco -i $Min50 -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_$(basename "$Min50") --out_path $STATS# --download_path $BUSCO

>&2 echo"### Step 1.D: Filter for 100kb and check stats"
seqkit seq -m 100000 -o $Min100 -t dna -j $CPU $SCAF
seqkit stats -a -o $STATS/100kb_filter_stats.txt -t dna -j $CPU $Min100
busco -i $Min100 -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_$(basename "$Min100") --out_path $STATS# --download_path $BUSCO

echo "### Step 2.A: generate HiC contact map"
#proceeded with the unfiltered output "_scaffolds_final.fa" which corresponds to "_scaffolds_final.agp"
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
##out_JBAT.txt             - BED format file for hic links
##out_JBAT.liftover.agp    - coordinate file between new and old contigs
##out_JBAT.assembly        - assembly annotation file for Juicebox
##out_JBAT.assembly.agp    - AGP file contains same information as the assembly annotation file. Not a real AGP file as there are no gaps.
##tmp_juicer_pre_JBAT.log  - the output log file
 
echo "### Step 4: need to run juicer_tools pre with 'out_JBAT.txt' (BED file for hic links)."
cat ${JBAT}/tmp_juicer_pre_JBAT.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${JBAT}/${PREF}_JBAT.chrom.sizes
(java -jar -Xmx32G $EBROOTJUICEBOX/juicer_tools.1.9.9.jar pre ${JBAT}/${PREF}_JBAT.txt ${JBAT}/${PREF}_JBAT.hic.part ${JBAT}/${PREF}_JBAT.chrom.sizes) \
&& (mv ${JBAT}/${PREF}_JBAT.hic.part ${JBAT}/${PREF}_JBAT.hic)
#load the JBAT.hic and the JBAT.assembly file into juicebox and manually edit
 
##copy the manually edited review.assembly file from Juicebox on local device into juicebox directory on cluster and
##run ~/NIH/juicer_post_yahs.sh