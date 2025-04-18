#!/bin/bash

REF='/scratch/ac05869/gese_final_yahs/assembly/gese_v1_organellar_filter.asm.fa' #the contig file
HIC='gel-an_1438201_S3HiC' #same as mapping pipeline script
PREF='gese_v1.org_filter.asm.yahs' #prefix of the output file
YAHS='/scratch/ac05869/gese_final_yahs/yahs' #same as mapping pipeline output directory
JBAT='/scratch/ac05869/gese_final_yahs/juicebox'
STATS='/scratch/ac05869/gese_final_yahs/yahs/filter'
CPU=32
Min50='/scratch/ac05869/gese_final_yahs/yahs/filter/gese_v1.org_50kb_filter.asm.yahs.fa'
Min100='/scratch/ac05869/gese_final_yahs/yahs/filter/gese_v1.org_100kb_filter.asm.yahs.fa'
SCAF='/scratch/ac05869/gese_final_yahs/yahs/gese_v1.org_filter_scaffolds_final.fa'
[ -d $JBAT ] || mkdir -p $JBAT
[ -d $STATS ] || mkdir -p $STATS 

module load YaHS/1.2.2-GCC-11.3.0
module load Juicebox/1.9.9
module load SAMtools/1.16.1-GCC-11.3.0
module load SeqKit/2.9.0
module load BUSCO/5.8.3-foss-2023a

>&2 echo "### Step 1: run yahs"
yahs -o ${YAHS}/${PREF} ${REF} ${YAHS}/${HIC}.bam
#####
#The outputs include several AGP format files and a FASTA format file.
#The *_inital_break_[0-9]{2}.agp AGP files are for initial assembly error corrections.
#The *_r[0-9]{2}.agp and related *_r[0-9]{2}_break.agp AGP files are for scaffolding results in each round.
#The *_scaffolds_final.agp and *_scaffolds_final.fa files are for the final scaffolding results.
#####

>&2 echo"### Step 1.1: check stats"
seqkit stats -a -o $STATS/prefilter_stats.txt -t dna -j $CPU $SCAF #finding statistics about your assembly, checking before/after filtering 
busco -i $SCAF -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_$(basename "$REF") --out_path $STATS #check for completeness of your assembly, before and after filtering 

>&2 echo"### Step 1.2: Filter for 50kb and check stats"
seqkit seq -m 50000 -o $Min50 -t dna -j $CPU $SCAF
seqkit stats -a -o $STATS/50kb_filter_stats.txt -t dna -j $CPU $Min50
busco -i $Min50 -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_$(basename "$Min50") --out_path $STATS

>&2 echo"### Step 1.2: Filter for 100kb and check stats"
seqkit seq -m 100000 -o $Min100 -t dna -j $CPU $SCAF
seqkit stats -a -o $STATS/100kb_filter_stats.txt -t dna -j $CPU $Min100
busco -i $Min100 -m genome -l eudicotyledons_odb12 -c $CPU -o BUSCO_$(basename "$Min100") --out_path $STATS

