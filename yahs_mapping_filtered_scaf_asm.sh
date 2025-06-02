#!/bin/bash
#SBATCH --job-name=mapping_filtered	                    # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=32		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=250gb			                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=BEGIN,END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#####
#The mapping pipeline will output a single binary alignment map file (BAM file) that contains paired and
#filtered Hi-C paired-end reads mapped to reference sequences. Workflow taken from Arima mapping pipeline
#####
WD='/scratch/ac05869/gese_final_yahs/scaf_assembly'
CPU=32
HIC='gel-an_1438201_S3HiC'
IN_DIR='/scratch/ac05869/gelsemium_yahs/gel-an_1438200'
REF='/scratch/ac05869/gese_final_yahs/scaf_assembly/gese_v1.org_filter.asm.yahs_scaffolds_final_100kb_filter.fa'
ASM='gese_v1.org_filter.asm.yahs_scaffolds_final_100kb_filter.fa'
FAIDX='/scratch/ac05869/gese_final_yahs/scaf_assembly/gese_v1.org_filter.asm.yahs_scaffolds_final_100kb_filter.fai'
RAW_DIR='/scratch/ac05869/gese_final_yahs/raw2'
FILT_DIR='/scratch/ac05869/gese_final_yahs/filtered2'
LABEL='gese_v1.org_filter.asm.yahs_scaffolds_final_100kb_filter'
FILTER='/home/ac05869/NIH/filter_five_end.pl'
COMBINER='/home/ac05869/NIH/two_read_bam_combiner.pl'
STATS='/home/ac05869/NIH/get_stats.pl'
TMP_DIR='/scratch/ac05869/gese_final_yahs/temp2'
PAIR_DIR='/scratch/ac05869/gese_final_yahs/paired2'
YAHS='/scratch/ac05869/gese_final_yahs/yahs2'
#REP_LABEL=${LABEL}_rep1
#MERGE_DIR='/path/to/final/merged/alignments/from/any/biological/replicates'
MAPQ_FILTER=10

#modules
module load BWA/0.7.17-GCCcore-11.3.0 #will create indeces for reference and map reads against reference 
module load SAMtools/1.16.1-GCC-11.3.0 #will also do index of the genome to be used in this pipeline and yahs
module load picard/3.2.0-Java-17

>&2 echo "### Step 0: Check output directories’ existence & create them as
needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $YAHS ] || mkdir -p $YAHS

cd $WD

##Run only once. Skip if this step has been completed
echo "### Step 0: Index reference"
bwa index -a bwtsw -p $ASM $REF

#####
#Next, we use BWA-MEM to align the Hi-C paired-end reads to reference sequences. Because Hi-C
#captures conformation via proximity-ligated fragments, paired-end reads are first mapped independently
#(assingle-ends)usingBWA-MEM and are subsequently paired in a later step.
#####
>&2 echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -t $CPU $REF $IN_DIR/${HIC}_R1.fastq.gz | samtools view -@ $CPU -Sb - > $RAW_DIR/${HIC}_1.bam

>&2 echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t $CPU $REF $IN_DIR/${HIC}_R2.fastq.gz | samtools view -@ $CPU -Sb - > $RAW_DIR/${HIC}_2.bam

#####
#Subsequent to mapping as single-ends, some of these single-end mapped reads can manifest a ligation
#junction and are therefore considered “chimeric” (i.e. they do not originate from a contiguous piece of
#DNA). When BWA-MEM maps these chimeric reads, there can be high quality mapping on both the 5’-
#side and 3’-side of the ligation junction within a given read. In such cases, only the 5’-side should be
#retained because the 3’-side can originate from the same contiguous DNA as the 5’-side of the reads
#mate-pair. Therefore, we retain only the portion of the chimeric read that maps in the 5’-orientation in
#relation to its read orientation. This is accomplished using the script “filter_five_end.pl.”
#####

>&2 echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h $RAW_DIR/${HIC}_1.bam | perl $FILTER | samtools view -Sb - > $FILT_DIR/${HIC}_1.bam

>&2 echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h $RAW_DIR/${HIC}_2.bam | perl $FILTER | samtools view -Sb - > $FILT_DIR/${HIC}_2.bam

#####
#After filtering, we pair the filtered single-end Hi-C reads using “two_read_bam_combiner.pl,” which
#outputs a sorted, mapping quality filtered, paired-end BAM file. We then add read groups to this BAM
#file using Picard Tools
#####

>&2 echo "### Step 3.A: Pair reads & mapping quality filter"
samtools faidx ${REF}
perl $COMBINER $FILT_DIR/${HIC}_1.bam $FILT_DIR/${HIC}_2.bam samtools $MAPQ_FILTER | samtools view -bS -t $FAIDX - | samtools sort -@ $CPU -o $TMP_DIR/${HIC}.bam

>&2 echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
INPUT=$TMP_DIR/${HIC}.bam OUTPUT=$PAIR_DIR/${HIC}.bam ID=$HIC LB=$HIC SM=$LABEL PL=ILLUMINA PU=none

##### No replicates in this case, so the replicate directory was renamed for yahs input file
#Note, that if you perform merging of technical replicates, then the file names and locations will
#change from the written flow of this pipeline. You will need to adjust the file names and locations that are
#used as input in the following step - PCR duplicate removal.
#####

>&2 echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar MarkDuplicates \
INPUT=$PAIR_DIR/${HIC}.bam OUTPUT=$YAHS/${HIC}.bam \
METRICS_FILE=$YAHS/metrics.${HIC}.txt TMP_DIR=$TMP_DIR \
ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index $YAHS/${HIC}.bam

perl $STATS $YAHS/${HIC}.bam > $YAHS/${HIC}.bam.stats
>&2 echo "Finished Mapping Pipeline through Duplicate Removal"

#####
#The final output of this pipeline is a single BAM
#file that contains the paired, 5’-filtered, and duplicate-removed Hi-C reads mapped to the reference
#sequences of choice. The resulting statistics file has a breakdown of the total number of intra-contig
#read-pairs, long-range intra-contig read-pairs, and inter-contig read-pairs in the final processed BAM
#file.
#####

#bash ~/NIH/yahs_filtered_asm.sh