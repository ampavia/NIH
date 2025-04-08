#!/bin/bash
#SBATCH --job-name=arima_mapping2	                    # Job name
#SBATCH --partition=highmem_p		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=64		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=250gb			                            # Total memory for job
#SBATCH --time=72:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/err_out/%x_%j.out	# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/err_out/%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#####
#The mapping pipeline will output a single binary alignment map file (BAM file) that contains paired and
#filtered Hi-C paired-end reads mapped to reference sequences. Workflow taken from Arima mapping pipeline
#####
CPU=64
HIC='gel-an_1438201_S3HiC'
IN_DIR='/scratch/ac05869/gelsemium_yahs/gel-an_1438200'
REF='/scratch/ac05869/gelsemium_yahs/gese_v1.asm.fa'
PREFIX='gese_v1.asm.fa'
FAIDX='$REF.fai'
RAW_DIR='/scratch/ac05869/gelsemium_yahs/raw2'
FILT_DIR='/scratch/ac05869/gelsemium_yahs/filtered'
LABEL='gelsemium_yahs'
FILTER='~/NIH/filter_five_end.pl'
COMBINER='~/NIH/two_read_bam_combiner.pl'
STATS='~/NIH/get_stats.pl'

TMP_DIR='/scratch/ac05869/gelsemium_yahs/temp'
PAIR_DIR='/scratch/ac05869/gelsemium_yahs/paired'
#REP_DIR='/scratch/ac05869/gelsemium_yahs/deduplicated'
#REP_LABEL=${LABEL}_rep1
#MERGE_DIR='/path/to/final/merged/alignments/from/any/biological/replicates'
MAPQ_FILTER=10

#module
module load BWA/0.7.17-GCCcore-11.3.0 #will create indeces for reference and map reads against reference 
module load SAMtools/1.16.1-GCC-11.3.0 #will sort the bam file that is being created, then index the bam file

echo "### Step 0: Check output directories’ existence & create them as
needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR

#Run only once. Skip if this step has been completed
echo "### Step 0: Index reference"
bwa index -a bwtsw -p $PREFIX $REF

#Next, we use BWA-MEM to align the Hi-C paired-end reads to reference sequences. Because Hi-C
#captures conformation via proximity-ligated fragments, paired-end reads are first mappedindependently
#(assingle-ends)usingBWA-MEMandaresubsequentlypaired in a later step.
echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -t $CPU $REF $IN_DIR/${HIC}_R1.fastq.gz | samtools view -@ $CPU -Sb - > $RAW_DIR/${HIC}_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t $CPU $REF $IN_DIR/${HIC}_R2.fastq.gz | samtools view -@ $CPU -Sb - > $RAW_DIR/${HIC}_2.bam

#Subsequent to mapping as single-ends, some of these single-end mapped reads can manifest a ligation
#junction and are therefore considered “chimeric” (i.e. they do not originate from a contiguous piece of
#DNA). When BWA-MEM maps these chimeric reads, there can be high quality mapping on both the 5’-
#side and 3’-side of the ligation junction within a given read. In such cases, only the 5’-side should be
#retained because the 3’-side can originate from the same contiguous DNA as the 5’-side of the reads
#mate-pair. Therefore, we retain only the portion of the chimeric read that maps in the 5’-orientation in
#relation to its read orientation. This is accomplished using the script “filter_five_end.pl.”
echo "### Step 2.A: Filter 5' end (1st)"
$SAMTOOLS view -h $RAW_DIR/${HIC}_1.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/${HIC}_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
$SAMTOOLS view -h $RAW_DIR/${HIC}_2.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/${SRA}_2.bam

#After filtering, we pair the filtered single-end Hi-C reads using “two_read_bam_combiner.pl,” which
#outputs a sorted, mapping quality filtered, paired-end BAM file. We then add read groups to this BAM
#file using Picard Tools
module load picard/3.2.0-Java-17

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/${HIC}_1.bam $FILT_DIR/${HIC}_2.bam samtools $MAPQ_FILTER | samtools view -bS -t $FAIDX - | samtools sort -@ $CPU -o $TMP_DIR/$HIC.bam

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar picard AddOrReplaceReadGroups INPUT=$TMP_DIR/$HIC.bam OUTPUT=$PAIR_DIR/$HIC.bam ID=$HIC LB=$HIC SM=$LABEL PL=ILLUMINA PU=none
