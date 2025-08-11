#!/bin/bash
#SBATCH --job-name=ATAC_alignment
#SBATCH --partition=highmem_p		# Partition name (batch, heighten_p, or gpu_p), _required_
#SBATCH --ntasks=1 		# Run job in single task or in paralelle, _required_
#SBATCH --cpus-per-task=16		# CPU cores per task
#SBATCH --mem=300G			# How much memory per node, _required_
#SBATCH --time=72:00:00		# Time Limit hrs:min:sec or day-hrs:min:sec 2-12:00:00 is 2.5 d, _required_
#SBATCH --export=NONE		# Don't export submit node variables to compute node
#SBATCH --output=/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis/err_out/%x_%j.out	# Standard output log
#SBATCH --error=/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu # Send an email when job is done or dead
#SBATCH --mail-type=ALL	# Mail events (BEGIN, END, FAIL, ALL)

###Post alignment QC: Remove duplicates & low-quality alignments; Shift read coordinates; Bam visualisation
# FASTQC did indicate adaptors were present. Did not process ATAC reads

cd /scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis
OUT='/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis/data'
OUT2='/scratch/ac05869/10X_Multiome/KRT_leaf/ATAC_analysis/deeptools_out'
bt2idx='/scratch/ac05869/10X_Multiome/KRT_leaf/genome'
ATAC1='../KRT_AC_AE/atac_libs'
ATAC2='../KRT_AD_AF/atac_libs'
[ -d $OUT ] || mkdir -p $OUT 
[ -d $OUT2 ] || mkdir -p $OUT2

ml Bowtie2/2.5.2-GCC-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0
ml picard/3.2.0-Java-17
ml deepTools/3.5.2-foss-2022a

#index
bowtie2-build --large-index $bt2idx/mitr_v1.asm.fa $bt2idx/mitr_v1.asm_bt2
#map
bowtie2 --local --very-sensitive --no-mixed --no-discordant \
-X 700 \
-x $bt2idx/mitr_v1.asm_bt2 \
-1 $ATAC1/KRT_AC_1_S8_R1_001.fastq.gz,$ATAC1/KRT_AC_2_S9_R1_001.fastq.gz,$ATAC1/KRT_AC_3_S10_R1_001.fastq.gz,$ATAC1/KRT_AC_4_S11_R1_001.fastq.gz \
-2 $ATAC1/KRT_AC_1_S8_R2_001.fastq.gz,$ATAC1/KRT_AC_2_S9_R2_001.fastq.gz,$ATAC1/KRT_AC_3_S10_R2_001.fastq.gz,$ATAC1/KRT_AC_4_S11_R2_001.fastq.gz \
| samtools view -bS - > $OUT/KRT_AC_all.bam

bowtie2 --local --very-sensitive --no-mixed --no-discordant \
-X 700 \
-x $bt2idx/mitr_v1.asm_bt2 \
-1 $ATAC2/KRT_AD_1_S12_R1_001.fastq.gz,$ATAC2/KRT_AD_2_S13_R1_001.fastq.gz,$ATAC2/KRT_AD_3_S14_R1_001.fastq.gz,$ATAC2/KRT_AD_4_S15_R1_001.fastq.gz \
-2 $ATAC2/KRT_AD_1_S12_R2_001.fastq.gz,$ATAC2/KRT_AD_2_S13_R2_001.fastq.gz,$ATAC2/KRT_AD_3_S14_R2_001.fastq.gz,$ATAC2/KRT_AD_4_S15_R2_001.fastq.gz \
| samtools view -bS - > $OUT/KRT_AD_all.bam

#Sort the output bam file by coordinate
for bam in *_all.bam
do
samtools sort $OUT/$bam -o $OUT/${bam/.bam/}_sorted.bam 
done

#Generate an index file
for bam in *_sorted.bam 
do
samtools index $OUT/$bam
done

# Mark duplicates
for bam in *_sorted.bam
do
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar MarkDuplicates \
QUIET=true INPUT=$OUT/$bam OUTPUT=$OUT/${bam/.bam/}.marked.bam \
METRICS_FILE=$OUT/${bam/.bam/}.dup.metrics \
REMOVE_DUPLICATES=false \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=.
done

#View the % of duplicates
for metrics in *_sorted.dup.metrics
do
head -n 8 $OUT/$metrics | cut -f 7,9 | grep -v ^ | tail -n 2 > $OUT/${metrics}.pct.dup
done

for bam in *_sorted.marked.bam
do
samtools view -q 30 -c $OUT/$bam >> $OUT/qscores.below.30.txt #divide this number by 2 to calculate the # of DNA fragments
done

# Remove low-quality alignments
for bam in *_sorted.marked.bam
do
samtools view -h -b -f 2 -F 1548 $OUT/$bam | samtools sort -n -o $OUT/${bam/.bam/}.filtered.bam
done

# Shift read coordinates
for bam in *_sorted.marked.filtered.bam
do
alignmentSieve --numberOfProcessors 16 --ATACshift --bam $OUT/$bam -o $OUT/${bam/.bam/}.shifted.bam
samtools index $OUT/KRT_A?_all_sorted.marked.filtered.shifted.bam
done

# Bam visualization
for bam in *_sorted.marked.filtered.shifted.bam
do
bamCoverage --numberOfProcessors 16 \
--normalizeUsing CPM \
--effectiveGenomeSize 1433922938 \
--bam $OUT/$bam \
-o $OUT2/${bam/_sorted.marked.filtered.shifted.bam/}_coverage_CPM.bw
--outFileFormat bigwig
done

#sbatch ~/NIH/multiome/ATAC_mapping_bowtie.sh

###Bowtie2 options:
# The local parameter is used to 'soft clip' the end of reads to allow the best possible alignment,
# including any remaining adapter sequences (e.g. 1 or 2bp).
# By using the --no-mixed and --no-discordant parameters, reads will only be aligned if 
# both reads align successfully as a pair (this avoids the need to later remove reads which are not properly paired,
# which is a common post-alignment QC step). 
# The -I 25 and -X 51 require fragments to be greater than I and less than X (multiome generated 50bp reads)

##Removing low-quality alignments
#Kept multimapping
#The following code uses the sam/bam flags to retain properly mapped pairs (-f 2) and to remove reads which fail the platform/vendor QC checks (-F 512), duplicate reads (-F 1024) and those which are unmapped (-F 12).
#The three flags to be removed can be combined into -F 1548, which will remove reads which meet any of the three individual flags

## Shifting read coordinates
#An optional step in analysing data generated using the Tn5 transposase
# (such as ATAC-seq, ChIPmentation etc.) is to account for a small DNA insertion,
# introducted as repair of the transposase-induced nick introduces a 9bp insertion.
# Reads aligning to the + strand should be offset by +4bp and reads aligned to the -ve strand
# should be offset by -5bp. For references, see the first ATAC-seq paper by Buenrostro et al., (2013)
# and the analysis by Adey et al., (2010) which showed this insertion bias.
# Shifting coordinates is only really important if single-base resolution is required,
# for example in the analysis of transcription factor motifs in ATAC-seq peak footprints.
# Be aware that some tools do this shifting themselves (so double check manuals!).