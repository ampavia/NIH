#!/bin/bash
#SBATCH --job-name=KRT_AF_starsolo	# Job name
#SBATCH --partition=highmem_p		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1			# Run job in single task, by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=24	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-2				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=200G			# Memory per node (30GB); by default using M as unit
#SBATCH --time=128:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=%x_%j.out		# Standard output log, e.g., testBowtie2_12345.out
#SBATCH --error=%x_%j.err		# Standard error log, e.g., testBowtie2_12345.err
#SBATCH --mail-user=ac05869@uga.edu # Send an email when job is done or dead
#SBATCH --mail-type=ALL   	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Single Cell - Align RNAseq Reads to Genome
#       Script function: Align RNAseq Reads to Genome
#       Input: trimmed_reads.fastq
#       Output: alignments.sam --> alignmnets.sorted.bam
################################################################################
LIB='KRT_AF'
WD='/scratch/ac05869/10X_Multiome/KRT_leaf/KRT_AD_AF'
[ -d $WD ] || mkdir -p $WD
cd $WD
R1='/scratch/ac05869/10X_Multiome/KRT_leaf/KRT_AD_AF/gex_paired/KRT_AF_S3_R1_001.fastq.gz'
R2='/scratch/ac05869/10X_Multiome/KRT_leaf/KRT_AD_AF/gex_paired/KRT_AF_S3_R2_001.trim.fastq.gz'

ml STAR/2.7.10b-GCC-11.3.0

#STARsolo alignment
STAR --genomeDir /scratch/ac05869/nih/kratom/star_index_mitr_v1/ \
--readFilesCommand zcat \
--readFilesIn $R2 $R1 \
--runThreadN 24 \
--alignIntronMax 5000 \
--soloUMIlen 12 \
--soloCellFilter EmptyDrops_CR \
--soloFeatures GeneFull \
--soloMultiMappers EM \
--soloType CB_UMI_Simple \
--soloCBwhitelist ../737K-arc-v1.txt \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 26843545600 \
--soloOutFileNames $LIB


#Parameters
#/path/to/STAR --genomeDir /path/to/genome/dir/ --readFilesIn ...  [...other parameters...] --soloType ... --soloCBwhitelist ...
#https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
#--readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read, i.e. R1
#--alignIntronMax 5000: default: 0. maximum intron size, if 0, max intron size will be determined by (2ˆwinBinNbits)*winAnchorDistNbins. This results in a default max of around 500,000
#--soloUMIlen 12: The default barcode lengths (CB=16b, UMI=10b) work for 10X Chromium V2. For V3, specify: --soloUMIlen 12
#--soloCellFilter  EmptyDrops_CR: CellRanger 3.0.0 use advanced filtering based on the EmptyDrop algorithm developed by Lun et al. This algorithm calls extra cells compared to the knee filtering, allowing for cells that have relatively fewer UMIs but are transcriptionally different from the ambient RNA. In STARsolo, this filtering can be activated by:
#--soloFeatures GeneFull: pre-mRNA counts, useful for single-nucleus RNA-seq. This counts all read that overlap gene loci, i.e. included both exonic and intronic reads:
#	10x now recommends using intronic counts going forward, so I've turned this option on
#The multi-gene read recovery options are specified with --soloMultiMappers. Several algorithms are implemented:
#	--soloMultiMappers Uniform: uniformly distributes the multi-gene UMIs to all genes in its gene set. Each gene gets a fractional count of 1/N_genes, where N_genes is the number of genes in the set. This is the simplest possible option, and it offers higher sensitivity for gene detection at the expense of lower precision.
#	--soloMultiMappers EM: uses Maximum Likelihood Estimation (MLE) to distribute multi-gene UMIs among their genes, taking into account other UMIs (both unique- and multi-gene) from the same cell (i.e. with the same CB). Expectation-Maximization (EM) algorithm is used to find the gene expression values that maximize the likelihood function. Recovering multi-gene reads via MLE-EM model was previously used to quantify transposable elements in bulk RNA-seq {TEtranscripts} and in scRNA-seq {Alevin; Kallisto-bustools}.
#--outSAMtype: Output in BAM sorted by coordinate
#--readFilesCommand the files are gzipped

#Additional parameters to consider based on https://cumulus.readthedocs.io/en/latest/starsolo.html
#--soloCBstart 1 #not included
#--soloCBlen 16 #not included
#--soloUMIstart 17 #not included
#--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts #not included
#--soloUMIfiltering MultiGeneUMI_CR #not included
#--soloUMIdedup 1MM_CR #not included
#--clipAdapterType CellRanger4 #not included
#--outFilterScoreMin 30 #not included
#--outSAMattributes CR UR CY UY CB UB #not included


