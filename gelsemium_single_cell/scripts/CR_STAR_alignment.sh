#!/bin/bash
#SBATCH --job-name=CR_STARsolo		# Job name 
#SBATCH --partition=batch		# Partition name (batch, highmem_p, or gpu_p)
#SBATCH --ntasks=1		# aka threads. Each task by default using 1 CPU core on a single node
#SBATCH --cpus-per-task=24	 	# CPU core count per task, by default 1 CPU core per task
#SBTACH --array=1-4				# Array element range from 0 to 1, i.e. 2 element jobs
#SBATCH --mem=100GB			# Memory per node (30GB); by default using M as unit
#SBATCH --time=48:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.out		# Standard output log
#SBATCH --error=/scratch/ac05869/gelsemium_sc/err_out/%x_%j.err		# Standard error log
#SBATCH --mail-user=ac05869@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL          	# Mail events (BEGIN, END, FAIL, ALL)

################################################################################
#Project: Single Cell processing for multiome
#       Script function: mapping R2 and R1 to genome index
#       Input: genome index, barcodes, and fastq reads 2 and 1 in that order
#       Output: BAM
################################################################################
LIB=`head -n ${SLURM_ARRAY_TASK_ID} ${INFILE} | cut -f 1 | tail -n 1`
IN='/scratch/ac05869/gelsemium_sc/GEL/raw'
OUT='/scratch/ac05869/gelsemium_sc/GEL/STAR'

if [ ! -d ${OUT} ]
then
    mkdir -p ${OUT}
fi

mkdir ${OUT}/${LIB}
cd ${OUT}/${LIB}

ml STAR/2.7.11b-GCC-13.3.0
ml CellRanger-ARC/2.0.2

#cp /apps/eb/CellRanger-ARC/2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz ${OUT}
#gunzip ${OUT}/737K-arc-v1.txt.gz

STAR --genomeDir /scratch/ac05869/gelsemium_sc/genomes/gese_mpi_v2 \
--readFilesCommand zcat \
--readFilesIn ${IN}/${LIB}*R2_001.fastq.gz \
${IN}/${LIB}*R1_001.fastq.gz \
--runThreadN 24 \
--alignIntronMax 5000 \
--soloUMIlen 12 \
--soloCellFilter EmptyDrops_CR \
--soloFeatures GeneFull \
--soloMultiMappers EM \
--soloType CB_UMI_Simple \
--soloCBwhitelist ${OUT}/737K-arc-v1.txt \
--outSAMtype BAM SortedByCoordinate \
--soloOutFileNames $LIB

#Parameters
#/path/to/STAR --genomeDir /path/to/genome/dir/ --readFilesIn ...  [...other parameters...] --soloType ... --soloCBwhitelist ...
#https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
#--readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read, i.e. R1
#--alignIntronMax 5000: default: 0. maximum intron size, if 0, max intron size will be determined by (2ˆwinBinNbits)*winAnchorDistNbins. This results in a default max of around 500,000
#--soloUMIlen 12: The default barcode lengths (CB=16b, UMI=10b) work for 10X Chromium V2. For V3, specify: --soloUMIlen 12
#--soloCellFilter  EmptyDrops_CR: CellRanger 3.0.0 use advanced filtering based on the EmptyDrop algorithm developed by Lun et al. This algorithm calls extra cells compared to the knee filtering, allowing for cells that have relatively fewer UMIs but are transcriptionally different from the ambient RNA. In STARsolo, this filtering can be activated by:
#--soloFeatures GeneFull: pre-mRNA counts, useful for single-nucleus RNA-seq. This counts all read that overlap gene loci, i.e. included both exonic and intronic reads:
#       10x now recommends using intronic counts going forward, so I've turned this option on
#The multi-gene read recovery options are specified with --soloMultiMappers. Several algorithms are implemented:
#       --soloMultiMappers Uniform: uniformly distributes the multi-gene UMIs to all genes in its gene set. Each gene gets a fractional count of 1/N_genes, where N_genes is the number of genes in the set. This is the simplest possible option, and it offers higher sensitivity for gene detection at the expense of lower precision.
#       --soloMultiMappers EM: uses Maximum Likelihood Estimation (MLE) to distribute multi-gene UMIs among their genes, taking into account other UMIs (both unique- and multi-gene) from the same cell (i.e. with the same CB). Expectation-Maximization (EM) algorithm is used to find the gene expression values that maximize the likelihood function. Recovering multi-gene reads via MLE-EM model was previously used to quantify transposable elements in bulk RNA-seq {TEtranscripts} and in scRNA-seq {Alevin; Kallisto-bustools}.
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

#sbatch --array 1-4 --export=INFILE=/scratch/ac05869/gelsemium_sc/GEL_libraries.txt /scratch/ac05869/gelsemium_sc/scripts/CR_STAR_alignment.sh