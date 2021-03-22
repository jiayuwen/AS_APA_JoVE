#! usr/bin/sh

#Title: Data downloading and pre-processing
#Running this file requires polyA_annotation.fwd.bed,  polyA_annotation.rev.bed 
#To run- nohup sh Data_PreProcessing.3PSeq.sh

#*************************Downloading Raw Data[SRA]*************************# 

## Download raw data from SRA using GNU parallel and SRA-toolkit

#Creat project directory to store raw data
PROJ_DIR=AS_analysis/
mkdir $PROJ_DIR
cd $PROJ_DIR

RAW_DATA=raw_data
mkdir -p $RAW_DATA
cd $RAW_DATA
# Download SRA files using accession numbers
seq 1553129 1553136 | parallel prefetch SRR{}
  
#Convert SRA to fastq
parallel -j 3 fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-e --clip --origfmt {} ::: SRR1553129/SRR1553129.sra SRR1553130/SRR1553130.sra SRR1553131/SRR1553131.sra SRR1553132/SRR1553132.sra SRR1553133/SRR1553133.sra SRR1553134/SRR1553134.sra SRR1553135/SRR1553135.sra SRR1553136/SRR1553136.sra

cd ..

#*************************PRE-PROCESSING*************************# 
# Unzip the fastq.gz data to gain fastq files
parallel “gunzip {}” :::$RAW_DATA/*.fastq.gz

# Create multiple directories to store data generated in different pre-processing steps
TRIM_DIR=adapter.trimming.results
mkdir $TRIM_DIR

TRIMRC_DIR=trimmed.rc.data
mkdir $TRIMRC_DIR

SAM_DIR=samfiles
mkdir $SAM_DIR

BAM_DIR=bamfiles
mkdir $BAM_DIR

BDG_DIR=bedgraphfiles
mkdir $BDG_DIR

BW_DIR=BigWigfiles
mkdir $BW_DIR

CM_DIR=compute.matrix.results
mkdir $CM_DIR

HM_DIR=plotHeatmap.results
mkdir $HM_DIR



# Download the bowtie Index of mm10
wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip
unzip mm10.zip

# Download the mm10.chrom.sizes
wget https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes



for fq in $RAW_DATA/*.fastq
do
OUTPUT=$(basename ${fq}| sed ’s/_pass.fastq//g’);

# Apply adapter trimming
cutadapt -a "AGATCGGAAGAGC"  -m 15 -n 3 ${fq} | fastx_reverse_complement -i > $TRIM_DIR\/${OUTPUT}_trimmed.fastq

# Get reverse complemented sequence
fastx_reverse_complement -i $TRIM_DIR\/${OUTPUT}_trimmed.fastq -o $TRIMRC_DIR\/${OUTPUT}_trimmed_rc.fastq

# Perform genome alignment of reads using bowtie
bowtie -n 2 -k 100 --best --strata -p 10  -x mm10/genome  --un unmapped_out $TRIMRC_DIR\/${OUTPUT}_trimmed_rc.fastq -S > SAM_DIR\/${OUTPUT}.sam

# Convert sam files to bam files
samtools view -bS SAM_DIR\/${OUTPUT}.sam -o BAM_DIR\/${OUTPUT}.bam
samtools sort BAM_DIR\/${OUTPUT}.bam -o BAM_DIR\/${OUTPUT}.sorted.bam
samtools index BAM_DIR\/${OUTPUT}.sorted.bam

# Generate genomeCoverageBed files for the coverage visualisation to decide the flanking regions included in the annotation file
genome_f = mm10.chrom.sizes.txt
genomeCoverageBed -split -strand + -bga -ibam BAM_DIR\/${OUTPUT}.sorted.bam | sort -k 1,1 -k2,2n > BDG_DIR\/${OUTPUT}.bedgraph_fwd
genomeCoverageBed -split -strand - -bga -ibam BAM_DIR\/${OUTPUT}.sorted.bam | sort -k 1,1 -k2,2n > BDG_DIR\/${OUTPUT}.bedgraph_rev

# Generate bigwig files for the check in the genome browser/IGV
wigToBigWig BDG_DIR\/${OUTPUT}.bedgraph_fwd mm10.chrom.sizes.txt BDG_DIR\/${OUTPUT}.fwd.bw
wigToBigWig BDG_DIR\/${OUTPUT}.bedgraph_rev mm10.chrom.sizes.txt BDG_DIR\/${OUTPUT}.rev.bw

# Visualisation of the genome coverage at polyA sites
computeMatrix reference-point -S BDG_DIR\/${OUTPUT}.fwd.bw -R polyA_annotation.fwd.bed -a 200 -b 200 -o CM_DIR\/${OUTPUT}.fwd.ud200.gz --referencePoint center --sortRegions no
plotHeatmap -m  CM_DIR\/${OUTPUT}.fwd.ud200.gz -o HM_DIR\/${OUTPUT}.fwd.ud200.pdf


computeMatrix reference-point -S BDG_DIR\/${OUTPUT}.rev.bw -R polyA_annotation.fwd.bed -a 200 -b 200 -o CM_DIR\/${OUTPUT}.rev.ud200.gz --referencePoint center --sortRegions no
plotHeatmap -m CM_DIR\/${OUTPUT}.rev.ud200.gz -o HM_DIR\/${OUTPUT}.rev.ud200.pdf


done

