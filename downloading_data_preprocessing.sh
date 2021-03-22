#! usr/bin/sh

#Title: Downloading data and preprocessing
#Running this file requires mapping.txt
#To run- nohup sh downloading_data_preprocessing.sh

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
seq 10261601 10261606 | parallel prefetch SRR{}
mv SRR*/*.sra . 
#Rename the files according to the mapping information from GEO entry of the dataset. The mapping file should contain two columns, the first one should be SRA Run numbers and the second column should contain sample names obtained from metadata.
awk -F',' 'system("mv " $1 " " $2)' mapping.txt
rm -rf SRR*
  
#Convert SRA to fastq
parallel -j 3 fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-e --clip --origfmt {} ::: Mbnl1KO_Thymus_1 Mbnl1KO_Thymus_2 Mbnl1KO_Thymus_3 WT_Thymus_1 WT_Thymus_2 WT_Thymus_3 

#Renaming fastq files for the paired reads to be identified by R1 and R2
for i in *.fastq.gz; 
  do 
  mv -v "$i" "${i/_pass/R}"; 
done

#*************************Downloading Reference genome*************************# 

#The reference genome and annotations for Mouse (Genome assembly GRCm39) were downloaded from www.ensembl.org.
cd ../
wget -nv -O annotation.gtf.gz http://ftp.ensembl.org/pub/release-103/gtf/mus_musculus/Mus_musculus.GRCm39.103.gtf.gz \
    && gunzip -f annotation.gtf.gz
wget -nv -O genome.fa.gz http://ftp.ensembl.org/pub/release-103/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
&& gunzip -f genome.fa.gz
GTF=$(readlink -f annotation.gtf)
GENOME=$(readlink -f genome.fa)

#*************************PRE-PROCESSING*************************# 
#Quality control using FASTQC
mkdir fastqc_out
parallel "fastqc {} -o fastqc_out" ::: $RAW_DATA/*.fastq.gz


#Read alignment using STAR
#Build STAR index
GDIR=STAR_indices
mkdir $GDIR
STAR --runMode genomeGenerate --genomeFastaFiles $GENOME --sjdbGTFfile $GTF --runThreadN 8 --genomeDir $GDIR

ODIR=results/mapping
mkdir -p $ODIR

#Align reads to genome
for fq1 in $RAW_DATA/*R1.fastq.gz; 
do
fq2=$(echo $fq1 | sed 's/1.fastq.gz/2.fastq.gz/g');
OUTPUT=$(basename ${fq1}| sed 's/R1.fastq.gz//g');

STAR --genomeDir $GDIR \
 --runThreadN 12 \
 --readFilesCommand zcat \
 --readFilesIn  ${fq1} ${fq2}\
 --outFileNamePrefix $ODIR\/${OUTPUT} \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMunmapped Within \
 --outSAMattributes Standard

done

#Rename bam files and move to another directory
for bm in $ODIR/*.bam; 
  do 
    mv -v "$bm" "${bm/_Aligned.sortedByCoord.out.bam/.bam}"; 
done;
mkdir bams
mv $ODIR/*.bam bams/ 

