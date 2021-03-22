#! usr/bin/sh

#Title: Annotations processing
#Running this file requires mm10.chrom.sizes.txt, polyA_annotation.bed
#To run- nohup sh downloading_data_preprocessing.sh

#*************************Annotation Processing*************************# 

# Download the mm10.chrom.sizes
wget https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# Change che added value based on the genomeCov results
bedtools slop -i polyA_annotation.bed -g mm10.chrom.sizes.txt -l 100 -r 0 -s > flanking100added.polyA_annotation.bed

bedtools slop -i polyA_annotation.bed -g mm10.chrom.sizes.txt -l 60 -r 0 -s > flanking60added.polyA_annotation.bed

