#! /usr/bin/env Rscript
#To run- Rscript prepare_mm10_annotation.R <GTF_file>
library(rtracklayer)
library(tidyverse)
args = commandArgs(trailingOnly = TRUE)
gff_mm10 <- args[1]
mgg <- import(gff_mm10)

anno = mgg %>% data.frame %>% 
  mutate(seqnames = as.character(seqnames)) %>% 
  dplyr::filter(!grepl("chr[1-19]|X|Y", seqnames)) %>% 
  filter(type %in% "exon") %>% 
  select(seqnames, start,   end ,width, strand, gene_name,  transcript_id)  %>% distinct() %>%
  group_by(gene_name, seqnames, start, end, width, strand) %>% 
  summarise(transcript_ids = paste(transcript_id, sep="",collapse=",")) %>% 
  ungroup()  %>%  
  mutate(number = 1) %>% 
  arrange( seqnames,   start ,    end) %>%  group_by(gene_name, seqnames) %>% #
  mutate(ticker = cumsum(number)) %>% 
  mutate(ExonID = paste(gene_name, ".",seqnames,".", ticker, sep="")) %>%
  dplyr::select(-number, ticker) %>%
  add_count(gene_name) %>%
  data.frame 
colnames(anno) <- c("GeneID", "Chr", "Start", "End", "Width", "Strand", "TranscriptIDs", "Ticker", "ExonID", "n")
save(anno, file=paste("mm10_exon_anno.RData", sep=""))
