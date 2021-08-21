#!/usr/bin/env Rscript
# Code to collect all SDMs counts
args = commandArgs(T)
library(tidyverse)
library(ggfortify)
library(ggforce)

inputf <- args[1]

## Get single breed inputs
allfiles <- list.files(pattern=paste(inputf, '*.txt', sep=""), recursive = T)
singlecounts <- allfiles[grepl('doubleCounts\\.', allfiles)]

## Load and collate all data
myfiles<-list()
for(f in singlecounts) {
  fname = str_match(f, "doubleCounts\\.*(.*?)\\.txt")[,2]
  myfiles[[fname]]<-read_tsv(f, col_names = TRUE)
}

## Combine files
dat<-bind_rows(myfiles, .id="Breed")

## Generate counts
counts = dat %>%
  select(Breed, Ind, Codon1, Codon2, Codon3, Count) %>%
  unite('IID', c('Breed', 'Ind'), sep = '-') %>%
  unite('Change', c('Codon1', 'Codon2', 'Codon3'), sep = '>') %>%
  group_by(IID, Change) %>%
  summarise(Count = sum(Count)) %>%
  pivot_wider(values_from = Count, names_from = Change)
counts[is.na(counts)] = 0

# Make PCA
counts[,-1]<-counts[,-1]/rowSums(counts[,-1], na.rm=T)
counts<-counts %>% separate(IID, c("Breed", "Id"), extra = "merge", sep = '-')
pca_res <- prcomp(counts[,-c(1,2)], scale. = TRUE)$x %>% as_tibble() 
pca_res <- cbind(counts[,c(1,2)], pca_res)
pdf("plot_sdm_mutSpectra.pdf", height = 8, width = 12)
ggplot(pca_res, aes(x = PC1, y = PC2, label = Breed, colour = Breed)) +
  geom_mark_ellipse() +
  geom_point()
dev.off()