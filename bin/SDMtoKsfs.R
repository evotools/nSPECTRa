#!/usr/bin/env Rscript
library(tidyverse)
library(reshape2)
args <- commandArgs(trailingOnly = TRUE)
popu<-args[1]

# Get input file list
temp = list.files(pattern=paste("sdm.",popu,"*", sep=""))

# Edit column headers 
cnames<-c("Chr", "Pos1", "Id1", "Pos2", "Id2", "Dist", "Null", "Cons1", "Cons2", "Anc1", "Anc2", "SameCodon", "StartF", "Mid1F", "Mid2F", "EndF", "StartC", "Mid1C", "Mid2C", "EndC", "Inds")
myfiles = lapply(temp, read_tsv, col_names = cnames)

# Combine all files
dat<-bind_rows(myfiles)
dat$Cons1[which((dat$Cons1 == "stop_gained") | (dat$Cons1 == "stop_lost") | (dat$Cons1 == "start_lost") | (dat$Cons1 == "missense_variant") | (dat$Cons1 == "synonymous_variant"))]<-"coding_variant"
dat$Cons2[which((dat$Cons2 == "stop_gained") | (dat$Cons2 == "stop_lost") | (dat$Cons2 == "start_lost") | (dat$Cons2 == "missense_variant") | (dat$Cons2 == "synonymous_variant"))]<-"coding_variant"

# Keep same non-coding codons
dat <- dat %>% 
  filter(SameCodon == 1) %>% 
  filter(Cons1 != 'coding_variant' & Cons2 != 'coding_variant')

# Fix misplaced 
dat$Rerrange <- unlist(lapply(strsplit(dat$StartC, ''), function(x){ x[2] %in% tolower(LETTERS)}))
dat[dat$Rerrange, c('Mid1C', 'Mid2C')] = dat[dat$Rerrange, c('Mid2C', 'Mid1C')]

# process SDMs
sdms <- dat %>% 
  filter( (Mid1F == 0 & Mid2F != 0) | (Mid1F != 0 & Mid2F == 0) ) %>%
  filter( EndF != 0 )
sdms$Mid = NA
sdms$Mid[sdms$Mid1F != 0] <- sdms$Mid1C[sdms$Mid1F != 0]
sdms$Mid[sdms$Mid2F != 0] <- sdms$Mid2C[sdms$Mid2F != 0]
sdms <- sdms %>% 
  unite('tmp1', c('StartC', 'Mid'), sep = '>') %>%
  unite('mutation', c('tmp1', 'EndC'), sep = '~') %>%
  mutate(mutation = toupper(mutation)) %>%
  select(Chr, Pos1, Pos2, Anc1, Anc2, mutation, EndF)
sdm_Ksfs <- sdms %>%
  select(mutation, EndF, Chr, Pos1, Pos2) %>%
  group_by(mutation, EndF) %>%
  summarise(value = n()) %>%
  pivot_wider(id_cols = EndF, names_from = mutation) %>%
  arrange(EndF) %>%
  rename(sample_frequency = EndF)
sdm_Ksfs[is.na(sdm_Ksfs)] = 0
write_tsv(sdm_Ksfs, paste('sdm_', popu, '.txt', sep = ''))

# process MNPs
mnps <- dat %>% 
  filter( Mid1F == 0 & Mid2F == 0 & EndF != 0 )

mnps <- mnps %>% 
  unite('mutation', c('StartC', 'EndC'), sep = '>') %>%
  mutate(mutation = toupper(mutation)) %>%
  select(Chr, Pos1, Pos2, Anc1, Anc2, mutation, EndF)

mnp_Ksfs <- mnps %>%
  select(mutation, EndF, Chr, Pos1, Pos2) %>%
  group_by(mutation, EndF) %>%
  summarise(value = n()) %>%
  pivot_wider(id_cols = EndF, names_from = mutation) %>%
  arrange(EndF) %>%
  rename(sample_frequency = EndF)
mnp_Ksfs[is.na(mnp_Ksfs)] = 0
write_tsv(mnp_Ksfs, paste('mnp_', popu, '.txt', sep = ''))

# process neighbouring SNPs
neighbouring = which( !((((dat$Mid1F == 0 & dat$Mid2F != 0) | (dat$Mid1F != 0 & dat$Mid2F == 0)) & dat$EndF !=0) | (dat$Mid1F == 0 & dat$Mid2F == 0 & dat$EndF != 0)) )
disnps <- dat[neighbouring, ]

disnps <- disnps %>% 
  unite('mutation1', c('StartC', 'Mid1C'), sep = '>', remove = F) %>%
  unite('mutation2', c('StartC', 'Mid2C'), sep = '>', remove = F) %>%
  unite('mutation3', c('StartC', 'EndC'), sep = '>', remove = F) %>%
  mutate(mutation1 = toupper(mutation1), mutation2 = toupper(mutation2), mutation3 = toupper(mutation3)) %>%
  rename(count1 = Mid1F, count2 = Mid2F, count3 = EndF) %>%
  select(mutation1, count1, mutation2, count2, mutation3, count3)

subs1 = disnps[,c(1,2)] %>% rename(mutation = mutation1, sample_frequency = count1)
subs2 = disnps[,c(3,4)] %>% rename(mutation = mutation2, sample_frequency = count2)
subs3 = disnps[,c(5,6)] %>% rename(mutation = mutation3, sample_frequency = count3)

disnps = rbind(rbind(subs1, subs2), subs3)

disnps_Ksfs <- disnps %>%
  group_by(mutation, sample_frequency) %>%
  summarise(value = n()) %>%
  pivot_wider(id_cols = sample_frequency, names_from = mutation) %>%
  arrange(sample_frequency) 
disnps_Ksfs[is.na(disnps_Ksfs)] = 0
write_tsv(disnps_Ksfs, paste('disnps_', popu, '.txt', sep = ''))


