library(tidyverse)
library(readr)
library(readxl)
library(picante)
library(phyloseq)


setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads corrected counts of each taxa in each sample 
#this excludes sequences that were determined to be contamination from nearby wells
seq <- read.table("seqtab_corrected.txt")

numCols <- ncol(seq)

cols <- str_c("col", 1:numCols)

dim(seq)
length(cols)

#gives a number to each sequence (column)
colnames(seq) <- cols

nrow(seq)

#gets list of samples
samples <- rownames(seq)

#makes a variable for the sample in seq instead of the samples just being
#the rownames
seq <- seq %>% mutate(sample = samples)

#there are 3170 sequences without excluding any samples
numCols

nrow(seq)

#this excludes the unrelated samples 
#seq <- seq %>% filter(!str_detect(sample, "Con|Epi|Pcn|Lin"))

#gathers seq into long form
seqG <- seq %>% gather(1:3170, key = "sequence", value = "abundance")


#gets the total number of sequences in each sample
#summ <- seqG %>% group_by(sample) %>% summarize(totalSeqs = sum(abundance))

#nrow(seqG)

#adds number of sequences in sample to each sequence
#seqG <- seqG %>% left_join(summ, by = c("sample"))

#nrow(seqG)

#seqG %>% filter(is.na(totalSeqs))

#gets rid of the rows for the samples that are unrelated (Sean told me which ones to exclude)
seqG <- seqG %>% filter(!str_detect(sample, "Con|Epi|Pcn|Lin"))

seqG %>% filter(str_detect(sample, "Mock")) %>% distinct(sample)

seqG <- seqG %>% filter(!str_detect(sample, "Mock"))

#gets rid of the rows for which a sequence is not present in a sample
#this excludes the sequences that are just present in the unrelated samples (Con|Epi|Pcn|Lin)
seqG <- seqG %>% filter(abundance > 0)

seqG %>% distinct(sequence) %>% nrow()

##I don't want to filter out sequences with relative abundances below .002 for rarefaction

#gets the proportion that each sequence is in each sample
#seqG <- seqG %>% mutate(propSample = abundance/totalSeqs)

#seqG %>% distinct(sequence) %>% nrow() 

#excludes the sequences in samples that have a low relative abundances in the sample
#seqG <- seqG %>% filter(propSample >= 0.002)

#seqG %>% distinct(sequence) %>% nrow()

#gets rid of unnecessary variables
#seqG <- seqG %>% select(-c(propSample, totalSeqs))


#makes sequence names in seqG match sequence names in repSeqs
seqG$sequence <- str_replace(seqG$sequence, "col", "contig_")

seqG %>% distinct(sequence) %>% nrow()

#excludes samples 1223 and JW3
seqG <- seqG %>% filter(sample != "1223")
seqG <- seqG %>% filter(sample != "JW3")


#loads the representative sequences fasta file so that I can get the 
#contig ID assigned to each feature ID
repSeqs <- read.table("representative_sameLengthSeqs.fasta", fill = TRUE)

#gets just the rows in repSeqs that correspond to IDs
repSeqs <- repSeqs %>% filter(str_detect(V1, ">"))

repSeqs %>% distinct(V1) %>% nrow()
repSeqs %>% distinct(V2) %>% nrow()
seqG %>% distinct(sequence) %>% nrow()

nrow(seqG)

#adds contig IDs to feature IDs in seqG
seqG <- seqG %>% left_join(repSeqs, by = c("sequence" = "V2"))

nrow(seqG)

seqG %>% filter(is.na(V1))

#gets rid of ">" in feature IDs
seqG <- seqG %>% mutate(V1 = str_replace(V1, ">", ""))

#loads cleaned up taxonomy
tax <- read.table("taxonomySameLengthSequences_cleanedUp.tsv")

seqG %>% distinct(V1) %>% nrow()
tax %>% distinct(Feature.ID) %>% nrow()

nrow(seqG)

#adds taxonomies to seqG
seqG <- seqG %>% left_join(tax, by = c("V1" = "Feature.ID"))

nrow(seqG)

seqG %>% filter(is.na(shortTaxa))

#ASN9601 corresponds to AS9601
seqG$sample <- ifelse(seqG$sample == "ASN9601", "AS9601", seqG$sample)

#S9503 corresponds to SYN9503
seqG$sample <- ifelse(seqG$sample == "S9503", "SYN9503", seqG$sample)

seqG <- seqG %>% mutate(sample = ifelse(sample == "1213", "1313", sample))
#seqG <- seqG %>% mutate(sample = ifelse(sample == "WH7803", "WH7801", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "SYN1320", "SYN1220", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "9201.2", "9311", sample))

#seqComb$sample <- ifelse(seqComb$sample == "8102", "SYNCLADEXVI", seqComb$sample)
seqG <- seqG %>% mutate(sample = ifelse(sample == "8102", "SYNCLADEXVI", sample))

seqG <- seqG %>% mutate(sample = ifelse(sample == "9201.1", "9201", sample))

seqG %>% filter(str_detect(shortTaxa, "Proch|Synec|proch|synech")) %>% distinct(shortTaxa)

#gets a dataframe without the Pro and Syn sequences
noProSyn <- seqG %>% filter(!str_detect(shortTaxa, "Proch|Synec|proch|synech"))

seqG %>% distinct(sequence) %>% nrow()
noProSyn %>% distinct(sequence) %>% nrow()

seqG %>% distinct(sample) %>% nrow()
noProSyn %>% distinct(sample) %>% nrow()


###rarefy

#seqG %>% group_by(sample) %>% summarize(totalReads = sum(abundance)) %>% arrange(totalReads)
#noProSyn %>% group_by(sample) %>% summarize(totalHetReads = sum(abundance)) %>% arrange(totalHetReads)

#gets rid of unnecessary variables
seqG <- seqG %>% select(sample, V1, abundance)
noProSyn <- noProSyn %>% select(sample, V1, abundance)

#spreads dataframes
seqG <- seqG %>% spread(key = V1, value = abundance, fill = 0)
noProSyn <- noProSyn %>% spread(key = V1, value = abundance, fill = 0)

#the samples in seqG and noProSyn are in the same order
unique(seqG$sample == noProSyn$sample)

#makes the rownames the sample IDs
rownames(seqG) <- seqG$sample
rownames(noProSyn) <- noProSyn$sample

seqG <- seqG %>% select(-sample)
noProSyn <- noProSyn %>% select(-sample)


##from https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/rarefy

#data(BCI)
#S <- specnumber(BCI) # observed number of species
#S
#(raremax <- min(rowSums(BCI)))
#raremax
#Srare <- rarefy(BCI, raremax)
#Srare
#plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
#abline(0, 1)
#arecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)

#rarefies seqG
(raremax_seqG <- min(rowSums(seqG)))
Srare_seqG <- rarefy(seqG, raremax_seqG)

#rarefies noProSyn
(raremax_noProSyn <- min(rowSums(noProSyn)))
Srare_noProSyn <- rarefy(noProSyn, raremax_noProSyn)

#makes rarefaction objects into data frames
Srare_seqG <- Srare_seqG %>% as.data.frame()
Srare_noProSyn <- Srare_noProSyn %>% as.data.frame()

#makes a variable for the sample IDs
Srare_seqG <- Srare_seqG %>% mutate(sample = rownames(Srare_seqG))
Srare_noProSyn <- Srare_noProSyn %>% mutate(sample = rownames(Srare_noProSyn))

colnames(Srare_seqG)[1] <- "Srare_seqG"
colnames(Srare_noProSyn)[1] <- "Srare_noProSyn"

#loads the number of non-Pro, non-Syn features in each sample
#from pcaPlot_andNumASVPlots.R
notRarefied_noProSyn <- read_csv("numNonProNonSynSeqsBySample.csv")

nrow(notRarefied_noProSyn)

#adds rarefied results to number of non-Pro, non-Syn features in each sample
merged <- notRarefied_noProSyn %>% left_join(Srare_seqG, by = c("colsample" = "sample")) %>% 
  left_join(Srare_noProSyn, by = c("colsample" = "sample"))

nrow(merged)
merged %>% filter(is.na(Srare_noProSyn))
merged %>% filter(is.na(Srare_seqG))

str(merged)

#gathers merged dataset
mergedG <- merged %>% gather(numSeqs:Srare_noProSyn, key = "type", value = "richness")

mergedG %>% group_by(colsample) %>% summarize(n = n()) %>% filter(n != 3)

mergedG %>% ggplot(aes(x = colsample, y = richness, fill = type)) + geom_bar(stat = 'identity', position = 'dodge')

merged %>% ggplot(aes(x = numSeqs, y = Srare_noProSyn)) + geom_point() + 
  labs(x = "Number of ASVS, excluding Pro. and Syn.", y = "Rarefied number of ASVs based on number of non-Pro., non-Syn. reads") + 
  annotate("text", x = 5, y = 20, size = 6, label = paste(("R^2="), as.character(round(cor(merged$numSeqs, merged$Srare_noProSyn)^2, 2))))

ggsave("numSeqs_vs_Srare_noProSyn.png", dpi = 600, height = 6, width = 6)
ggsave("numSeqs_vs_Srare_noProSyn.svg", dpi = 600, height = 6, width = 6)


merged %>% ggplot(aes(x = numSeqs, y = Srare_seqG)) + geom_point() + 
  labs(x = "Number of ASVS, excluding Pro. and Syn.", y = "Rarefied number of ASVs based on number of reads") + 
  annotate("text", x = 5, y = 20, size = 6, label = paste(("R^2="), as.character(round(cor(merged$numSeqs, merged$Srare_seqG)^2, 2))))

ggsave("numSeqs_vs_Srare_seqG.png", dpi = 600, height = 6, width = 6)
ggsave("numSeqs_vs_Srare_seqG.svg", dpi = 600, height = 6, width = 6)


merged %>% ggplot(aes(x = Srare_noProSyn, y = Srare_seqG)) + geom_point() + 
  labs(x = "Rarefied number of ASVs based on number of non-Pro., non-Syn. reads",  y = "Rarefied number of ASVs based on number of reads") + 
  annotate("text", x = 5, y = 20, size = 6, label = paste(("R^2="), as.character(round(cor(merged$Srare_noProSyn, merged$Srare_seqG)^2, 2))))

ggsave("Srare_noProSyn_vs_Srare_seqG.png", dpi = 600, height = 6, width = 6)
ggsave("Srare_noProSyn_vs_Srare_seqG.svg", dpi = 600, height = 6, width = 6)

              








                  