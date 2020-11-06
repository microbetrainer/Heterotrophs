library(tidyverse)
library(readr)
library(readxl)
library(picante)
library(phyloseq)
library(ggplot2)


setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads corrected counts of each taxa in each sample 
#this excludes sequences that were determined to be contamination from nearby wells
seq <- read.table("seqtab_corrected.txt")

numCols <- ncol(seq)

cols <- str_c("col", 1:numCols)

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

#now there are 82 samples instead of 96
#nrow(seq)

#gathers seq into long form
seqG <- seq %>% gather(1:3170, key = "sequence", value = "abundance")


#gets the total number of reads in each sample
summ <- seqG %>% group_by(sample) %>% summarize(totalSeqs = sum(abundance))

nrow(seqG)

#adds total number of reads in each sample to seqG
seqG <- seqG %>% left_join(summ, by = c("sample"))

nrow(seqG)

seqG %>% filter(is.na(totalSeqs))

#gets rid of the rows for the samples that are unrelated (Sean told me which ones to exclude)
seqG <- seqG %>% filter(!str_detect(sample, "Con|Epi|Pcn|Lin"))

#gets rid of the rows for the mock samples
seqG <- seqG %>% filter(!str_detect(sample, "Mock"))

#gets rid of the rows for which a sequence is not present in a sample
#this excludes the sequences that are just present in the unrelated samples (Con|Epi|Pcn|Lin) and/or 
#mock samples
seqG <- seqG %>% filter(abundance > 0)

seqG %>% distinct(sequence) %>% nrow()

str(seqG)

#makes a variable for the relative abundance for each sequence in each sample
seqG <- seqG %>% mutate(propSample = abundance/totalSeqs)

seqG %>% distinct(sequence) %>% nrow() 

#excludes the sequences in samples that have relative abundances below 0.002
#this is what I did for the trees and heatmaps
seqG <- seqG %>% filter(propSample >= 0.002)

#this includes Pro and Syn sequences
seqG %>% distinct(sequence) %>% nrow()

#makes sequence names in seqG match sequence names in repSeqs
seqG$sequence <- str_replace(seqG$sequence, "col", "contig_")

seqG %>% distinct(sequence) %>% nrow()


#loads cleaned up taxonomy of the features
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv")

tax %>% nrow()

nrow(seqG)

#adds taxonomy to seqG
seqG <- seqG %>% left_join(tax, by = c("sequence" = "V1"))

nrow(seqG)

seqG %>% filter(is.na(shortTaxa))

#this includes Pro and Syn sequences
seqG %>% distinct(sequence) %>% nrow()

#JW3 only has one sequence which was identified as Pro
seqG %>% filter(sample == "JW3")

#excludes this sample because it only has one sequence
seqG <- seqG %>% filter(sample != "JW3")

#excludes 1223 because it doesn't have Pro or Syn
seqG %>% filter(sample == "1223")
seqG <- seqG %>% filter(sample != "1223")

seqG %>% distinct(sample) %>% nrow()
seqG %>% distinct(sequence) %>% nrow()

seqG %>% filter(abundance == 0)
seqG %>% filter(propSample < 0.002)

#excludes the Pro and Syn sequences 
#they have NA taxonomies because they weren't included in tax or repSeqs
seqG <- seqG %>% filter(!is.na(taxa))

seqG %>% distinct(sample) %>% nrow()
seqG %>% distinct(sequence) %>% nrow()


###extracts each taxonomic level

#kingdom
seqG <- seqG %>% mutate(kingdom = str_extract(Taxon, "k__[a-z,A-Z,0-9]*;{0,1}"))
seqG <- seqG %>% mutate(kingdom = str_replace_all(kingdom, " ", ""))
seqG <- seqG %>% mutate(kingdom = str_replace_all(kingdom, ";", ""))
seqG %>% distinct(kingdom)

#phylum
seqG <- seqG %>% mutate(phylum = str_extract(Taxon, "p__[a-z,A-Z,0-9]*;{0,1}"))
seqG <- seqG %>% mutate(phylum = str_replace_all(phylum, " ", ""))
seqG <- seqG %>% mutate(phylum = str_replace_all(phylum, ";", ""))
seqG %>% distinct(phylum)
seqG %>% filter(is.na(phylum))
seqG$phylum <- ifelse(is.na(seqG$phylum), "p__", seqG$phylum)

#class
seqG <- seqG %>% mutate(class = str_extract(Taxon, "c__[a-z,A-Z,0-9,\\[,\\]]*;{0,1}"))
seqG <- seqG %>% mutate(class = str_replace_all(class, " ", ""))
seqG <- seqG %>% mutate(class = str_replace_all(class, ";", ""))
seqG %>% distinct(class)
seqG %>% filter(is.na(class))
seqG$class <- ifelse(is.na(seqG$class), "c__", seqG$class)

#order
seqG <- seqG %>% mutate(order = str_extract(Taxon, "o__[a-z,A-Z,0-9,\\[,\\]]*;{0,1}"))
seqG <- seqG %>% mutate(order = str_replace_all(order, " ", ""))
seqG <- seqG %>% mutate(order = str_replace_all(order, ";", ""))
seqG %>% distinct(order)
seqG %>% filter(is.na(order)) %>% distinct(Taxon)
seqG %>% filter(order == "o__") %>% distinct(Taxon)
seqG$order <- ifelse(is.na(seqG$order), "o__", seqG$order)

#family
seqG <- seqG %>% mutate(family = str_extract(Taxon, "f__[a-z,A-Z,0-9,\\[,\\]]*;{0,1}"))
seqG <- seqG %>% mutate(family = str_replace_all(family, " ", ""))
seqG <- seqG %>% mutate(family = str_replace_all(family, ";", ""))
seqG %>% distinct(family)
seqG %>% filter(is.na(family)) %>% distinct(Taxon)
seqG %>% filter(family == "f__") %>% distinct(Taxon)
seqG$family <- ifelse(is.na(seqG$family), "f__", seqG$family)

#genus
seqG <- seqG %>% mutate(genus = str_extract(Taxon, "g__[a-z,A-Z,0-9,\\[,\\]]*;{0,1}"))
seqG <- seqG %>% mutate(genus = str_replace_all(genus, " ", ""))
seqG <- seqG %>% mutate(genus = str_replace_all(genus, ";", ""))
seqG %>% distinct(genus)
seqG %>% filter(is.na(genus)) %>% distinct(Taxon)
seqG %>% filter(genus == "g__") %>% distinct(Taxon)
seqG$genus <- ifelse(is.na(seqG$genus), "g__", seqG$genus)

#species
seqG <- seqG %>% mutate(species = str_extract(Taxon, "s__[a-z,A-Z,0-9,\\[,\\]]*;{0,1}"))
seqG <- seqG %>% mutate(species = str_replace_all(species, " ", ""))
seqG <- seqG %>% mutate(species = str_replace_all(species, ";", ""))
seqG %>% distinct(species)
seqG %>% filter(is.na(species)) %>% distinct(Taxon)
seqG %>% filter(species == "s__") %>% distinct(Taxon)
seqG$species <- ifelse(is.na(seqG$species), "s__", seqG$species)


###makes tables of the number of samples each taxonomy is found at
###at each taxonomic level

#all of the samples are either Pro or Syn
seqG %>% filter(!str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916|SYN|WH|9220|S9503|8102|RS9916|0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223|JW7|B7|DV|C12B|JW2|JW4|1213")) %>% 
  distinct(sample) %>% nrow()

#makes a variable for whether sample is Pro or Syn
seqG <- seqG %>% mutate(culture = ifelse(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916"), "Syn", "Pro"))

colnames(seqG)

seqG %>% group_by(culture) %>% distinct(sample) %>% summarize(n = n())
56+18

#kingdom
seqG %>% distinct(kingdom)
seqG %>% group_by(kingdom) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(desc(n))

#phylum
seqG %>% distinct(phylum)
seqG %>% filter(phylum != "p__") %>% group_by(phylum) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(phylum)
seqG %>% filter(phylum != "p__") %>% group_by(phylum) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) %>% 
  write_csv("numSamplesEachPhylumIsFoundIn.csv")


seqG %>% filter(phylum != "p__") %>% group_by(phylum) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples))
seqG %>% filter(phylum == "p__Proteobacteria") %>% distinct(class)
seqG %>% filter(phylum == "p__Proteobacteria") %>% group_by(class) %>% distinct(sample) %>% 
  summarize(n = n()) %>% arrange(desc(n))
seqG %>% filter(phylum != "p__") %>% group_by(phylum, culture) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(phylum)
seqG %>% filter(phylum != "p__") %>% group_by(phylum, culture) %>% distinct(sample) %>% 
  summarize(numSamples = n()) %>% arrange(phylum) %>% mutate(prop = ifelse(culture == "Pro", numSamples/56, numSamples/18))


order <- seqG %>% filter(phylum != "p__") %>% mutate(phylum = str_replace(phylum, "p__", "")) %>% group_by(phylum) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) 
order <- order %>% select(phylum)
order <- as.list(order)
order <- order$phylum
toPlot <- seqG %>% filter(phylum != "p__") %>% mutate(phylum = str_replace(phylum, "p__", ""))
toPlot$phylum <- factor(toPlot$phylum, levels = order)
toPlot <- toPlot %>% group_by(phylum, culture) %>% distinct(sample) %>% summarize(numSamples = n()) 
toPlot <- toPlot %>% spread(key = culture, value = numSamples, fill = 0)
toPlot <- toPlot %>% gather(Pro:Syn, key = "culture", value = "numSamples")
toPlot %>% ggplot(aes(x = phylum, y = numSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Number of samples", x = "Phylum") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("numSamplesEachPhylumIsFoundIn.png", dpi = 600, height = 4, width = 4)
toPlot %>% mutate(propSamples = ifelse(culture == "Pro", numSamples/56, numSamples/18)) %>% ggplot(aes(x = phylum, y = propSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Proportion of samples", x = "Phylum") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("propSamplesEachPhylumIsFoundIn.png", dpi = 600, height = 4, width = 4)
ggsave("propSamplesEachPhylumIsFoundIn.svg", dpi = 600, height = 4, width = 4)

toPlot
toPlot <- toPlot %>% ungroup()
toPlot

#makes dataframe with the prop of pro and syn samples phylums are found in
binom_phylum <- toPlot %>% mutate(propSamples = ifelse(culture == "Pro", numSamples/56, numSamples/18))
binom_phylum

#class
seqG %>% distinct(class)
seqG %>% filter(class != "c__") %>% group_by(class) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(class) 
seqG %>% filter(class != "c__") %>% group_by(class) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) %>% 
  write_csv("numSamplesEachClassIsFoundIn.csv")

seqG %>% filter(class != "c__") %>% group_by(class, culture) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(class)
seqG %>% filter(class != "c__") %>% group_by(class, culture) %>% distinct(sample) %>% 
  summarize(numSamples = n()) %>% arrange(class) %>% mutate(prop = ifelse(culture == "Pro", numSamples/56, numSamples/18))


order <- seqG %>% filter(class != "c__") %>% mutate(class = str_replace(class, "c__", "")) %>% group_by(class) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) 
order <- order %>% select(class)
order <- as.list(order)
order <- order$class
toPlot <- seqG %>% filter(class != "c__") %>% mutate(class = str_replace(class, "c__", ""))
toPlot$class <- factor(toPlot$class, levels = order)
toPlot <- toPlot %>% group_by(class, culture) %>% distinct(sample) %>% summarize(numSamples = n()) 
toPlot <- toPlot %>% spread(key = culture, value = numSamples, fill = 0)
toPlot <- toPlot %>% gather(Pro:Syn, key = "culture", value = "numSamples")
toPlot %>% ggplot(aes(x = class, y = numSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Number of samples", x = "Class") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("numSamplesEachClassIsFoundIn.png", dpi = 600, height = 4, width = 4)
toPlot %>% mutate(propSamples = ifelse(culture == "Pro", numSamples/56, numSamples/18)) %>% ggplot(aes(x = class, y = propSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Proportion of samples", x = "Class") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("propSamplesEachClassIsFoundIn.png", dpi = 600, height = 4, width = 4)
ggsave("propSamplesEachClassIsFoundIn.svg", dpi = 600, height = 4, width = 4)

toPlot
toPlot <- toPlot %>% ungroup()
toPlot

#makes dataframe with the prop of pro and syn samples classes are found in
binom_class <- toPlot %>% mutate(propSamples = ifelse(culture == "Pro", numSamples/56, numSamples/18))
binom_class

#order
seqG %>% distinct(order)
seqG %>% filter(order != "o__") %>% group_by(order) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(order) 
seqG %>% filter(order != "o__") %>% group_by(order) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) %>% 
  write_csv("numSamplesEachOrderIsFoundIn.csv")

order <- seqG %>% filter(order != "o__") %>% mutate(order = str_replace(order, "o__", "")) %>% group_by(order) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) 
order <- order %>% select(order)
order <- as.list(order)
order <- order$order
toPlot <- seqG %>% filter(order != "o__") %>% mutate(order = str_replace(order, "o__", ""))
toPlot$order <- factor(toPlot$order, levels = order)
toPlot <- toPlot %>% group_by(order, culture) %>% distinct(sample) %>% summarize(numSamples = n())
toPlot <- toPlot %>% spread(key = culture, value = numSamples, fill = 0)
toPlot <- toPlot %>% gather(Pro:Syn, key = "culture", value = "numSamples")
toPlot %>% ggplot(aes(x = order, y = numSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Number of samples", x = "Order") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("numSamplesEachOrderIsFoundIn.png", dpi = 600, height = 4, width = 4)
toPlot %>% mutate(propSamples = ifelse(culture == "Pro", numSamples/56, numSamples/18)) %>% ggplot(aes(x = order, y = propSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Proportion of samples", x = "Order") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("propSamplesEachOrderIsFoundIn.png", dpi = 600, height = 4, width = 4)
ggsave("propSamplesEachOrderIsFoundIn.svg", dpi = 600, height = 4, width = 4)

#family
seqG %>% distinct(family)
seqG %>% filter(family != "f__") %>% group_by(family) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(family)
seqG %>% filter(family != "f__") %>% group_by(family) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) %>% 
  write_csv("numSamplesEachFamilyIsFoundIn.csv")

order <- seqG %>% filter(family != "f__") %>% mutate(family = str_replace(family, "f__", "")) %>% group_by(family) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) 
order <- order %>% select(family)
order <- as.list(order)
order <- order$family
toPlot <- seqG %>% filter(family != "f__") %>% mutate(family = str_replace(family, "f__", ""))
toPlot$family <- factor(toPlot$family, levels = order)
toPlot <- toPlot %>% group_by(family, culture) %>% distinct(sample) %>% summarize(numSamples = n())
toPlot <- toPlot %>% spread(key = culture, value = numSamples, fill = 0)
toPlot <- toPlot %>% gather(Pro:Syn, key = "culture", value = "numSamples")
toPlot %>% ggplot(aes(x = family, y = numSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Number of samples", x = "Family") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("numSamplesEachFamilyIsFoundIn.png", dpi = 600, height = 4, width = 4)
toPlot %>% mutate(propSamples = ifelse(culture == "Pro", numSamples/56, numSamples/18)) %>% ggplot(aes(x = family, y = propSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Proportion of samples", x = "Family") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("propSamplesEachFamilyIsFoundIn.png", dpi = 600, height = 4, width = 4)
ggsave("propSamplesEachFamilyIsFoundIn.svg", dpi = 600, height = 4, width = 4)

#genus
seqG %>% distinct(genus)
seqG %>% filter(genus != "g__") %>% group_by(genus) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(genus)
seqG %>% filter(genus != "g__") %>% group_by(genus) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) %>% 
  write_csv("numSamplesEachGenusIsFoundIn.csv")

order <- seqG %>% filter(genus != "g__") %>% mutate(genus = str_replace(genus, "g__", "")) %>% group_by(genus) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) 
order <- order %>% select(genus)
order <- as.list(order)
order <- order$genus
toPlot <- seqG %>% filter(genus != "g__") %>% mutate(genus = str_replace(genus, "g__", ""))
toPlot$genus <- factor(toPlot$genus, levels = order)
toPlot <- toPlot %>% group_by(genus, culture) %>% distinct(sample) %>% summarize(numSamples = n())
toPlot <- toPlot %>% spread(key = culture, value = numSamples, fill = 0)
toPlot <- toPlot %>% gather(Pro:Syn, key = "culture", value = "numSamples")
toPlot %>% ggplot(aes(x = genus, y = numSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Number of samples", x = "Genus") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("numSamplesEachGenusIsFoundIn.png", dpi = 600, height = 4, width = 8)
toPlot %>% mutate(propSamples = ifelse(culture == "Pro", numSamples/56, numSamples/18)) %>% ggplot(aes(x = genus, y = propSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Proportion of samples", x = "Genus") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("propSamplesEachGenusIsFoundIn.png", dpi = 600, height = 4, width = 8)
ggsave("propSamplesEachGenusIsFoundIn.svg", dpi = 600, height = 4, width = 8)


seqG %>% group_by(genus) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) %>% 
  write_csv("numSamplesEachGenusIsFoundIn_includingNAGenus.csv")


#species
#contig_7 is Alt Mac MIT1002 (the lab strain)
seqG <- seqG %>% mutate(species = ifelse(sequence == "contig_7", "macleodii MIT1002", species))
seqG %>% distinct(species)
seqG %>% filter(species != "s__") %>% group_by(species) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(species)
seqG %>% filter(species != "s__") %>% group_by(species) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) %>% 
  write_csv("numSamplesEachSpeciesIsFoundIn.csv")

order <- seqG %>% filter(species != "s__") %>% mutate(species = str_replace(species, "s__", "")) %>% group_by(species) %>% distinct(sample) %>% summarize(numSamples = n()) %>% arrange(desc(numSamples)) 
order <- order %>% select(species)
order <- as.list(order)
order <- order$species
toPlot <- seqG %>% filter(species != "s__") %>% mutate(species = str_replace(species, "s__", ""))
toPlot$species <- factor(toPlot$species, levels = order)
toPlot <- toPlot %>% group_by(species, culture) %>% distinct(sample) %>% summarize(numSamples = n())
toPlot <- toPlot %>% spread(key = culture, value = numSamples, fill = 0)
toPlot <- toPlot %>% gather(Pro:Syn, key = "culture", value = "numSamples")
toPlot %>% ggplot(aes(x = species, y = numSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Number of samples", x = "Species") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("numSamplesEachSpeciesIsFoundIn.png", dpi = 600, height = 4, width = 4)
toPlot %>% mutate(propSamples = ifelse(culture == "Pro", numSamples/56, numSamples/18)) %>% ggplot(aes(x = species, y = propSamples, fill = culture)) + geom_bar(stat = 'identity', position = "dodge") + labs(y = "Proportion of samples", x = "Species") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("propSamplesEachSpeciesIsFoundIn.png", dpi = 600, height = 4, width = 4)
ggsave("propSamplesEachSpeciesIsFoundIn.svg", dpi = 600, height = 4, width = 4)


###runs fisher test to see if proportion samples is different between
###pro and syn samples for each heterotroph phylum and class

#phylum

#spreads number of samples
binom_phylum
binom_phylum <- binom_phylum %>% ungroup() %>% select(-propSamples) %>% spread(key = culture, value = numSamples) 
binom_phylum

#makes a dataframe tbat also has a variable for the fisher test p-value
phylumDF <- binom_phylum %>% mutate(`Fisher test p-value` = NA)
phylumDF

matrix(c(binom_phylum$Pro[1], 56-binom_phylum$Pro[1], binom_phylum$Syn[1], 18-binom_phylum$Syn[1]), 2,2)

#https://stats.stackexchange.com/questions/28232/fisher-test-in-r
#runs fisher test to see if proportion samples is different between
#pro and syn samples for each heterotroph phylum
for (i in 1:nrow(phylumDF)) {
  phylumDF$`Fisher test p-value`[i] <- fisher.test(matrix(c(binom_phylum$Pro[i], 56-binom_phylum$Pro[i], binom_phylum$Syn[i], 18-binom_phylum$Syn[i]), 2,2), alternative = "two.sided")$p.value
}

#makes a variable for the proportion of pro and syn samples that each heterotroph
#phylum is found in
phylumDF
phylumDF <- phylumDF %>% mutate(`Prop. of Pro. samples` = Pro/56)
phylumDF <- phylumDF %>% mutate(`Prop. of Syn. samples` = Syn/18)
phylumDF

#renames variables
colnames(phylumDF)[2:3]
colnames(phylumDF)[2:3] <- c("Num. of Pro. samples", "Num. of Syn. samples")
phylumDF

#reorder columns 
phylumDF <- phylumDF %>% select(1:3, 5:6, 4)
phylumDF



#class

#spreads number of samples
binom_class <- binom_class %>% ungroup() %>% select(-propSamples) %>% spread(key = culture, value = numSamples)
binom_class

#makes a dataframe tbat also has a variable for the fisher test p-value
classDF <- binom_class %>% mutate(`Fisher test p-value` = NA)

matrix(c(binom_class$Pro[1], 56-binom_class$Pro[1], binom_class$Syn[1], 18-binom_class$Syn[1]), 2,2)

#https://stats.stackexchange.com/questions/28232/fisher-test-in-r
#runs fisher test to see if proportion samples is different between
#pro and syn samples for each heterotroph class
for (i in 1:nrow(classDF)) {
  classDF$`Fisher test p-value`[i] <- fisher.test(matrix(c(binom_class$Pro[i], 56-binom_class$Pro[i], binom_class$Syn[i], 18-binom_class$Syn[i]), 2,2), alternative = "two.sided")$p.value
}

#makes a variable for the proportion of pro and syn samples that each heterotroph
#class is found in
classDF
classDF <- classDF %>% mutate(`Prop. of Pro. samples` = Pro/56)
classDF <- classDF %>% mutate(`Prop. of Syn. samples` = Syn/18)
classDF

#renames variables
colnames(classDF)[2:3]
colnames(classDF)[2:3] <- c("Num. of Pro. samples", "Num. of Syn. samples")
classDF

#reorder columns 
classDF <- classDF %>% select(1:3, 5:6, 4)
classDF

#calculated adjusted p-value 
#from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099145/
phylumDF$`Adjusted p-value` <- p.adjust(phylumDF$`Fisher test p-value`, method="fdr")
classDF$`Adjusted p-value` <- p.adjust(classDF$`Fisher test p-value`, method="fdr")

write_csv(phylumDF, "fisherTestPropOfSamplesProSyn_byPhylum.csv")
write_csv(classDF, "fisherTestPropOfSamplesProSyn_byClass.csv")



