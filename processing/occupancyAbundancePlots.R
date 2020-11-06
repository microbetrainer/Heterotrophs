library(tidyverse)
library(readr)
library(ggplot2)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads proportion of samples of each culture each sequence is found in
#this came from ubiquityAbundance_withoutSample1223.R (checked)
#not including Pro, Syn
ubi <- read_csv("ubiquity_noContamination_noMock_noProSyn_proportionByCultureType_withoutSample1223.csv")

ubi %>% nrow()
ubi %>% distinct(sequence) %>% nrow()
ubi %>% group_by(sequence) %>% summarize(n = n()) %>% filter(n != 2)

#loads relative abundance of each sequence in each sample
#from randomForestProSyn.R (checked)
#not including Pro, Syn
seqG <- read_csv("seqG_fromRandomForestProSyn.csv")

head(ubi)
head(seqG)

ubi %>% distinct(sequence) %>% nrow()
seqG %>% distinct(sequence) %>% nrow()

seqG %>% distinct(sample) %>% nrow()

#gets rid of variables so I can spread data
seqG <- seqG %>% select(sample, sequence, propSample, culture)

head(seqG)

#spreads seqG so that the sequences that are not in a sample will have 0 abundance for that sample
seqSpread <- seqG %>% spread(key = sequence, value = propSample, fill = 0)

dim(seqSpread)

#gathers seqG again, now with rows for sequences that are not present in a sample
seqG <- seqSpread %>% gather(3:237, key = "V1", value = "propSample")

head(seqG)

seqG %>% distinct(sample) %>% nrow()
seqG %>% distinct(V1) %>% nrow()
seqG %>% group_by(V1) %>% distinct(sample) %>% summarize(n = n()) %>% filter(n != 74)
seqG %>% group_by(sample) %>% distinct(V1) %>% summarize(n = n()) %>% filter(n != 235)
seqG %>% group_by(sample) %>% distinct(culture) %>% summarize(n = n()) %>% filter(n != 1)

#some of these are low because Pro and Syn have been excluded
seqG %>% group_by(sample) %>% summarize(sumProp = sum(propSample)) %>% arrange(desc(sumProp))
seqG %>% group_by(sample) %>% summarize(sumProp = sum(propSample)) %>% arrange(sumProp)

#I am eventually going to take log of the prop so I replaced prop 0 with prop .002
seqG %>% head()
seqG <- seqG %>% mutate(propSample = ifelse(propSample == 0, .002, propSample))
seqG %>% filter(propSample == 0)

head(seqG)
seqGSummaryNotSplitByCulture <- seqG %>% group_by(V1) %>% summarize(meanProp = mean(log10(propSample)), n = n())
seqGSummaryNotSplitByCulture %>% head()
seqGSummaryNotSplitByCulture %>% filter(n != 74)

head(seqG)

#calculates the mean relative abundance across samples of each culture for each sequence
seqGSummary <- seqG %>% group_by(culture, V1) %>% summarize(meanProp = mean(log10(propSample)), n = n())

head(seqGSummary)

seqGSummary <- seqGSummary %>% ungroup()

seqGSummary %>% distinct(culture, n)
seqGSummary %>% group_by(V1) %>% summarize(n = n()) %>% filter(n != 2)

head(ubi)

nrow(ubi)
ubi %>% group_by(sequence) %>% summarize(n = n()) %>% filter(n != 2)

head(seqGSummary)
head(ubi)

nrow(ubi)
nrow(seqGSummary)

#adds ubi to seqGSummary
seqGSummary <- seqGSummary %>% left_join(ubi, by = c("V1" = "sequence", "culture"))

nrow(seqGSummary)
head(seqGSummary)

seqGSummary %>% filter(is.na(proportionOfSamplesOfCulture))

seqGSummary$culture <- str_replace(seqGSummary$culture, "Pro", "Across Pro. samples")
seqGSummary$culture <- str_replace(seqGSummary$culture, "Syn", "Across Syn. samples")

head(seqGSummary)


seqGSummary %>% filter(meanProp == 0) %>% distinct(proportionOfSamplesOfCulture)
seqGSummary %>% filter(proportionOfSamplesOfCulture == 0) %>% distinct(meanProp)

#exclude rows for sequences that are not present in any of samples of culture 
#because this messes up plot
seqGSummary <- seqGSummary %>% filter(proportionOfSamplesOfCulture != 0)
seqGSummary %>% nrow()

#how does a positive relationship show selection?? 
#if a sequence has high relative abundance, doesn't that mean that it will be more likely to contaminate/establish in other cultures?? 
seqGSummary %>% ggplot(aes(x = meanProp, y = proportionOfSamplesOfCulture, color = culture)) + 
  geom_point() + labs(x = "Mean log10(relative abundance in sample)", 
  y = "Proportion of samples ASV was found in", color = "") 

ggsave("occupancyAbundance_ASV.png", dpi = 600, height = 13, width = 10)
ggsave("occupancyAbundance_ASV.svg", dpi = 600, height = 13, width = 10)



#loads proportion of samples each sequence (not split by culture) 
#is found in this came from ubiquityAbundance_withoutSample1223.R (checked)
#not including Pro, Syn
ubiNotByCulture <- read_csv("ubiquity_noContamination_noMock_noProSyn_proportionAscrossAllSamples_withoutSample1223.csv")

nrow(seqGSummaryNotSplitByCulture)
seqGSummaryNotSplitByCulture <- seqGSummaryNotSplitByCulture %>% left_join(ubiNotByCulture, by = c("V1" = "sequence"))
nrow(seqGSummaryNotSplitByCulture)

seqGSummaryNotSplitByCulture %>% filter(is.na(numSamples))

seqGSummaryNotSplitByCulture %>% filter(numSamples == 0)

seqGSummaryNotSplitByCulture %>% ggplot(aes(x = meanProp, y = numSamples/74)) + 
  geom_point() + labs(x = "Mean log10(relative abundance in sample)", 
  y = "Proportion of samples ASV was found in", color = "") 

ggsave("occupancyAbundance_ASV_notSplitByCulture.png", dpi = 600, height = 13, width = 10)
ggsave("occupancyAbundance_ASV_notSplitByCulture.svg", dpi = 600, height = 13, width = 10)




###clusters


#loads biom tsv of clustered sean and seawater samples
#includes Pro and Syn
#feature-table-fromBiom-sequences_noContamination_noMock_withProSyn_no1223_plusSeawater_dn-1Edited.tsv (checked)
biom <- read.table("feature-table-fromBiom-sequences_noContamination_noMock_withProSyn_no1223_plusSeawater_dn-1Edited.tsv", colClasses = "character")

#makes column names the first row of biom
colnames(biom) <- biom[1,]
biom <- biom[-1,]  

ncol(biom)
colnames(biom)[1:3]
colnames(biom)[2060:2063]

#gathers the biom dataset into long format
biomG <- biom %>% gather(2:2063, key = "sample", value = "present")

#there are 759 clustered features and 2062 sequences/samples
biomG %>% distinct(OTU_ID) %>% nrow()
biomG %>% distinct(sample) %>% nrow()

str(biomG)

#makes variable indicating whether a feature is present in a sample numeric
biomG$present <- as.numeric(biomG$present)

str(biomG)

biomG %>% filter(present > 0) %>% distinct(OTU_ID) %>% nrow()
biomG %>% filter(present > 0) %>% distinct(sample) %>% nrow()

biomG %>% filter(present > 1)

#gets just the rows for which a clustered feature is present in a sample
biomG <- biomG %>% filter(present > 0)

biomG %>% group_by(OTU_ID) %>% summarize(n = n()) %>% filter(n > 1)
biomG %>% group_by(OTU_ID) %>% distinct(sample) %>% summarize(n = n()) %>% filter(n > 1)


#seqtab_corrected.txt (checked) came from CC16SAna_edited.R (checked)
#loads corrected counts of each taxa in each sample 
#this excludes sequences that were determined to be contamination from nearby wells
seq <- read.table("seqtab_corrected.txt")

numCols <- ncol(seq)

cols <- str_c("col", 1:numCols)

#gives a number to each sequence (column)
colnames(seq) <- cols

#there are 96 samples to begin with
nrow(seq)

#gets list of samples
samples <- rownames(seq)

#makes a variable for the sample in seq instead of the samples just being
#the rownames
seq <- seq %>% mutate(sample = samples)

#there are 3170 sequences without excluding any samples
numCols

dim(seq)
colnames(seq)[3170:3171]

#gathers seq into long form
seqG <- seq %>% gather(1:3170, key = "sequence", value = "abundance")

str(seqG)

#gets rid of the rows for which a sequence is not present in a sample
seqG <- seqG %>% filter(abundance > 0)

#excludes unrelated samples
seqG <- seqG %>% filter(!str_detect(sample, "Epi|Pcn|Lin"))

#excludes Mock samples
seqG <- seqG %>% filter(!str_detect(sample, "Mock"))

##excludes samples JW3 and 1223 
seqG <- seqG %>% filter(sample != "JW3")
seqG <- seqG %>% filter(sample != "1223")

#excludes seawater sequences
seqG %>% filter(str_detect(sample, "Con")) %>% distinct(sample)
seqG <- seqG %>% filter(!str_detect(sample, "Con"))

seqG %>% distinct(sample) %>% nrow()

head(seqG)
seqG

#gets the total number of sequences in each sample
summ <- seqG %>% group_by(sample) %>% summarize(totalSeqs = sum(abundance))

nrow(summ)

nrow(seqG)

#adds number of sequences in sample to each sequence
seqG <- seqG %>% left_join(summ, by = c("sample"))

nrow(seqG)

seqG %>% filter(is.na(totalSeqs))

#calculates the relative abundance of each sequence in each sample
seqG <- seqG %>% mutate(prop = abundance/totalSeqs)

seqG %>% group_by(sample) %>% summarize(totalProp = sum(prop)) %>% filter(totalProp != 1)

#excludes the sequences in a sample that have less than .002 relative abundance
seqG <- seqG %>% filter(prop >= 0.002) 

#includes Pro and Syn
seqG %>% distinct(sequence) %>% nrow()

head(seqG)

#makes sequence IDs match biomG
seqG$sequence <- str_replace(seqG$sequence, "col", "")
seqG$sequence <- str_c(seqG$sequence, "seq")

head(seqG)

#gets rid of unnecessary variables
#seqG <- seqG %>% select(-totalSeqs, -prop)

biomG
head(biomG)
nrow(biomG)

#adds abundance of each sequence in each sample to biomG (the clustering results)
biomG <- biomG %>% left_join(seqG, by = c("sample" = "sequence"))

nrow(biomG)

head(biomG)

#excludes sequences that were from seawater samples
biomG <- biomG %>% filter(!is.na(prop))

#includes Pro and Syn
biomG %>% distinct(sample) %>% nrow()

#gets the total number of reads in each sample
biomG
totalReadsBySamples <- biomG %>% distinct(sample.y, totalSeqs)
head(totalReadsBySamples)
nrow(totalReadsBySamples)

head(biomG)
biomG %>% group_by(sample.y, OTU_ID, sample) %>% summarize(n = n()) %>% filter(n != 1)


#loads taxonomy of clusters which came from analysisSeanAndSeawaterClustered.txt (checked)
tax <- read_tsv("sequences_noContamination_noMock_withProSyn_no1223_plusSeawaterClusteredTaxonomy.tsv")

#gets rid of the first row which just has variable names
tax <- tax[-1,]

#correct number of clusters
nrow(tax)

head(biomG)
head(tax)

nrow(biomG)

#adds taxonomies to seawater
biomG <- biomG %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))

nrow(biomG)

head(biomG)
biomG %>% filter(is.na(Taxon))

#excludes Pro and Syn clusters
biomG %>% filter(str_detect(Taxon, "Proch|Synech")) %>% distinct(Taxon)
biomG <- biomG %>% filter(!str_detect(Taxon, "Proch|Synech"))

#correct number of sequences
biomG %>% distinct(sample) %>% nrow()

biomG %>% distinct(sample.y) %>% nrow()

biomG
head(biomG)
#makes a variable for whether sample is a Pro or Syn culture
biomG <- biomG %>% mutate(culture = ifelse(str_detect(sample.y, "SYN|WH|9220|S9503|8102|RS9916"), "Synechococcus", "Prochlorococcus"))

#correct
biomG %>% group_by(culture) %>% distinct(sample.y) %>% summarize(n = n())

head(biomG)
biomG %>% nrow()
biomG %>% distinct(OTU_ID, sample, sample.y) %>% nrow()

head(biomG)

#gets the total abundance of each cluster in each sample
biomG
abun <- biomG %>% group_by(sample.y, OTU_ID) %>% summarize(sumAbundance = sum(abundance))

head(abun)
abun <- abun %>% ungroup()
head(abun)

head(biomG)
biomG %>% filter(present != 1)
biomG %>% arrange(abundance) %>% select(abundance) %>% slice(1:10)
biomG %>% arrange(abundance) %>% select(prop) %>% slice(1:10)

head(abun)

nrow(abun)
abun <- abun %>% left_join(totalReadsBySamples, by = c("sample.y"))
nrow(abun)

head(abun)
abun %>% filter(is.na(totalSeqs))

abun %>% distinct(sample.y, totalSeqs)

#number of clusters each sample has
abun %>% group_by(sample.y) %>% summarize(n = n()) %>% arrange(desc(n))

head(abun)

#calculates the relative abundance of each cluster in each sample
abun <- abun %>% mutate(propSample = sumAbundance/totalSeqs)

head(abun)

#some of the relative abundances are so low because Pro and Syn 
#are excluded
abun %>% group_by(sample.y) %>% summarize(sumProp = sum(propSample)) %>% arrange(sumProp)

abun %>% group_by(sample.y) %>% summarize(sumProp = sum(propSample)) %>% arrange(desc(sumProp))

head(abun)

#number of clusters each sample has
abun %>% group_by(sample.y) %>% summarize(n = n())

head(abun)

#there are 119 clusters total
abun %>% ungroup() %>% distinct(OTU_ID) %>% summarize(n = n())

abun <- abun %>% ungroup() 

head(abun)

#gets rid of unnecessary variables
abun <- abun %>% select(-c(sumAbundance, totalSeqs))
abun %>% head()

#spreads data so there is a relative abundance for each cluster in each sample
abun <- abun %>% spread(key = sample.y, value = propSample, fill = 0)

head(abun)
dim(abun)

#gathers data
abun <- abun %>% gather(2:75, key = "sample", value = "propSample")

head(abun)

abun %>% group_by(sample) %>% summarize(n = n()) %>% filter(n != 119)
abun %>% group_by(OTU_ID) %>% summarize(n = n()) %>% filter(n != 74)

head(abun)

#makes a variable for whether sample is a Pro or Syn culture
abun <- abun %>% mutate(culture = ifelse(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916"), "Synechococcus", "Prochlorococcus"))

head(abun)

#correct
abun %>% group_by(culture) %>% distinct(sample) %>% summarize(n = n())

abun %>% nrow()
abun %>% distinct(OTU_ID, sample) %>% nrow()

head(abun)
#each OTU is in 74 samples
abun %>% group_by(OTU_ID) %>% summarize(n = n()) %>% filter(n != 74)

#each sample has 119 clusters
abun %>% group_by(sample) %>% summarize(n = n()) %>% filter(n != 119)

head(abun)
abun %>% group_by(culture) %>% distinct(OTU_ID) %>% summarize(n = n())

abun %>% filter(propSample == 0)
abun %>% filter(propSample != 0) %>% arrange(propSample)

abun %>% group_by(OTU_ID, sample) %>% summarize(n = n()) %>% filter(n != 1)

#replaces 0 relative abundance of clusters with .002 so that I can take log10 later
abun <- abun %>% mutate(propSample = ifelse(propSample == 0, .002, propSample))

abun %>% group_by(OTU_ID) %>% summarize(n = n()) %>% filter(n != 74)
abun %>% group_by(sample) %>% summarize(n = n()) %>% filter(n != 119)

#calculates the mean log10(relative abundance) of each cluster across all the samples
abun
abunSummaryNotGroupedByCulture <- abun %>% group_by(OTU_ID) %>% summarize(meanProp = mean(log10(propSample)), n = n())
head(abunSummaryNotGroupedByCulture)
abunSummaryNotGroupedByCulture %>% filter(n != 74)

#calculates the mean log10(relative abundance) of each cluster across Pro and Syn cultures separately
head(abun)
abunSummary <- abun %>% group_by(culture, OTU_ID) %>% summarize(meanProp = mean(log10(propSample)), n = n())

head(abunSummary)
abunSummary <- abunSummary %>% ungroup()
abunSummary %>% group_by(culture) %>% summarize(n = n())
abunSummary %>% distinct(culture, n)

##calculates the number of samples of each culture type that each 
##sequence is found in

head(biomG)
biomG %>% distinct(present)

head(biomG)
biomG %>% nrow()
biomG %>% distinct(OTU_ID, sample, sample.y) %>% nrow()

biomG %>% arrange(prop) %>% slice(1)

head(biomG)
biomG <- biomG %>% ungroup()
ubi <- biomG %>% group_by(OTU_ID, culture) %>% distinct(sample.y) %>% summarize(numCultures = n())
head(ubi)

biomG <- biomG %>% ungroup()
ubiNotGroupedByCulture <- biomG %>% group_by(OTU_ID) %>% distinct(sample.y) %>% summarize(numCultures = n())

ubi <- ubi %>% ungroup()

#calculates the proportion of samples of each culture type that each sequence is found in  
ubi <- ubi %>% mutate(propCultures = ifelse(culture == "Prochlorococcus", numCultures/56, numCultures/18))

head(ubi)
head(abunSummary)

ubi %>% ungroup() %>% distinct(OTU_ID) %>% nrow()
abunSummary %>% ungroup() %>% distinct(OTU_ID) %>% nrow()

ubi <- ubi %>% ungroup() 
abunSummary <- abunSummary %>% ungroup()

#gets rid of unnecessary variables
head(ubi)
ubi <- ubi %>% select(-numCultures)

#spreads ubi so there is a row for each cluster and culture type
head(ubi)
ubi <- ubi %>% spread(key = culture, value = propCultures, fill = 0)
head(ubi)
nrow(ubi)
ubi <- ubi %>% gather(Prochlorococcus:Synechococcus, key = "culture", value = "propCultures")
head(ubi)
nrow(ubi)

nrow(abunSummary)
nrow(ubi)

head(abunSummary)
head(ubi)

abunSummary %>% group_by(OTU_ID) %>% summarize(n = n()) %>% filter(n != 2)
ubi %>% group_by(OTU_ID) %>% summarize(n = n()) %>% filter(n != 2)

head(abunSummary)
head(ubi)

nrow(abunSummary)
abunSummary <- abunSummary %>% left_join(ubi, by = c("OTU_ID", "culture"))
nrow(abunSummary)

head(abunSummary)
abunSummary %>% filter(is.na(propCultures))

abunSummary
head(abunSummary)

abunSummary$culture <- str_replace(abunSummary$culture, "Prochlorococcus", "Across Pro. samples")
abunSummary$culture <- str_replace(abunSummary$culture, "Synechococcus", "Across Syn. samples")

head(abunSummary)
abunSummary %>% distinct(culture, n)
abunSummary %>% nrow()

abunSummary %>% filter(propCultures == 0) %>% distinct(meanProp)
abunSummary %>% filter(meanProp == 0) %>% distinct(propCultures)

#exclude rows for sequences that are not present in any of samples of culture 
#because this messes up plot
abunSummary <- abunSummary %>% filter(propCultures != 0)
abunSummary %>% filter(propCultures == 0)

head(abunSummary)
nrow(abunSummary)

#how does a positive relationship show selection?? 
#if an OTU has high relative abundance, doesn't that mean that it will be more likely to contaminate/establish in other cultures?? 
abunSummary %>% ggplot(aes(x = meanProp, y = propCultures, color = culture)) + 
  geom_point() + labs(x = "Mean log10(relative abundance in sample)", 
  y = "Proportion of samples OTU was found in", color = "")

ggsave("occupancyAbundance_OTU.png", dpi = 600, height = 13, width = 10)
ggsave("occupancyAbundance_OTU.svg", dpi = 600, height = 13, width = 10)


abunSummaryNotGroupedByCulture %>% nrow()
ubiNotGroupedByCulture %>% nrow()

#adds number of samples clusters are found to mean log10(relative abundance)
#of clusters data
abunSummaryNotGroupedByCulture
ubiNotGroupedByCulture
abunSummaryNotGroupedByCulture <- abunSummaryNotGroupedByCulture %>% left_join(ubiNotGroupedByCulture, by = c("OTU_ID"))

abunSummaryNotGroupedByCulture %>% nrow()

head(abunSummaryNotGroupedByCulture)
abunSummaryNotGroupedByCulture %>% filter(is.na(numCultures))

head(abunSummaryNotGroupedByCulture)

abunSummaryNotGroupedByCulture %>% filter(numCultures == 0)
abunSummaryNotGroupedByCulture %>% arrange(numCultures)

abunSummaryNotGroupedByCulture

abunSummaryNotGroupedByCulture %>% ggplot(aes(x = meanProp, y = numCultures/74)) + 
  geom_point() + labs(x = "Mean log10(relative abundance in sample)", 
                      y = "Proportion of samples OTU was found in", color = "")

ggsave("occupancyAbundance_OTU_notSplitByCulture.png", dpi = 600, height = 13, width = 10)
ggsave("occupancyAbundance_OTU_notSplitByCulture.svg", dpi = 600, height = 13, width = 10)


