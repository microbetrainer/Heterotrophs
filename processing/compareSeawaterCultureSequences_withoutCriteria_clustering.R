library(tidyverse)
library(readr)
library(readxl)
library(picante)
library(phyloseq)
library(svglite)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#feature-table-fromBiom-sequences_noContamination_noMock_withProSyn_no1223_plusSeawater_dn-1.tsv came 
#from analysisSeanAndSeawaterClustered.txt (checked)
#feature-table-fromBiom-sequences_noContamination_noMock_withProSyn_no1223_plusSeawater_dn-1Edited.tsv is identical to 
#feature-table-fromBiom-sequences_noContamination_noMock_withProSyn_no1223_plusSeawater_dn-1.tsv except I uncommented the comments
unedited <- read.table("feature-table-fromBiom-sequences_noContamination_noMock_withProSyn_no1223_plusSeawater_dn-1.tsv", colClasses = "character")
edited <- read.table("feature-table-fromBiom-sequences_noContamination_noMock_withProSyn_no1223_plusSeawater_dn-1Edited.tsv", colClasses = "character")
edited <- edited[-1,]
unique(edited$V1 == unedited$V1)

unedited <- unedited %>% select(-V1)
edited <- edited %>% select(-V1)

edited <- lapply(edited, function(x) as.numeric(x))
unedited <- lapply(unedited, function(x) as.numeric(x))

identical(edited, unedited)


#feature-table-fromBiom-sequences_noContamination_noMock_withProSyn_no1223_plusSeawater_dn-1Edited.tsv (checked)
biom <- read.table("feature-table-fromBiom-sequences_noContamination_noMock_withProSyn_no1223_plusSeawater_dn-1Edited.tsv", colClasses = "character")

#makes column names the first row of biom
colnames(biom) <- biom[1,]
biom <- biom[-1,]  

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

#gathers seq into long form
seqG <- seq %>% gather(1:3170, key = "sequence", value = "abundance")

#gets rid of the rows for which a sequence is not present in a sample
seqG <- seqG %>% filter(abundance > 0)

#excludes unrelated samples, doesn't exlude seawater samples ("Con")
seqG <- seqG %>% filter(!str_detect(sample, "Epi|Pcn|Lin"))

#excludes Mock samples
seqG <- seqG %>% filter(!str_detect(sample, "Mock"))

##excludes samples JW3 and 1223 
seqG <- seqG %>% filter(sample != "JW3")
seqG <- seqG %>% filter(sample != "1223")

#includes just the seawater samples that weren't incubated (T16)
#I didn't select ConCt0 because it had under 10,000 reads and Sean 
#told me to exclude it if it had under 10,000 reads

seqG %>% filter(str_detect(sample, "ConCt0")) %>% distinct(sample)

seqG <- seqG %>% filter(!str_detect(sample, "ConCt0"))

seqG %>% filter(str_detect(sample, "Con")) %>% distinct(sample)
seqG %>% filter(str_detect(sample, "t16")) %>% distinct(sample)

seqG <- seqG %>% filter(!str_detect(sample, "t16"))

seqG %>% distinct(sample) %>% nrow()
seqG %>% filter(!str_detect(sample, "Con")) %>% distinct(sample) %>% nrow()
seqG %>% filter(str_detect(sample, "Con")) %>% distinct(sample)


head(seqG)

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

seqG %>% filter(str_detect(sample, "Con")) %>% distinct(sample)

#excludes the sequences in a sample that have less than .002 relative abundance for sequences 
#in cultures but not sequences in seawater
seqG <- seqG %>% filter(prop >= 0.002 | str_detect(sample, "Con")) 

seqG %>% filter(prop < 0.002) %>% distinct(sample)

#this is the correct number of sequences, including Pro and Syn sequences, for the 
#culture samples
seqG %>% filter(!str_detect(sample, "Con")) %>% distinct(sequence) %>% nrow()

#makes sequence IDs match biomG
seqG$sequence <- str_replace(seqG$sequence, "col", "")
seqG$sequence <- str_c(seqG$sequence, "seq")

#gets rid of unnecessary variables
seqG <- seqG %>% select(-totalSeqs, -prop)

nrow(biomG)

#adds abundance of each sequence in each sample to biomG (the clustering results)
biomG <- biomG %>% left_join(seqG, by = c("sample" = "sequence"))

nrow(biomG)

biomG %>% filter(is.na(abundance))

#gets just the seawater samples
seawater <- biomG %>% filter(str_detect(sample.y, "Con"))

seawater %>% distinct(sample.y)

#gets the total abundance of each cluster across the two seawater samples
seawater <- seawater %>% group_by(OTU_ID) %>% summarize(totalAcrossSamples = sum(abundance), n = n())

seawater %>% nrow()
seawater %>% distinct(OTU_ID) %>% nrow()

#get the total number of sequences in the two seawater samples
seawaterSumm <- summ %>% filter(str_detect(sample, "Con"))

seawaterSumm %>% distinct(sample)

#calculates the sum of the total number of sequences across the two seawater samples
seawaterSumm <- seawaterSumm %>% summarize(totalSeawaterSeqs = sum(totalSeqs))

nrow(seawater)

#makes a variable in seawater for the sum of the total number of sequences across the two seawater samples
seawater <- seawater %>% mutate(totalSeawaterSeqs = seawaterSumm$totalSeawaterSeqs)

nrow(seawater)

#makes a variable for the relative abundance of each sequence in the seawater samples
seawater <- seawater %>% mutate(prop = totalAcrossSamples/totalSeawaterSeqs)

seawater %>% summarize(totalProp = sum(prop))

#gets just the culture sequences
culture <- biomG %>% filter(!str_detect(sample.y, "Con"))

#correct number of sequences and samples
culture %>% distinct(sample) %>% summarize(n = n())
culture %>% distinct(sample.y) %>% summarize(n = n())

nrow(culture)

#27 clusters are in culture samples and also found in seawater samples
culture %>% semi_join(seawater, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()
seawater %>% semi_join(culture, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

#gets the sequences that are in culture samples and also found in seawater samples
cultureSeawater <- culture %>% semi_join(seawater, by = c("OTU_ID"))

cultureSeawater %>% nrow()
cultureSeawater %>% distinct(OTU_ID) %>% nrow()

#gets the distinct clusters that are in both culture and seawater samples
cultureSeawater <- cultureSeawater %>% distinct(OTU_ID)

cultureSeawater %>% nrow()

#makes variable for whether sequence in seawater is also present in at least one culture sample
seawater <- seawater %>% mutate(presentInCulture = ifelse(OTU_ID %in% cultureSeawater$OTU_ID, TRUE, FALSE))

seawater %>% filter(presentInCulture == TRUE) %>% nrow()

#loads taxonomy which came from analysisSeanAndSeawaterClustered.txt
tax <- read_tsv("sequences_noContamination_noMock_withProSyn_no1223_plusSeawaterClusteredTaxonomy.tsv")

#gets rid of the first row which just has variable names
tax <- tax[-1,]

#correct number of clusters
nrow(tax)

nrow(seawater)

#adds taxonomies to seawater
seawater <- seawater %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))

nrow(seawater)

seawater %>% filter(is.na(Taxon))


#these are the taxons that will be labeled on the plot
seawater %>% filter(presentInCulture == TRUE) %>% distinct(Taxon)

#makes new cleaned up taxon variable
seawater <- seawater %>% mutate(newTaxon = Taxon)

seawater %>% filter(presentInCulture == TRUE) %>% distinct(Taxon, newTaxon)

#gets rid "k__", "p__"... in newTaxon
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace(newTaxon, "k__", ""), NA))
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace(newTaxon, "; p__", " "), NA))
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace(newTaxon, "; c__", " "), NA))
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace(newTaxon, "; o__", " "), NA))
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace(newTaxon, "; f__", " "), NA))
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace(newTaxon, "; g__", " "), NA))
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace(newTaxon, "; s__", " "), NA))

#sequences that are not also in culture sample now have NA newTaxon values

seawater %>% filter(presentInCulture == TRUE) %>% distinct(Taxon, newTaxon)

#gets rid of extra spaces at the end of newTaxon
str_extract(seawater$newTaxon, " *$")
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace_all(newTaxon, " $", ""), NA))
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace_all(newTaxon, " $", ""), NA))
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace_all(newTaxon, " $", ""), NA))
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace_all(newTaxon, " $", ""), NA))
seawater <- seawater %>% mutate(newTaxon = ifelse(presentInCulture == TRUE, str_replace_all(newTaxon, " $", ""), NA))

str_extract(seawater$newTaxon, " *$")


#makes a variable for the shortened taxonomy
seawater <- seawater %>% mutate(shortTaxon = str_extract(newTaxon, "[A-Z][a-z\\s0-9I]*$"))
seawater <- seawater %>% mutate(shortTaxon = ifelse(str_detect(newTaxon, "GMD14H09"), str_extract(newTaxon, "GMD14H09"), shortTaxon))
seawater <- seawater %>% mutate(shortTaxon = ifelse(str_detect(newTaxon, "HTCC"), str_extract(newTaxon, "HTCC"), shortTaxon))

#seawater %>% filter(presentInCulture == TRUE) %>% distinct(Taxon, newTaxon, shortTaxon) %>% View()

seawater %>% filter(prop == 0)

head(seawater)

seawater %>% distinct(OTU_ID) %>% nrow()

#arranges seawater by prop so that I can put the sequences in order on the figure
order <- seawater %>% arrange(prop) %>% distinct(OTU_ID)

#puts the sequences in order by prop for figure
seawater$OTU_ID <- factor(seawater$OTU_ID, levels = rev(c(order$OTU_ID)))

seawater %>% distinct(OTU_ID) %>% nrow()

seawater %>% mutate(prop = prop*200000) %>% ggplot(aes(x = OTU_ID, y = prop, fill = presentInCulture)) + geom_bar(stat = 'identity') +
  scale_fill_manual(values=c("grey", "purple")) + labs(y = "Reads per 200,000 reads in seawater", x = "16S rRNA sequence",
                                                       fill = "Present in at least one culture") + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=5)) + 
  geom_point(data = seawater %>% mutate(prop = prop*200000) %>% filter(presentInCulture == TRUE), aes(x = OTU_ID, y = prop), color = "purple", show.legend = FALSE, size = .1) +
  geom_bar(data = seawater %>% mutate(prop = prop*200000) %>% filter(presentInCulture == TRUE), aes(x = OTU_ID, y = prop), color = "purple", show.legend = FALSE, stat = 'identity') + 
  scale_y_log10() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank()) + theme(axis.ticks.x=element_blank())

ggsave("seawaterSequenceRelativeAbundanceComparedToCultures_logScale_clustered.png", dpi = 600, height = 6, width = 8)


#makes variable for whether taxonomy of sequence is Syn and/or Pro so that I can color the Syn/Pro labels
#a different color
seawater <- seawater %>% mutate(synPro = ifelse(str_detect(shortTaxon, "Syn|Pro"), TRUE, FALSE))

seawater %>% mutate(prop = prop*200000) %>% ggplot(aes(x = OTU_ID, y = prop, fill = presentInCulture)) + geom_bar(stat = 'identity') +
  scale_fill_manual(values=c("grey", "purple")) + labs(y = "Reads per 200,000 reads in seawater", x = "16S rRNA sequence",
                                                       fill = "Present in at least one culture") + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=5)) + 
  geom_point(data = seawater %>% mutate(prop = prop*200000) %>% filter(presentInCulture == TRUE), aes(x = OTU_ID, y = prop), color = "purple", show.legend = FALSE, size = .1) + 
  geom_bar(data = seawater %>% mutate(prop = prop*200000) %>% filter(presentInCulture == TRUE), aes(x = OTU_ID, y = prop), color = "purple", show.legend = FALSE, stat = 'identity') + 
  geom_text(data = seawater %>% filter(presentInCulture == TRUE), aes(x = OTU_ID, y = (prop*200000)+20, label = shortTaxon, color = synPro), size = 3, angle = 90) +
  scale_y_log10() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

ggsave("seawaterSequenceRelativeAbundanceComparedToCultures_logScale_withLabels_clustered.svg", dpi = 600, height = 10, width = 12)


seawater %>% filter(presentInCulture == TRUE)
seawater %>% filter(presentInCulture == TRUE) %>% distinct(Taxon)
seawater %>% filter(presentInCulture == TRUE) %>% filter(str_detect(Taxon, "\\[\\]"))


#seawater %>% filter(presentInCulture == TRUE) %>% arrange(desc(prop)) %>% distinct(OTU_ID, shortTaxon) %>% View()


seawater %>% filter(presentInCulture == TRUE) %>% distinct(OTU_ID) %>% nrow()

#there are two Syn/Pro clusters that are in both seawater and culture samples
seawater %>% filter(presentInCulture == TRUE) %>% filter(synPro == TRUE)

#for the two Syn/Pro clusters that are in both seawater and culture samples, add the culture info
df <- seawater %>% filter(presentInCulture == TRUE) %>% filter(synPro == TRUE) %>% left_join(culture, by = c("OTU_ID"))

head(df)

#makes a variable for whether sample is a Pro or Syn culture
df <- df %>% mutate(culture = ifelse(str_detect(sample.y, "SYN|WH|9220|S9503|8102|RS9916"), "Synechococcus", "Prochlorococcus"))

#the features in 0602 and 0603 got identified as Synechococcaceae rather than Synechococcus
#like they have previously
df %>% filter(sample.y %in% c("0602", "0603")) %>% select(shortTaxon)

#I cannot identify 00591adaa51f7754dda8bee8349eb9b11c22f3b6 more than Synechococcaceae because 
#this cluster is found in both Pro and Syn cultures, so it is composed of both Pro and Syn cultures
df %>% distinct(OTU_ID, shortTaxon, culture)

#loads blast results of clustered features compared to Alt Mac MIT1002 (lab strain of Alt Mac)
#this came from analysisSeanAndSeawaterClustered.txt
blast <- read.table("MIT1002_vs_sequences_noContamination_noMock_withProSyn_no1223_plusSeawater_clusteredFeatures")

str(blast)

#none of the clustered features matched perfectly to Alt Mac Mit1002
blast %>% filter(V3 == 100)
blast %>% arrange(desc(V3)) %>% head()
seawater %>% filter(OTU_ID == "00bd38f89dba01002db00431e5c5ea66d06d6190") %>% View()


