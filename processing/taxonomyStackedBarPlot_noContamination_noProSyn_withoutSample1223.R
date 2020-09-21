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
nrow(seq)

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

seqG %>% distinct(sequence) %>% nrow() 

write_csv(seqG, "seqG_forDiatomAndSeanCombinedAnalysis.csv")

#this includes Pro and Syn sequences
seqG %>% distinct(sequence) %>% nrow()


#makes sequence names in seqG match sequence names in repSeqs
seqG$sequence <- str_replace(seqG$sequence, "col", "contig_")


seqG %>% distinct(sequence) %>% nrow()


#loads cleaned up taxonomy of the features
#this doesn't include Pro and Syn
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv")

tax %>% nrow()

#loads the representative sequence IDs for the features
#this doesn't include Pro and Syn and corresponds to 
#taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv
repSeqs <- read.table("representative_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.fasta", fill = TRUE)

#gets just the rows that correspond to IDs
repSeqs <- repSeqs %>% filter(str_detect(V1, ">"))

#gets rid of ">" in ID names
repSeqs$V1 <- str_replace(repSeqs$V1, ">", "")

#same number of sequences as tax which is good
repSeqs %>% nrow()


nrow(tax)

#adds "contig_..." format ID names to tax so tax can be added to seqG
tax <- tax %>% left_join(repSeqs, by = c("Feature ID" = "V1"))

nrow(tax)

tax %>% filter(is.na(V2))

nrow(seqG)

#adds taxonomy to seqG
seqG <- seqG %>% left_join(tax, by = c("sequence" = "V2"))

nrow(seqG)

#the rows with NA taxonomies are the Pro and Syn sequences
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

#gets rid of unnecessary variables from the tax dataframe
seqG <- seqG %>% select(c(1:5, Taxon, taxa))

#the non-Pro, non-Syn sequences are the rows that have taxonomies
seqG %>% filter(!is.na(taxa)) %>% distinct(sequence) %>% nrow()

nrow(tax)

#this is from analysisSameLengthSequences.txt
#loads taxonomy dataframe for which I did not exclude Pro and Syn sequences
taxFull <- read_tsv("taxonomySameLengthSequences.tsv")

#gets rid of first row that just has notes
taxFull <- taxFull[-1,]

#there are 649 features because this is from when I ran qiime on the uncorrected seqtab 
#so now, there are fewer features because now I have excluded the low relative abundance 
#sequences in samples and the sequences I determined to be contamination
nrow(taxFull)

nrow(seqG)

#this is from analysisSameLengthSequences.txt
#loads the representative sequences fasta file so that I can get the 
#contig ID assigned to each feature ID
#this is based on running qiime on uncorrected seqtab
repSeqsFull <- read.table("representative_sameLengthSeqs.fasta", fill = TRUE)

#gets just the rows in repSeqs that correspond to IDs
repSeqsFull <- repSeqsFull %>% filter(str_detect(V1, ">"))

#gets rid of ">" in feature IDs
repSeqsFull$V1 <- str_replace(repSeqsFull$V1, ">", "")

#makes V2 match the sequence IDs in seqG
repSeqsFull$V2 <- str_replace(repSeqsFull$V2, "contig_", "col")

#there are 649 features because this is from when I ran qiime on the uncorrected seqtab 
#so now, there are fewer features because now I have excluded the low relative abundance 
#sequences in samples and the sequences I determined to be contamination
nrow(repSeqsFull)


#adds "col.." IDs to tax so that I can join tax to seqG
nrow(taxFull)

taxFull <- taxFull %>% left_join(repSeqsFull, by = c("Feature ID" = "V1"))

nrow(taxFull)

taxFull %>% filter(is.na(V2))

#makes V2 values match seqG sequence values
taxFull$V2 <- str_replace(taxFull$V2, "col", "contig_")

#adds taxonomies to seqG
nrow(seqG)

seqG <- seqG %>% left_join(taxFull, by = c("sequence" = "V2"))

nrow(seqG)

seqG %>% filter(is.na(Taxon.y))

#still 268 sequences which is good
seqG %>% distinct(sequence) %>% nrow()

#gets rid of unnecessary variables that came from taxFull
seqG <- seqG %>% select(c(1:7,9))


seqG %>% filter(Taxon.x != Taxon.y) %>% select(Taxon.x, Taxon.y)

#the only rows with a NA Taxon.x value are Pro and Syn
seqG %>% filter(is.na(Taxon.x)) %>% distinct(Taxon.y)

#makes a new taxonomy variable
#if Taxon.x is NA because it is a Pro or Syn sequence, makes the new taxonomy variable Taxon.y
seqG$newTaxon <- ifelse(is.na(seqG$Taxon.x), seqG$Taxon.y, seqG$Taxon.x)

seqG %>% filter(is.na(newTaxon))

seqG %>% distinct(sample) %>% nrow()
seqG %>% distinct(sequence) %>% nrow()

seqG %>% filter(abundance == 0)
seqG %>% filter(propSample < 0.002)

#makes a variable for the class
seqG <- seqG %>% mutate(class = str_extract(newTaxon, "c__[a-z,A-Z,0-9,\\[,\\]]*;{0,1}"))
seqG %>% distinct(class)
seqG <- seqG %>% mutate(class = str_replace_all(class, " ", ""))
seqG <- seqG %>% mutate(class = str_replace_all(class, ";", ""))
seqG %>% distinct(class)
seqG %>% filter(is.na(class))
seqG$class <- ifelse(is.na(seqG$class), "Unclassified at class level", seqG$class)
seqG %>% filter(class == "c__")
seqG %>% distinct(class)

#gets the total number of reads in each sample
#this is after excluding the low propSample reads from each sample that I determined to be contamination
summ <- seqG %>% group_by(sample) %>% summarize(totalSeqsWithoutLowProp = sum(abundance))

nrow(summ)

nrow(seqG)

#adds total number of reads in each sample (after excluding the low propSample reads from each sample) 
#to seqG
seqG <- seqG %>% left_join(summ, by = c("sample"))

nrow(seqG)

seqG %>% filter(is.na(totalSeqs))

#makes a variable for the relative abundance for each sequence in each sample
seqG <- seqG %>% mutate(propSampleWithoutLowProp = abundance/totalSeqsWithoutLowProp)

#the difference in the total seqs in each sample, with and without excluding low prop sequences is 
#not large relative to the total number of seqs
seqG %>% mutate(totalDiff = abs(totalSeqsWithoutLowProp - totalSeqs)) %>% distinct(sample, totalSeqsWithoutLowProp, totalDiff) %>% arrange(desc(totalDiff)) %>% slice(1:6)

#the difference in the relative abundance of each sequence in a sample, with and without excluding low prop 
#sequences iin the sample is not large relative to relative abundance of the sequence in the sample
seqG %>% mutate(propDiff = abs(propSampleWithoutLowProp - propSample)) %>% select(sample, sequence, propSampleWithoutLowProp, propDiff) %>% arrange(desc(propDiff)) %>% slice(1:6)


#there are 12 classes
seqG %>% distinct(class) %>% nrow()

#gets rid of "c__" in class names
seqG$class <- str_replace(seqG$class, "c__", "")

seqG %>% distinct(class) %>% nrow()
seqG %>% distinct(class)

#gets rid of "[" and "]" in class names
#seqG$class <- str_replace(seqG$class, "\\[", "")
#seqG$class <- str_replace(seqG$class, "\\]", "")

seqG %>% distinct(class) %>% nrow()
seqG %>% distinct(class)

seqG %>% distinct(sample) %>% nrow()

#Pro, Synechococcaceae, and Synechococcus have the same class
seqG %>% filter(is.na(taxa)) %>% distinct(newTaxon, class)

#all of the samples are either Pro or Syn cultures
seqG %>% filter(!str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916|SYN|WH|9220|S9503|8102|RS9916|0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223|JW7|B7|DV|C12B|JW2|JW4|1213")) %>% 
  distinct(sample) %>% nrow()

#gets samples that are either Pro or Syn, which is all of them
seqG <- seqG %>% filter(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916|SYN|WH|9220|S9503|8102|RS9916|0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223|JW7|B7|DV|C12B|JW2|JW4|1213"))

seqG %>% distinct(sample) %>% nrow()

#makes a variable for whether sample is a Pro or Syn culture
seqG <- seqG %>% mutate(culture = ifelse(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916"), "Synechococcus", "Prochlorococcus"))

#theses samples have more than one Pro/Syn sequence
#I look at these sequences in getSequencesForTree_noContamination_noProSyn_withoutSample1223.R
seqG %>% group_by(sample) %>% filter(class == "Synechococcophycideae") %>% summarize(n = n()) %>% filter(n != 1)
seqG %>% group_by(sample) %>% filter(class == "Synechococcophycideae") %>% filter(sample %in% c("0602", "1300", "ASN9601", "C12B")) %>% 
  select(sample, sequence, abundance, newTaxon) %>% arrange(sample)

#if the class is Synechococcophycideae, adds whether the sequence is Pro or Syn
seqG <- seqG %>% mutate(class = ifelse(class == "Synechococcophycideae", str_c(class, culture, sep = " "), class))

seqG %>% distinct(class)

##from comparingSampleIDs.R
#1213 is 1313
#WH7803 is WH7801
#SYN1320 is SYN1220
#9201.2 is 9311

#AS9601 is correct ID rather than ASN9601
#SYN9503 is correct ID rather than S9503


seqG$sample <- ifelse(seqG$sample == "1213", "1313", seqG$sample)
#seqG$sample <- ifelse(seqG$sample == "WH7803", "WH7801", seqG$sample)
seqG$sample <- ifelse(seqG$sample == "SYN1320", "SYN1220", seqG$sample)
seqG$sample <- ifelse(seqG$sample == "9201.2", "9311", seqG$sample)
seqG$sample <- ifelse(seqG$sample == "9201.1", "9201", seqG$sample)

seqG$sample <- ifelse(seqG$sample == "ASN9601", "AS9601", seqG$sample)
seqG$sample <- ifelse(seqG$sample == "S9503", "SYN9503", seqG$sample)

seqG$sample <- ifelse(seqG$sample == "8102", "SYNCLADEXVI", seqG$sample)

seqG %>% distinct(sample) %>% nrow()


#from heatmap_noContamination_noLowAbundance_noMock_noProSyn_without1223_propReads.R
sampleOrder <- read_csv("orderOfClusteredSamples.csv")

list <- sampleOrder$V1 %>% as.list()
list

seqG$sample <- factor(seqG$sample, levels = list)

levels(seqG$sample)
seqG %>% distinct(class)

#write_csv(seqG, "taxonomyStackedBarPlot_noContamination_noMock_withProSyn_withoutSample1223_data.csv")

#http://sape.inf.usi.ch/quick-reference/ggplot2/colour
colors <- c("cyan3", "mediumvioletred", "gray58", "springgreen4", "hotpink1", "red", "tan3", "blueviolet", "black", "blue1", "springgreen", "pink", "darkorange1")


length(colors)

seqG %>% distinct(class) %>% arrange(class) %>% nrow()

seqG %>% distinct(class) %>% arrange(class) %>% as.list()

order <-  c("A712011", "Alphaproteobacteria", "Betaproteobacteria", "Cytophagia", "Deltaproteobacteria",                  
            "Flavobacteriia", "Gammaproteobacteria", "[Leptospirae]", "Phycisphaerae", "[Rhodothermi]",                          
            "Synechococcophycideae Prochlorococcus", "Synechococcophycideae Synechococcus", "Unclassified at class level")          

length(order)

seqG$class <- factor(seqG$class, levels = order)

levels(seqG$class)

seqG %>% ggplot(aes(x = sample, y = propSampleWithoutLowProp, fill = class)) + geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = colors) + labs(x = "", y = "Relative abundance", fill = "")

ggsave("taxonomyStackedBarPlot_noContamination_noMock_withProSyn_withoutSample1223.png", dpi = 600, height = 6, width = 12)
ggsave("taxonomyStackedBarPlot_noContamination_noMock_withProSyn_withoutSample1223.svg", dpi = 600, height = 6, width = 12)


#makes dataframe without Pro or Syn
noProSyn <- seqG %>% filter(class != "Synechococcophycideae Prochlorococcus") %>% filter(class != "Synechococcophycideae Synechococcus")

noProSyn %>% distinct(class)

#gets rid of unnecessary variables
noProSyn <- noProSyn %>% select(-c(totalSeqs, propSample, totalSeqsWithoutLowProp, propSampleWithoutLowProp))

#gets the total number of reads in each sample
#this is after excluding the Pro and Syn sequences
summ <- noProSyn %>% group_by(sample) %>% summarize(totalSeqsWithoutProSyn = sum(abundance))

nrow(summ)

nrow(noProSyn)

#adds total number of reads in each sample after excluding 
#Pro and Syn sequences to noProSyn
noProSyn <- noProSyn %>% left_join(summ, by = c("sample"))

nrow(noProSyn)

noProSyn %>% filter(is.na(totalSeqsWithoutProSyn))

#calculates the relative abundance of each sequence in each sample
#after excluding Pro and Syn sequences
noProSyn <- noProSyn %>% mutate(propSampleWithoutProSyn = abundance/totalSeqsWithoutProSyn)


#http://sape.inf.usi.ch/quick-reference/ggplot2/colour
colors <- c("cyan3", "mediumvioletred", "gray58", "springgreen4", "hotpink1", "red", "tan3", "blueviolet", "black", "blue1", "darkorange1")


length(colors)

noProSyn %>% distinct(class) %>% arrange(class) %>% nrow()
noProSyn %>% distinct(class) %>% arrange(class)

order <-  c("A712011", "Alphaproteobacteria", "Betaproteobacteria", "Cytophagia", "Deltaproteobacteria",                  
            "Flavobacteriia", "Gammaproteobacteria", "[Leptospirae]", "Phycisphaerae", "[Rhodothermi]",                          
            "Unclassified at class level")          

length(order)

noProSyn$class <- factor(noProSyn$class, levels = order)

levels(noProSyn$class)

noProSyn %>% ggplot(aes(x = sample, y = propSampleWithoutProSyn, fill = class)) + geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = colors) + labs(x = "", y = "Relative abundance", fill = "")

ggsave("taxonomyStackedBarPlot_noContamination_noMock_noProSyn_withoutSample1223.png", dpi = 600, height = 6, width = 12)
ggsave("taxonomyStackedBarPlot_noContamination_noMock_noProSyn_withoutSample1223.svg", dpi = 600, height = 6, width = 12)


