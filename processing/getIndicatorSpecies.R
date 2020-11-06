library(tidyverse)
library(readr)
library(readxl)
library(picante)
library(phyloseq)
library(randomForest)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#https://cran.r-project.org/web/packages/indicspecies/indicspecies.pdf
#http://ichthyology.usm.edu/courses/multivariate/apr_11.pdf

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

#gathers seq into long form
seqG <- seq %>% gather(1:3170, key = "sequence", value = "abundance")



#gets the total number of sequences in each sample
summ <- seqG %>% group_by(sample) %>% summarize(totalSeqs = sum(abundance))

nrow(seqG)

#adds number of sequences in sample to each sequence
seqG <- seqG %>% left_join(summ, by = c("sample"))

nrow(seqG)

seqG %>% filter(is.na(totalSeqs))

#gets rid of the rows for the samples that are unrelated (Sean told me which ones to exclude)
seqG <- seqG %>% filter(!str_detect(sample, "Con|Epi|Pcn|Lin"))

seqG %>% filter(str_detect(sample, "Mock")) %>% distinct(sample)

#gets rid of Mock samples
seqG <- seqG %>% filter(!str_detect(sample, "Mock"))

#gets rid of the rows for which a sequence is not present in a sample
#this excludes the sequences that are just present in the unrelated samples (Con|Epi|Pcn|Lin)
seqG <- seqG %>% filter(abundance > 0)

seqG %>% distinct(sequence) %>% nrow()

#gets the proportion that each sequence is in each sample
seqG <- seqG %>% mutate(propSample = abundance/totalSeqs)



seqG %>% distinct(sequence) %>% nrow() 

#excludes the sequences in samples that have a low relative abundances in the sample
seqG <- seqG %>% filter(propSample >= 0.002)

seqG %>% distinct(sequence) %>% nrow()

#gets rid of unnecessary variables
seqG <- seqG %>% select(-c(propSample, totalSeqs))

head(seqG)

#makes sequence names in seqG match sequence names in repSeqs
seqG$sequence <- str_replace(seqG$sequence, "col", "contig_")

seqG %>% distinct(sequence) %>% nrow()

seqG %>% distinct(sample) %>% nrow()


#loads the representative sequences fasta file so that I can get the 
#contig ID assigned to each feature ID
repSeqs <- read.table("representative_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.fasta", fill = TRUE)

#gets just the rows in repSeqs that correspond to IDs
repSeqs <- repSeqs %>% filter(str_detect(V1, ">"))

head(repSeqs)

repSeqs %>% distinct(V2) %>% nrow()

seqG %>% distinct(sequence) %>% nrow()

nrow(repSeqs)

head(seqG)
head(repSeqs)

nrow(seqG)

#adds contig IDs to feature IDs in seqG
seqG <- seqG %>% left_join(repSeqs, by = c("sequence" = "V2"))

nrow(seqG)

head(seqG)

#these are the Pro and Syn sequences
seqG %>% filter(is.na(V1))

#excludes the Pro and Syn sequences 
seqG <- seqG %>% filter(!is.na(V1))

seqG %>% distinct(sequence) %>% nrow()

head(seqG)

#gets rid of ">" in feature IDs
seqG <- seqG %>% mutate(V1 = str_replace(V1, ">", ""))

head(seqG)

#makes feature ID, contig ID combined variable that matches 
#the IDs used in the tree file
seqG <- seqG %>% mutate(featureContig = str_c(V1, sequence, sep = ""))

head(seqG)

seqG %>% distinct(featureContig) %>% nrow()


#loads cleaned up taxonomy of the same length features
#corresponds to representative_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.fasta
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv") %>% select(`Feature ID`, shortTaxa)

tax %>% distinct(`Feature ID`) %>% nrow()

head(seqG)
head(tax)

nrow(seqG)

#adds taxonomy to seqG
seqG <- seqG %>% left_join(tax, by = c("V1" = "Feature ID"))

nrow(seqG)

seqG %>% filter(is.na(shortTaxa))


#JW3 only has Pro sequences so it gets excluded which is good because it only has 
#one sequence
seqG <- seqG %>% filter(sample != "JW3")

#excludes sample 1223
seqG <- seqG %>% filter(sample != "1223")

seqG %>% distinct(sample) %>% nrow()

#seqG %>% filter(!str_detect(sample, "SYN|WH|9220|S9503")) %>% 
#filter(!str_detect(sample, "0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223")) %>% 
#filter(str_detect(shortTaxa, "Synecho")) %>% distinct(sample)

#seqG %>% filter(!str_detect(sample, "SYN|WH|9220|S9503")) %>% 
#filter(!str_detect(sample, "0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223")) %>% 
#filter(str_detect(shortTaxa, "Prochlo")) %>% distinct(sample)



#seqG %>% filter(!str_detect(sample, "SYN|WH|9220|S9503")) %>% 
#filter(!str_detect(sample, "0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223")) %>% 
#distinct(sample)


#according to Sean:
#pro: JW7, B7, DV, C12B, JW2, JW4, 1213
#syn: 8102, RS9916

#there are 74 samples
seqG %>% distinct(sample) %>% nrow()

seqG %>% distinct(sequence) %>% nrow()

#all of the samples are either Pro or Syn
seqG %>% filter(!str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916|SYN|WH|9220|S9503|8102|RS9916|0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223|JW7|B7|DV|C12B|JW2|JW4|1213")) %>% 
  distinct(sample) %>% nrow()


#gets samples that are either Pro or Syn, which is all of them
seqG <- seqG %>% filter(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916|SYN|WH|9220|S9503|8102|RS9916|0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223|JW7|B7|DV|C12B|JW2|JW4|1213"))

#makes a variable for whether sample is Pro or Syn
seqG <- seqG %>% mutate(culture = ifelse(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916"), "Syn", "Pro"))

seqG %>% distinct(sample) %>% nrow()
seqG %>% distinct(sequence) %>% nrow()

#correct
seqG %>% filter(culture == "Syn") %>% distinct(sample)

#correct
seqG %>% filter(culture == "Pro") %>% distinct(sample) %>% nrow()

#excludes Pro sequences
seqG <- seqG %>% filter(!str_detect(shortTaxa, "Prochlo"))

#exludes Synechococcaceae and Synechococcus sequences
seqG <- seqG %>% filter(!str_detect(shortTaxa, "Syn"))

head(seqG)

seqG %>% group_by(featureContig) %>% distinct(V1) %>% summarize(n = n()) %>% filter(n > 1)
seqG %>% group_by(featureContig) %>% distinct(sequence) %>% summarize(n = n()) %>% filter(n > 1)

head(seqG)

#gets rid of variables so I can spread data
seqG <- seqG %>% select(-c(V1, sequence))
seqG <- seqG %>% select(-shortTaxa)

head(seqG)

#doesn't include Pro or Syn
seqG %>% distinct(sample) %>% nrow()
seqG %>% distinct(featureContig) %>% nrow()

#spreads seqG so that the sequences that are not in a sample will have 0 abundance for that sample
seqSpread <- seqG %>% spread(key = featureContig, value = abundance, fill = 0)

dim(seqSpread)

cultureCol <- seqSpread$culture

seqSpread <- seqSpread %>% select(-c(sample, culture))

#makes column names not start with numbers
colnames(seqSpread) <- str_c("col", colnames(seqSpread))

dim(seqSpread)

library(indicspecies)

set.seed(7)
wetpt <- multipatt(seqSpread, cultureCol, control = how(nperm=999))

summary(wetpt)

#makes seqSpread presence, absence dataframe instead of abundance
seqSpread_pa <- seqSpread
seqSpread_pa[seqSpread_pa > 0] <- 1


#I copied results from here to indicatorSpeciesResults.txt

#I set func to "r.g" because https://cran.r-project.org/web/packages/indicspecies/vignettes/indicspeciesTutorial.pdf
#says "It is a good practice to correct the phi coefficient for the fact that some
#groups have more sites than others [Tich´y and Chytr´y, 2006]. To do that,
#we need to use func = "r.g" instead of func = "r""
set.seed(7)
summary(multipatt(seqSpread_pa, cultureCol, control = how(nperm=999), func = "r.g"))



res <- read.table("indicatorSpeciesResults.txt", fill = TRUE)

#makes variable for whether ASV is indicator in Pro or Syn culture
res %>% filter(V1 == "col300f0019e52c003e5defa26315520fd164ae5f8ccontig_17")
res <- res %>% mutate(indicCulture = ifelse(V1 == "col300f0019e52c003e5defa26315520fd164ae5f8ccontig_17", "Pro", "Syn"))

#gets rid of unnecessary rows 
res %>% filter(str_detect(V1, "col")) %>% nrow()
res <- res %>% filter(str_detect(V1, "col"))

head(res)

#gets rid of column with asterisks
res <- res %>% select(-V4)

head(res)

#renames variables
colnames(res) <- c("ASV", "Association value", "Unadjusted P-Value", "Culture type with this ASV as an indicator")

nrow(res)

head(res)

#gets rid of "col" at the beginning of sequence IDs
res$ASV <- str_replace(res$ASV, "col", "")

#gets rid of contig... at the end of sequence IDs
res$ASV <- str_replace(res$ASV, "contig.*", "")

head(res)


#loads cleaned up taxonomy of the same length features
#corresponds to representative_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.fasta
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv") %>% select(`Feature ID`, shortTaxa)

head(res)
head(tax)

nrow(res)

#adds taxonomy to indicspecies results
res <- res %>% left_join(tax, by = c("ASV" = "Feature ID"))

nrow(res)

res %>% filter(is.na(shortTaxa)) %>% nrow()

res

#contig_7 matches the lab strain Alt Mac 1002

repSeqs %>% filter(V2 == "contig_7")

#Alt Mac MIT1002 is not in res
res %>% filter(ASV == "f1e1d1742851920183773cc5710788dbe75225e8")

res %>% filter(str_detect(shortTaxa, "Alte"))

write_csv(res, "indicatorSpeciesResults.csv")




