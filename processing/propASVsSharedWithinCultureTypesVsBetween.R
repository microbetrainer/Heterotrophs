library(tidyverse)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads corrected counts of each taxa in each sample 
#this excludes sequences that were determined to be contamination from nearby wells
seq <- read.table("seqtab_corrected.txt")

numCols <- ncol(seq)

cols <- str_c("col", 1:numCols)

#gives a number to each sequence (column)
colnames(seq) <- cols

#gets list of samples
samples <- rownames(seq)

#makes a variable for the sample in seq instead of the samples just being
#the rownames
seq <- seq %>% mutate(sample = samples)

#gathers seq into long form
seqG <- seq %>% gather(1:3170, key = "sequence", value = "abundance")

#gets the total number of sequences in each sample
summ <- seqG %>% group_by(sample) %>% summarize(totalSeqs = sum(abundance))

#adds number of sequences in sample to each sequence
seqG <- seqG %>% left_join(summ, by = c("sample"))

#gets rid of the rows for the samples that are unrelated (Sean told me which ones to exclude)
seqG <- seqG %>% filter(!str_detect(sample, "Con|Epi|Pcn|Lin"))

#gets rid of Mock samples
seqG <- seqG %>% filter(!str_detect(sample, "Mock"))

#gets rid of the rows for which a sequence is not present in a sample
#this excludes the sequences that are just present in the unrelated samples (Con|Epi|Pcn|Lin)
seqG <- seqG %>% filter(abundance > 0)

#gets the proportion that each sequence is in each sample
seqG <- seqG %>% mutate(propSample = abundance/totalSeqs)

#excludes the sequences in samples that have a low relative abundances in the sample
seqG <- seqG %>% filter(propSample >= 0.002)

#gets rid of unnecessary variables
seqG <- seqG %>% select(-c(abundance, totalSeqs))

#makes sequence names in seqG match sequence names in repSeqs
seqG$sequence <- str_replace(seqG$sequence, "col", "contig_")

#loads the representative sequences fasta file so that I can get the 
#contig ID assigned to each feature ID
#this is with excluding Pro and Syn
repSeqs <- read.table("representative_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.fasta", fill = TRUE)

#gets just the rows in repSeqs that correspond to IDs
repSeqs <- repSeqs %>% filter(str_detect(V1, ">"))

#adds contig IDs to feature IDs in seqG
seqG <- seqG %>% left_join(repSeqs, by = c("sequence" = "V2"))

#gets rid of ">" in feature IDs
seqG <- seqG %>% mutate(V1 = str_replace(V1, ">", ""))

#makes feature ID, contig ID combined variable that matches 
#the IDs used in the tree file
seqG <- seqG %>% mutate(featureContig = str_c(V1, sequence, sep = ""))

#loads cleaned up taxonomy of the same length features
#corresponds to representative_sameLengthSeqs.fasta
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv")

#adds taxonomy to seqG
seqG <- seqG %>% left_join(tax, by = c("V1" = "Feature ID"))

#gets rid of pro and syn sequences
seqG <- seqG %>% filter(!is.na(shortTaxa))

#JW3 only has Pro sequences so it gets excluded which is good because it only has 
#one sequence
seqG <- seqG %>% filter(sample != "JW3")

#excludes sample 1223
seqG <- seqG %>% filter(sample != "1223")

#according to Sean:
#pro: JW7, B7, DV, C12B, JW2, JW4, 1213
#syn: 8102, RS9916

#gets samples that are either Pro or Syn, which is all of them
seqG <- seqG %>% filter(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916|SYN|WH|9220|S9503|8102|RS9916|0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223|JW7|B7|DV|C12B|JW2|JW4|1213"))

#makes a variable for whether sample is Pro or Syn
seqG <- seqG %>% mutate(culture = ifelse(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916"), "Synechococcus", "Prochlorococcus"))

##from comparingSampleIDs.R
#1213 is 1313
#WH7803 is WH7803
#SYN1320 is SYN1220
#9201.2 is 9311

#AS9601 is correct ID rather than ASN9601
#SYN9503 is correct ID rather than S9503

seqG <- seqG %>% mutate(sample = ifelse(sample == "1213", "1313", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "9201.2", "9311", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "ASN9601", "AS9601", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "9201.1", "9201", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "S9503", "SYN9503", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "SYN1320", "SYN1220", sample))

#seqG$sample <- ifelse(seqG$sample == "8102", "SYNCLADEXVI", seqG$sample)
seqG <- seqG %>% mutate(sample = ifelse(sample == "8102", "SYNCLADEXVI", sample))

noProSyn <- seqG

#doesn't include pro, syn sequences
head(noProSyn)
noProSyn %>% distinct(sequence) %>% nrow()

#sequences that are found in both pro and syn cultures
noProSyn %>% group_by(sequence) %>% distinct(culture) %>% summarize(n = n()) %>% filter(n == 2) %>% nrow()

#sequences that are found in just one type of cultures
noProSyn %>% group_by(sequence) %>% distinct(culture) %>% summarize(n = n()) %>% filter(n == 1) %>% nrow()

54+181

#proportion of sequences that are found in both pro and syn cultures
54/235

noProSyn %>% filter(str_detect(shortTaxa, "Alteromonas")) %>% distinct(sample) %>% nrow()
noProSyn %>% filter(str_detect(shortTaxa, "Marinobacter")) %>% distinct(sample) %>% nrow()

#number of samples of each culture type
noProSyn %>% group_by(culture) %>% distinct(sample) %>% summarize(n = n()) 

#gets the number of samples of each culture type that sequences are found in
noProSyn %>% group_by(sequence, culture) %>% distinct(sample) %>% summarize(n = n()) 

noProSyn %>% group_by(sequence, culture) %>% distinct(sample) %>% summarize(n = n()) %>% 
  group_by(culture) %>% summarize(max = max(n), min = min(n))

str(noProSyn)
noProSyn %>% summarize(min = min(propSample))

noProSyn %>% group_by(sample) %>% summarize(totalProp = sum(propSample))

str(noProSyn)


#makes a list of the samples
samples <- noProSyn %>% distinct(sample)
samples
samples <- samples %>% as.list()
samples
length(samples)
samples <- samples$sample
samples
length(samples)

#initializes lists
propSharedList = {}
samplesList = {}

set.seed(7)

for (i in 1:100){
  
  #choose two random samples with replacement
  random2Samples <- sample(samples, 2, replace = FALSE)
  #makes dataframe of these samples
  random2Samples <- random2Samples %>% as.data.frame()
  #fixes column name
  colnames(random2Samples) <- "sample"
  
  #gets the rows in noProSyn for the two selected samples
  noProSyn_subset <- noProSyn %>% semi_join(random2Samples, by = c("sample"))
  #gets the number of distinct sequences between the two samples
  totalSeq <- noProSyn_subset %>% distinct(sequence) %>% nrow()
  #gets the number of sequence that are present in both samples
  sharedSeq <- noProSyn_subset %>% group_by(sequence) %>% distinct(sample) %>% summarize(n = n()) %>% 
    filter(n == 2) %>% nrow()
  sharedSeq/totalSeq
  random2Samples
  
  #adds the culture types to the dataframe of samples
  random2Samples <- random2Samples %>% left_join(noProSyn %>% distinct(sample, culture), by = c("sample"))
  
  #if the samples are of the same cultures, make TOrF variable TRUE
  random2Samples %>% distinct(culture) %>% summarize(n = n())
  TOrF <- random2Samples %>% distinct(culture) %>% summarize(n = n()) == 1
  TOrF
  #makes sameCulture variable Pro if both samples are Pro cultures, Syn if both 
  #samples are Syn cultures, and Diff cultures if they are different cultures
  random2Samples <- random2Samples %>% mutate(sameCulture = ifelse(TOrF & culture[1] == "Prochlorococcus", "Pro", 
                                                                   ifelse(TOrF == TRUE, "Syn", "Diff cultures")))
  
  #adds proportion of sequence that are present in both samples to list
  propSharedList[[i]] <- sharedSeq/totalSeq
  #adds the pair of samples to list
  samplesList[[i]] <- random2Samples
  
}

propSharedList
samplesList

#propSharedList
#samplesList

#initializes list for extracting sameCulture values from dataframes
samplesListValue = {}

for (i in 1:100){
  samplesListValue[[i]] <- as.list(samplesList[[i]] %>% distinct(sameCulture))$sameCulture
}

#samplesListValue

#propSharedList
#samplesListValue

#makes dataframe with prop of sequences shared, and culture type of pair of samples
df <- data.frame(samplePair = samplesListValue, propSeqsShared = propSharedList)

#it doesn't seem that there is a difference in the prop of sequences shared between 
#samples within Syn samples, Pro samples, and between the two groups because the results 
#change so much depending on the number I put in set.seed and how many pairs I sample
str(df)
df %>% group_by(samplePair) %>% summarize(meanPropSeqsShared = mean(propSeqsShared)*100, numPairs = n())

78*77
56*55
18*17
