library(tidyverse)
library(readr)
library(readxl)
library(picante)
library(phyloseq)
library(randomForest)

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

head(seqG)

#gets rid of unnecessary variables
seqG <- seqG %>% select(-c(abundance, totalSeqs))

head(seqG)

#makes sequence names in seqG match sequence names in repSeqs
seqG$sequence <- str_replace(seqG$sequence, "col", "contig_")

head(seqG)

#loads the representative sequences fasta file so that I can get the 
#contig ID assigned to each feature ID
repSeqs <- read.table("representative_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.fasta", fill = TRUE)

#gets just the rows in repSeqs that correspond to IDs
repSeqs <- repSeqs %>% filter(str_detect(V1, ">"))

head(seqG)
head(repSeqs)

repSeqs %>% distinct(V2) %>% nrow()

seqG %>% distinct(sequence) %>% nrow()


nrow(repSeqs)

nrow(seqG)

#adds contig IDs to feature IDs in seqG
seqG <- seqG %>% left_join(repSeqs, by = c("sequence" = "V2"))

nrow(seqG)

head(seqG)

#these are the Pro and Syn sequences
seqG %>% filter(is.na(V1))

#excludes the Pro and Syn sequences
seqG <- seqG %>% filter(!is.na(V1))

head(seqG)

seqG %>% distinct(V1) %>% nrow()
seqG %>% distinct(sequence) %>% nrow()

#gets rid of ">" in feature IDs
seqG <- seqG %>% mutate(V1 = str_replace(V1, ">", ""))

head(seqG)

#makes feature ID, contig ID combined variable that matches 
#the IDs used in the tree file
seqG <- seqG %>% mutate(featureContig = str_c(V1, sequence, sep = ""))

seqG %>% distinct(featureContig) %>% nrow()

head(seqG)



#loads cleaned up taxonomy of the same length features
#corresponds to representative_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.fasta
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv") %>% select(`Feature ID`, shortTaxa)

#from analyzeAltBlast.R
#contig_7 matches the lab strain Alt Mac 1002
#this is the f1e1d1742851920183773cc5710788dbe75225e8 feature
repSeqs %>% filter(V2 == "contig_7")

nrow(tax)

head(tax)

tax <- tax %>% mutate(shortTaxa = ifelse(`Feature ID` == "f1e1d1742851920183773cc5710788dbe75225e8", "Alteromonas macleodii MIT1002", shortTaxa))

nrow(tax)

tax %>% filter(`Feature ID` == "f1e1d1742851920183773cc5710788dbe75225e8")

tax %>% filter(str_detect(shortTaxa, "Alter"))

tax %>% distinct(`Feature ID`) %>% nrow()

head(seqG)

nrow(seqG)

#adds taxonomy to seqG
seqG <- seqG %>% left_join(tax, by = c("V1" = "Feature ID"))

nrow(seqG)

seqG %>% filter(is.na(shortTaxa))

head(seqG)

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

#all of the samples are either Pro or Syn
seqG %>% filter(!str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916|SYN|WH|9220|S9503|8102|RS9916|0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223|JW7|B7|DV|C12B|JW2|JW4|1213")) %>% 
  distinct(sample) %>% nrow()


#gets samples that are either Pro or Syn, which is all of them
seqG <- seqG %>% filter(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916|SYN|WH|9220|S9503|8102|RS9916|0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223|JW7|B7|DV|C12B|JW2|JW4|1213"))

#makes a variable for whether sample is Pro or Syn
seqG <- seqG %>% mutate(culture = ifelse(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916"), "Syn", "Pro"))

seqG %>% distinct(sample) %>% nrow()

#correct
seqG %>% filter(culture == "Syn") %>% distinct(sample) %>% nrow()

#correct
seqG %>% filter(culture == "Pro") %>% distinct(sample) %>% nrow()


#excludes Pro sequences
seqG %>% filter(str_detect(shortTaxa, "Prochlo"))
seqG <- seqG %>% filter(!str_detect(shortTaxa, "Prochlo"))

#exludes Synechococcaceae and Synechococcus sequences
seqG %>% filter(str_detect(shortTaxa, "Syn"))
seqG <- seqG %>% filter(!str_detect(shortTaxa, "Syn"))


seqG %>% filter(propSample == 0)
seqG %>% filter(propSample < .002)

head(seqG)

seqG %>% group_by(featureContig) %>% distinct(V1) %>% summarize(n = n()) %>% filter(n > 1)
seqG %>% group_by(featureContig) %>% distinct(sequence) %>% summarize(n = n()) %>% filter(n > 1)

head(seqG)

write_csv(seqG, "seqG_fromRandomForestProSyn.csv")

#gets rid of variables so I can spread data
seqG <- seqG %>% select(-c(V1, sequence))
seqG <- seqG %>% select(-shortTaxa)

head(seqG)

#spreads seqG so that the sequences that are not in a sample will have 0 abundance for that sample
seqSpread <- seqG %>% spread(key = featureContig, value = propSample, fill = 0)

dim(seqSpread)

#gathers seqG again, now with rows for sequences that are not present in a sample
seqG <- seqSpread %>% gather(3:237, key = "V1", value = "propSample")

head(seqG)

seqG %>% distinct(sample) %>% nrow()

seqG %>% distinct(V1) %>% nrow()

seqG %>% group_by(V1) %>% distinct(sample) %>% summarize(n = n()) %>% filter(n != 74)

seqG %>% group_by(sample) %>% distinct(culture) %>% summarize(n = n()) %>% filter(n != 1)

set.seed(7)

#gets rid of sample variable so that I can run random forest
seqSpread <- seqSpread %>% select(-sample)

#makes column names not start with numbers
colnames(seqSpread) <- str_c("col", colnames(seqSpread))

#makes culture variable a factor so that I can run random forest
seqSpread$colculture <- factor(seqSpread$colculture)


seqSpread %>% filter(colculture == "Pro") %>% nrow()
seqSpread %>% filter(colculture == "Syn") %>% nrow()

head(seqSpread)

set.seed(7)

#I am sampling 18 at a time (with replacement I believe) because I want the same 
#number of samples for Pro and Syn and there are 18 Syn samples
#runs random forest
#changing mtry makes error rate worse
#https://stats.stackexchange.com/questions/90347/unbalanced-samples-random-forests
#https://stats.stackexchange.com/questions/24330/is-there-a-formula-or-rule-for-determining-the-correct-sampsize-for-a-randomfore
rf <- randomForest(colculture ~ ., data=seqSpread, sampsize = c(18,18))

rf

#        OOB estimate of  error rate: 8.11%
#Confusion matrix:
#  Pro Syn class.error
#Pro  53   3  0.05357143
#Syn   3  15  0.16666667



importance(rf) %>% head()

#makes dataframe of the importance of the variables
impDF <- importance(rf) %>% as.data.frame()

head(impDF)
nrow(impDF)

#makes variable for the sequence name
impDF <- impDF %>% mutate(V1 = rownames(impDF))

head(impDF)

#gets rid of "col" at the beginning of sequence name
impDF <- impDF %>% mutate(V1 = str_replace(V1, "^col", ""))

head(impDF)

#makes a variable for the sequence name without "contig..." at the end
impDF <- impDF %>% mutate(Feature.ID = str_replace(V1, "contig.*", ""))

head(impDF)

impDF %>% nrow()
impDF %>% distinct(Feature.ID) %>% nrow()

#arranges impDF from most important to least important
#higher Mean Decrease in Gini indicates higher variable importance
impDF <- impDF %>% arrange(desc(MeanDecreaseGini))

head(tax)

nrow(impDF)

#adds taxonomy to impDf
impDF <- impDF %>% left_join(tax, by = c("Feature.ID" = "Feature ID"))

nrow(impDF)

impDF %>% filter(is.na(shortTaxa))

head(seqG)
head(impDF)

#arranges impDF from most important to least important
#higher Mean Decrease in Gini indicates higher variable importance
impDF <- impDF %>% arrange(desc(MeanDecreaseGini))

head(impDF)

#gets the most important sequences from seqG
toPlot <- seqG %>% inner_join(impDF[1:20,], by = c("V1"))

head(toPlot)

toPlot %>% filter(is.na(shortTaxa))

toPlot %>% distinct(V1) %>% nrow()

toPlot %>% distinct(V1) %>% arrange(V1)

impDF[1:20,] %>% distinct(V1) %>% arrange(V1)

#makes combined variable for sequence ID and taxonomy
toPlot <- toPlot %>% mutate(featureTaxa = str_c(Feature.ID, shortTaxa, sep = ":\n"))

impDF <- impDF %>% mutate(featureTaxa = str_c(Feature.ID, shortTaxa, sep = ":\n"))

head(toPlot)
head(impDF)

impDF[1:20,] %>% distinct(featureTaxa)

order <- impDF[1:20,] %>% distinct(featureTaxa)

order <- order$featureTaxa %>% as.list()

toPlot$featureTaxa <- factor(toPlot$featureTaxa, levels = order)

levels(toPlot$featureTaxa)

head(toPlot)

toPlot %>% group_by(featureTaxa) %>% distinct(MeanDecreaseGini) %>% summarize(n = n()) %>% filter(n != 1)

summToPlot <- toPlot %>% group_by(culture, featureTaxa) %>% summarize(m = mean(propSample), sd = sd(propSample), meanDecGini = min(MeanDecreaseGini)) 

summToPlot %>% ggplot() +
  geom_bar(aes(x = culture, y = m), stat = 'identity') + labs(y = "Mean relative abundance", x = "") +
  geom_errorbar(aes(x = culture, ymin=m-sd, ymax=m+sd), width=.2, position=position_dodge(0.05)) + 
  geom_text(data = summToPlot %>% filter(culture == "Syn"), aes(label = round(meanDecGini, 2), x = culture, y = m + .02), position = position_jitter()) + 
  theme(text = element_text(size=6)) +
  facet_wrap(~featureTaxa, scales = 'free')

ggsave("randomForest20MostImportantVariablesToDetermineProOrSynCulture.png", dpi = 600, height = 10, width = 14)
ggsave("randomForest20MostImportantVariablesToDetermineProOrSynCulture.svg", dpi = 600, height = 10, width = 14)

###I am assuming that all Synechococcaceae are Synechococcus 






