library(tidyverse)
library(readr)
library(readxl)
library(picante)
library(phyloseq)
library(randomForest)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")


###I didn't fix this to have the correct repSeqs and taxonomy files

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
seqG <- seqG %>% select(-c(abundance, totalSeqs))


#makes sequence names in seqG match sequence names in repSeqs
seqG$sequence <- str_replace(seqG$sequence, "col", "contig_")

#loads the representative sequences fasta file so that I can get the 
#contig ID assigned to each feature ID
#this is without excluding Pro and Syn
repSeqs <- read.table("representative_sameLengthSeqs.fasta", fill = TRUE)

#gets just the rows in repSeqs that correspond to IDs
repSeqs <- repSeqs %>% filter(str_detect(V1, ">"))

repSeqs %>% distinct(V2) %>% nrow()

seqG %>% distinct(sequence) %>% nrow()


nrow(repSeqs)

nrow(seqG)

#adds contig IDs to feature IDs in seqG
seqG <- seqG %>% left_join(repSeqs, by = c("sequence" = "V2"))

nrow(seqG)

seqG %>% filter(is.na(V1))

#gets rid of ">" in feature IDs
seqG <- seqG %>% mutate(V1 = str_replace(V1, ">", ""))




#makes feature ID, contig ID combined variable that matches 
#the IDs used in the tree file
seqG <- seqG %>% mutate(featureContig = str_c(V1, sequence, sep = ""))



seqG %>% distinct(featureContig) %>% nrow()



#loads cleaned up taxonomy of the same length features
#corresponds to representative_sameLengthSeqs.fasta
tax <- read.table("taxonomySameLengthSequences_cleanedUp.tsv")

tax %>% distinct(Feature.ID) %>% nrow()


nrow(seqG)

#adds taxonomy to seqG
seqG <- seqG %>% left_join(tax, by = c("V1" = "Feature.ID"))

nrow(seqG)

seqG %>% filter(is.na(shortTaxa))



#JW3 only has Pro sequences so it gets excluded which is good because it only has 
#one sequence
seqG <- seqG %>% filter(sample != "JW3")

seqG %>% distinct(sample)

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

#there are 75 samples
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
seqG %>% filter(culture == "Syn") %>% distinct(sample)

#correct
#seqG %>% filter(culture == "Pro") %>% distinct(sample) %>% View()


#excludes Pro sequences
seqG <- seqG %>% filter(!str_detect(shortTaxa, "Prochlo"))

#exludes Synechococcaceae and Synechococcus sequences
seqG <- seqG %>% filter(!str_detect(shortTaxa, "Syn"))


seqG %>% filter(propSample == 0)
seqG %>% filter(propSample < .002)


#excludes 1223 because it doesn't have Pro, Synechococcaceae, or Synechococcus
#I figured this out in getSamplesWithoutProSyn.R
seqG <- seqG %>% filter(sample != "1223")


seqG %>% group_by(featureContig) %>% distinct(V1) %>% summarize(n = n()) %>% filter(n > 1)
seqG %>% group_by(featureContig) %>% distinct(sequence) %>% summarize(n = n()) %>% filter(n > 1)

#gets rid of variables so I can spread data
seqG <- seqG %>% select(-c(V1, sequence))
seqG <- seqG %>% select(-shortTaxa)


#spreads seqG so that the sequences that are not in a sample will have 0 abundance for that sample
seqSpread <- seqG %>% spread(key = featureContig, value = propSample, fill = 0)

dim(seqSpread)

#gathers seqG again, now with rows for sequences that are not present in a sample
seqG <- seqSpread %>% gather(3:237, key = "V1", value = "propSample")

seqG %>% distinct(sample) %>% nrow()

seqG %>% distinct(V1) %>% nrow()

seqG %>% group_by(V1) %>% distinct(sample) %>% summarize(n = n())

seqG %>% group_by(sample) %>% distinct(culture) %>% summarize(n = n())

set.seed(7)

#gets rid of sample variable so that I can run random forest
seqSpread <- seqSpread %>% select(-sample)

#makes column names not start with numbers
colnames(seqSpread) <- str_c("col", colnames(seqSpread))

#makes culture variable a factor so that I can run random forest
seqSpread$colculture <- factor(seqSpread$colculture)

set.seed(7)

#runs random forest
#changing mtry makes error rate worse
rf <- randomForest(colculture ~ ., data=seqSpread)

rf

seqSpread %>% filter(colculture == "Pro") %>% nrow()
seqSpread %>% filter(colculture == "Syn") %>% nrow()

set.seed(7)

model_sameSampleSize <- randomForest(colculture ~ ., data=seqSpread, sampsize = c(18,18))

model_sameSampleSize


#OOB estimate of  error rate: 8.11%
#Confusion matrix:
#  Pro Syn class.error
#Pro  55   1  0.01785714
#Syn   5  13  0.27777778


train <- sample(nrow(seqSpread), 0.5*nrow(seqSpread), replace = FALSE)
TrainSet <- seqSpread[train,]
ValidSet <- seqSpread[-train,]
set.seed(7)
trainModel <- randomForest(colculture ~ ., data = TrainSet, importance = TRUE)
trainModel
predTrain <- predict(trainModel, TrainSet, type = "class")
table(predTrain, TrainSet$colculture)  
predValid <- predict(trainModel, ValidSet, type = "class")
mean(predValid == ValidSet$colculture)                    
table(predValid,ValidSet$colculture)

importance(rf)

#makes dataframe of the importance of the variables
impDF <- importance(rf) %>% as.data.frame()

#makes variable for the sequence name
impDF <- impDF %>% mutate(V1 = rownames(impDF))

#gets rid of "col" at the beginning of sequence name
impDF <- impDF %>% mutate(V1 = str_replace(V1, "^col", ""))

#makes a variable for the sequence name without "contig..." at the end
impDF <- impDF %>% mutate(Feature.ID = str_replace(V1, "contig.*", ""))

#arranges impDF from most important to least important
#higher Mean Decrease in Gini indicates higher variable importance
impDF <- impDF %>% arrange(desc(MeanDecreaseGini))

nrow(impDF)

#adds taxonomy to impDf
impDF <- impDF %>% left_join(tax, by = c("Feature.ID"))

nrow(impDF)

impDF %>% filter(is.na(shortTaxa))

#gets the most important sequences from seqG
toPlot <- seqG %>% inner_join(impDF[1:20,], by = c("V1"))

toPlot %>% filter(is.na(shortTaxa))

toPlot %>% distinct(V1)

impDF[1:20,] %>% distinct(V1)

#makes combined variable for sequence ID and taxonomy
toPlot <- toPlot %>% mutate(featureTaxa = str_c(Feature.ID, shortTaxa, sep = ":\n"))

impDF <- impDF %>% mutate(featureTaxa = str_c(Feature.ID, shortTaxa, sep = ":\n"))

impDF[1:20,] %>% distinct(featureTaxa)



#makes dataframe of the importance of the variables
impDF_train <- importance(trainModel) %>% as.data.frame()

#makes variable for the sequence name
impDF_train <- impDF_train %>% mutate(V1 = rownames(impDF_train))

#gets rid of "col" at the beginning of sequence name
impDF_train <- impDF_train %>% mutate(V1 = str_replace(V1, "^col", ""))

#makes a variable for the sequence naem without "contig..." at the end
impDF_train <- impDF_train %>% mutate(Feature.ID = str_replace(V1, "contig.*", ""))

#arranges impDF_train from most important to least important
impDF_train <- impDF_train %>% arrange(desc(MeanDecreaseGini))

nrow(impDF_train)

#adds taxonomy to impDF_train
impDF_train <- impDF_train %>% left_join(tax, by = c("Feature.ID"))

nrow(impDF_train)

impDF_train %>% filter(is.na(shortTaxa))

#gets the most important sequences from seqG
toPlot_train <- seqG %>% inner_join(impDF_train[1:20,], by = c("V1"))

toPlot_train %>% filter(is.na(shortTaxa))

toPlot_train %>% distinct(V1)
toPlot %>% distinct(V1)

#15/20 of the sequences that are called the most important using model based on training 
#data set are also called the most important by the model based on the full dataset
toPlot_train %>% semi_join(toPlot, by = c("V1")) %>% distinct(V1)




importance(model_sameSampleSize)

#makes dataframe of the importance of the variables
impDF_sameSampleSize <- importance(model_sameSampleSize) %>% as.data.frame()

#makes variable for the sequence name
impDF_sameSampleSize <- impDF_sameSampleSize %>% mutate(V1 = rownames(impDF_sameSampleSize))

#gets rid of "col" at the beginning of sequence name
impDF_sameSampleSize <- impDF_sameSampleSize %>% mutate(V1 = str_replace(V1, "^col", ""))

#makes a variable for the sequence name without "contig..." at the end
impDF_sameSampleSize <- impDF_sameSampleSize %>% mutate(Feature.ID = str_replace(V1, "contig.*", ""))

#arranges impDF_sameSampleSize from most important to least important
#higher Mean Decrease in Gini indicates higher variable importance
impDF_sameSampleSize <- impDF_sameSampleSize %>% arrange(desc(MeanDecreaseGini))

nrow(impDF_sameSampleSize)

#adds taxonomy to impDF_sameSampleSize
impDF_sameSampleSize <- impDF_sameSampleSize %>% left_join(tax, by = c("Feature.ID"))

nrow(impDF_sameSampleSize)

impDF_sameSampleSize %>% filter(is.na(shortTaxa))

#gets the most important sequences from seqG
toPlot_sameSampleSize <- seqG %>% inner_join(impDF_sameSampleSize[1:20,], by = c("V1"))

#17/20 of the sequences that are called the most important using model based on same 
#sample sizes are also called the most important by the model based on different sample sizes
toPlot_sameSampleSize %>% semi_join(toPlot, by = c("V1")) %>% distinct(V1)

toPlot %>% distinct(V1)
toPlot_sameSampleSize %>% distinct(V1)


#16/20 of the sequences that are called the most important using model based on same 
#sample sizes are also called the most important by the model based on different sample sizes
#and the training dataset
toPlot_sameSampleSize %>% semi_join(toPlot_train, by = c("V1")) %>% distinct(V1)





