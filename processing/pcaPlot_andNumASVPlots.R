library(tidyverse)
library(readr)
library(readxl)
library(picante)
library(phyloseq)
library(randomForest)
library(ggfortify); library(ggplot2)
library(lubridate)

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

repSeqs %>% distinct(V2) %>% nrow()

seqG %>% distinct(sequence) %>% nrow()


nrow(repSeqs)

nrow(seqG)

#adds contig IDs to feature IDs in seqG
seqG <- seqG %>% left_join(repSeqs, by = c("sequence" = "V2"))

nrow(seqG)

seqG %>% filter(!is.na(V1)) %>% distinct(sequence) %>% nrow()
seqG %>% filter(is.na(V1))

#gets rid of ">" in feature IDs
seqG <- seqG %>% mutate(V1 = str_replace(V1, ">", ""))


#makes feature ID, contig ID combined variable that matches 
#the IDs used in the tree file
seqG <- seqG %>% mutate(featureContig = str_c(V1, sequence, sep = ""))

head(seqG)


#loads cleaned up taxonomy of the same length features
#corresponds to representative_sameLengthSeqs.fasta
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv")

tax %>% distinct(`Feature ID`) %>% nrow()

head(tax)

seqG %>% distinct(sample) %>% nrow()
seqG %>% distinct(sequence) %>% nrow()
seqG_withAllSamples_withProSyn <- seqG


nrow(seqG)

#adds taxonomy to seqG
seqG <- seqG %>% left_join(tax, by = c("V1" = "Feature ID"))

nrow(seqG)

seqG %>% filter(is.na(shortTaxa))
seqG <- seqG %>% filter(!is.na(shortTaxa))

seqG %>% distinct(sequence) %>% nrow()
seqG %>% distinct(sample) %>% nrow()

#JW3 only has Pro sequences so it gets excluded which is good because it only has 
#one sequence
seqG %>% filter(sample == "JW3")
seqG <- seqG %>% filter(sample != "JW3")


seqG %>% filter(sample == "1223")
seqG <- seqG %>% filter(sample != "1223")

seqG %>% distinct(sequence) %>% nrow()
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
seqG <- seqG %>% mutate(culture = ifelse(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916"), "Synechococcus", "Prochlorococcus"))

seqG %>% distinct(sample) %>% nrow()

#correct
seqG %>% filter(culture == "Synechococcus") %>% distinct(sample) %>% nrow()

#correct
seqG %>% filter(culture == "Prochlorococcus") %>% distinct(sample) %>% nrow()

seqG %>% filter(propSample == 0)
seqG %>% filter(propSample < .002)

seqG %>% group_by(featureContig) %>% distinct(V1) %>% summarize(n = n()) %>% filter(n > 1)
seqG %>% group_by(featureContig) %>% distinct(sequence) %>% summarize(n = n()) %>% filter(n > 1)

seqG %>% distinct(sample) %>% nrow()
seqG %>% distinct(sequence) %>% nrow()

head(seqG)

#gets rid of variables so I can spread data
seqG <- seqG %>% select(sample, culture, featureContig, propSample)

head(seqG)

#spreads seqG so that the sequences that are not in a sample will have 0 abundance for that sample
seqSpread <- seqG %>% spread(key = featureContig, value = propSample, fill = 0)

dim(seqSpread)

colnames(seqSpread)

#gathers seqG again, now with rows for sequences that are not present in a sample
seqG <- seqSpread %>% gather(3:237, key = "V1", value = "propSample")

head(seqG)

#correct number of samples and features
seqG %>% distinct(sample) %>% nrow()
seqG %>% distinct(V1) %>% nrow()

seqG %>% group_by(V1) %>% distinct(sample) %>% summarize(n = n()) %>% filter(n != 74)
seqG %>% group_by(sample) %>% distinct(culture) %>% summarize(n = n()) %>% filter(n != 1)

rownames(seqSpread) <- seqSpread$sample

#makes column names not start with numbers
colnames(seqSpread) <- str_c("col", colnames(seqSpread))

seqSpread %>% filter(colculture == "Prochlorococcus") %>% nrow()
seqSpread %>% filter(colculture == "Synechococcus") %>% nrow()


#loads cleaned up environmental variables for each culture
#this came from cleanUpEnvVariablesOfCulturesForPCA.R
env <- read_csv("cleanedUpCultureEnvVariablesForPCA.csv")

env %>% distinct(`METHOD (for isolation)`)

#these are the samples that have mismatched IDs or are not 
#present in cleanedUpCultureEnvVariablesForPCA.csv
env %>% anti_join(seqSpread, by = c("CULTURE_originalName" = "colsample"))
seqSpread %>% anti_join(env, by = c("colsample" = "CULTURE_originalName")) %>% distinct(colsample)


##from comparingSampleIDs.R
#1213 is 1313
#WH7803 is WH7803
#SYN1320 is SYN1220
#9201.2 is 9311

#AS9601 is correct ID rather than ASN9601
#SYN9503 is correct ID rather than S9503

seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "1213", "1313", colsample))
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "9201.2", "9311", colsample))
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "ASN9601", "AS9601", colsample))
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "9201.1", "9201", colsample))
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "S9503", "SYN9503", colsample))
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "SYN1320", "SYN1220", colsample))

#seqG$sample <- ifelse(seqG$sample == "8102", "SYNCLADEXVI", seqG$sample)
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "8102", "SYNCLADEXVI", colsample))

env %>% anti_join(seqSpread, by = c("CULTURE_originalName" = "colsample")) %>% nrow()
seqSpread %>% anti_join(env, by = c("colsample" = "CULTURE_originalName")) %>% nrow()

#excludes 1223 and JW3 from the env dataset
env %>% nrow()
seqSpread %>% nrow()

env$`DATE ISOLATED_Edited`
env %>% filter(is.na(`DATE ISOLATED_Edited`)) %>% nrow()

#makes a new variable that has the same month and day but has 2019 as the year 
env$dayMonth <- str_c("2019", str_extract(env$`DATE ISOLATED_Edited`, "-.*"))

env$dayMonth

env %>% filter(is.na(dayMonth)) %>% nrow()

#makes day month variable time format
env$dayMonth <- ymd(env$dayMonth, tz="UTC")

env %>% filter(is.na(dayMonth)) %>% nrow()

env$dayMonth

env %>% anti_join(seqSpread, by = c("CULTURE_originalName" = "colsample")) %>% nrow()
seqSpread %>% anti_join(env, by = c("colsample" = "CULTURE_originalName")) %>% nrow()

nrow(seqSpread)

#adds environmental variables to seqSpread
seqSpread <- seqSpread %>% left_join(env, by = c("colsample" = "CULTURE_originalName"))

nrow(seqSpread)

colnames(seqSpread)

env %>% distinct(`METHOD (for isolation)`)
seqSpread %>% distinct(`METHOD (for isolation)`)


#gets just the columns with relative abundances of sequences
df <- seqSpread[c(3:237)]

#replaces spaces with "_" in column names
colnames(seqSpread) <- str_replace_all(colnames(seqSpread), " ", "_")

df_forVector <- df
colnames(df_forVector)
colnames(df_forVector) <- str_extract(colnames(df_forVector), "contig.*")
colnames(df_forVector)
dim(df_forVector)

#this is so I can get the vector labels
autoplot(prcomp(df_forVector), data = seqSpread, colour = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3) + labs(colour = "sample") 

ggsave("pca_vectorLabels.svg", dpi = 600, height = 6, width = 8)


###do I have to check if these variables are correlated??

autoplot(prcomp(df), data = seqSpread, colour = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3) + labs(colour = "Culture")

ggsave("pca_culture.png", dpi = 600, height = 6, width = 8)

autoplot(prcomp(df), data = seqSpread, colour = 'ECOTYPE', shape = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3) + labs(shape = "Culture", colour = "Ecotype")

ggsave("pca_ecotype.png", dpi = 600, height = 6, width = 8)

cladeOrder <- seqSpread %>% arrange(ECOTYPE, CLADE) %>% distinct(CLADE)
cladeOrder
cladeOrder <- cladeOrder$CLADE
cladeOrder

seqSpread$CLADE <- factor(seqSpread$CLADE, levels = cladeOrder)
levels(seqSpread$CLADE)

autoplot(prcomp(df), data = seqSpread, colour = 'CLADE', shape = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3) + labs(shape = "Culture", colour = "Clade")

ggsave("pca_clade.png", dpi = 600, height = 6, width = 8)

autoplot(prcomp(df), data = seqSpread, colour = 'PLACE_OF_ORIGIN', shape = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3) + labs(shape = "Culture", colour = "Place of origin")

ggsave("pca_originLocation.png", dpi = 600, height = 6, width = 8)

str_extract(seqSpread$CRUISE, "EqPac/IRONEX\\?")

seqSpread$CRUISE <- str_replace(seqSpread$CRUISE, "EqPac/IRONEX\\?", "Possibly EqPac/IRONEX")

seqSpread %>% distinct(CRUISE)

autoplot(prcomp(df), data = seqSpread, colour = 'CRUISE', shape = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3) + labs(shape = "Culture", colour = "Cruise")

ggsave("pca_cruise.png", dpi = 600, height = 6, width = 8)

autoplot(prcomp(df), data = seqSpread, colour = 'DEPTH', shape = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3) + labs(shape = "Culture", colour = "Depth") + 
  scale_colour_distiller(trans = 'reverse')

ggsave("pca_depth.png", dpi = 600, height = 6, width = 8)

autoplot(prcomp(df), data = seqSpread, colour = 'ISOLATOR', shape = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3)  + labs(shape = "Culture", colour = "Isolator") + 
  theme(legend.text = element_text(size = 4)) 

ggsave("pca_isolator.png", dpi = 600, height = 6, width = 8)

colnames(seqSpread)[247]
colnames(seqSpread)[247] <- "isolationMethod"

seqSpread %>% distinct(isolationMethod)

autoplot(prcomp(df), data = seqSpread, colour = 'isolationMethod', shape = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3)  + labs(shape = "Culture", colour = "Isolation method")

ggsave("pca_isolationMethod.png", dpi = 600, height = 6, width = 8)


#seqSpread %>% distinct(Notes) %>% View()

seqSpread <- seqSpread %>% mutate(Notes = ifelse(Notes == "The seawater was amended with nitrogen, phosphorous and trace metals (PRO2 nutrient additions, except all nitrogen sources were replaced by 0.217 mM sodium nitrate)",
                                                 "The seawater was amended with nitrogen, phosphorous\nand trace metals (PRO2 nutrient additions,\nexcept all nitrogen sources were replaced by\n0.217 mM sodium nitrate)", Notes))


seqSpread %>% filter(str_detect(Notes, "was diluted in PRO3V m")) %>% distinct(Notes)

seqSpread <- seqSpread %>% mutate(Notes = ifelse(str_detect(Notes, "was diluted in PRO3V m"),
                                                 "The seawater used for isolations was first filtered\nthrough a 1 μm filter with no amendments and kept in\nthe dark at 18–20 °C for 21 days. The total red fluorescing\nphytoplankton population (1×105 cells ml−1determined with a\nGuava EasyCyte flow cytometer) was diluted in PRO3V media37\nmade with the same South Atlantic water that had been\nfiltered through a 0.1 μm Supor 142 mm filter, then\nautoclaved to sterilize. This media contained 100 μM NH4Cl,\n10 μM NaH2PO4, PRO2 trace metals37 and f/2 vitamins (0.1 μg l−1\ncyanocobalamin, 20 g l−1thiamin and 1 μg l−1 biotin38,39).\nTen cells were dispensed into 1 ml volumes in a\n48-well polystyrene multiwell culture plate and\nincubated at 20 °C in ~20 μmol Q m−2 s−1 (14:10 light:dark)\nfor 2 months.",
                                                 Notes))

seqSpread %>% filter(str_detect(Notes, "was diluted in PRO3V m")) %>% distinct(Notes)


#seqSpread %>% distinct(Notes) %>% View()

autoplot(prcomp(df), data = seqSpread, colour = 'Notes', shape = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3) + labs(shape = "Culture", colour = "Notes") + theme(legend.text = element_text(size = 4)) 

ggsave("pca_notes.png", dpi = 600, height = 6, width = 8)

#remember that if the day wasn't available, I made it the first of the month
#and if neither the month or day were available, I made it the first day of the year
autoplot(prcomp(df), data = seqSpread, colour = 'DATE_ISOLATED_Edited', shape = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3) + labs(shape = "Culture", colour = "Isolation date")

ggsave("pca_dateIsolated.png", dpi = 600, height = 6, width = 8)


#seqSpread %>% filter(dayMonth == ymd("2019-01-01"))

#remember that if the day wasn't available, I made it the first of the month
#and if neither the month or day were available, I made it the first day of the year
autoplot(prcomp(df), data = seqSpread, colour = 'dayMonth', shape = 'colculture',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = FALSE, loadings.label.size = 3) + labs(colour = "Month, day of isolation", shape = "Culture")

ggsave("pca_dayMonthIsolated.png", dpi = 600, height = 6, width = 8)

colnames(seqSpread)
colnames(seqSpread)[3:237]

seqGEnv <- seqSpread %>% gather(3:237, key = "sequence", value = "propSample")

head(seqGEnv)

seqGEnv %>% group_by(colsample) %>% summarize(sumPropSample = sum(propSample), numSeqs = n()) %>% filter(numSeqs != 235)

seqGEnv %>% distinct(colsample) %>% nrow()
seqGEnv %>% distinct(sequence) %>% nrow()

colnames(seqGEnv)

#there is only one isolation method for each sample
seqGEnv %>% group_by(colsample) %>% distinct(isolationMethod) %>% summarize(n = n()) %>% filter(n != 1)

#gets the number of nonPro, nonSyn sequences in each sample
summ <- seqGEnv %>% group_by(colsample, isolationMethod) %>% filter(propSample > 0) %>% summarize(numSeqs = n())

head(summ)

summ %>% group_by(isolationMethod) %>% summarize(meanNumSeqs = mean(numSeqs), sd = sd(numSeqs), n = n())

summ %>% filter(is.na(isolationMethod)) %>% nrow()

summ <- summ %>% mutate(isolationMethod = ifelse(is.na(isolationMethod), "NA", isolationMethod))

summ %>% group_by(isolationMethod) %>% summarize(meanNumSeqs = mean(numSeqs), sd = sd(numSeqs), n = n())
summ %>% filter(isolationMethod == "NA") %>% nrow()

#gets list of the samples, arranged from most seqs to least seqs 
order <- summ %>% group_by(isolationMethod) %>% summarize(meanNumSeqs = mean(numSeqs), sd = sd(numSeqs)) %>% arrange(desc(meanNumSeqs))
order
order <- order %>% ungroup() %>% distinct(isolationMethod)
order
order <- order %>% as.list()
order
order <- order$isolationMethod
order

#makes the isolation variable a factor, in order of number of seqs
summ$isolationMethod <- factor(summ$isolationMethod, levels = order)

levels(summ$isolationMethod)

#plots the mean number of nonPro, nonSyn sequences by isolation method
summ %>% group_by(isolationMethod) %>% summarize(meanNumSeqs = mean(numSeqs), sd = sd(numSeqs)) %>% 
  ggplot() + geom_bar(aes(x = isolationMethod, y = meanNumSeqs), stat = 'identity') + 
  geom_errorbar(aes(x = isolationMethod, ymin=meanNumSeqs-sd, ymax=meanNumSeqs+sd), width=.2, position=position_dodge(0.05)) + 
  labs(x = "Isolation method", y = "Mean number of ASVs in sample, excluding Pro. and Syn.") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("meanNumberOfASVsByIsolationMethod.png", dpi = 600, height = 10, width = 10)
ggsave("meanNumberOfASVsByIsolationMethod.svg", dpi = 600, height = 10, width = 10)



head(seqGEnv)

#these are the sample IDs that I fixed
seqGEnv %>% filter(colsample != CULTURE) %>% distinct(colsample, CULTURE)

#gets the number of nonPro, nonSyn sequences in each sample
seqGEnv %>% group_by(colsample, colculture) %>% filter(propSample > 0) %>% summarize(numSeqs = n()) %>% arrange(numSeqs)

seqGEnv %>% distinct(colculture)

toPlot <- seqGEnv %>% group_by(colsample, colculture) %>% filter(propSample > 0) %>% summarize(numSeqs = n()) %>% 
  mutate(colculture = str_replace(colculture, "Prochlorococcus", "Pro. sample")) %>% 
  mutate(colculture = str_replace(colculture, "Synechococcus", "Syn. sample"))

toPlot %>% nrow()
toPlot %>% distinct(colsample) %>% nrow()

#gets list of the samples, arranged from most seqs to least seqs 
orderNumSeqs <- toPlot %>% arrange(desc(numSeqs))
orderNumSeqs <- orderNumSeqs %>% ungroup() %>% distinct(colsample)
orderNumSeqs <- orderNumSeqs %>% as.list()
orderNumSeqs <- orderNumSeqs$colsample

#makes the samples variable a factor, in order of number of seqs
toPlot$colsample <- factor(toPlot$colsample, levels = orderNumSeqs)

write_csv(toPlot, "numNonProNonSynSeqsBySample.csv")

#plots the number of nonPro, nonSyn sequences in each sample, color coding by whether it is a Pro or Syn culture  
toPlot %>% ggplot(aes(x = colsample, y = numSeqs, fill = colculture)) + geom_bar(stat = 'identity', width = .5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "Number of ASVs, excluding Pro. and Syn.", x = "Sample", 
                                                                  fill = "")                                                                 

ggsave("numNonProNonSynSeqsBySample.png", dpi = 600, height = 6, width = 10)
ggsave("numNonProNonSynSeqsBySample.svg", dpi = 600, height = 6, width = 10)



#loads tree of all Pro, Syn sequences, not just ones that are dominant in a culture
#these sequences are redundant because the same Pro, Syn ASV can be found in multiple samples
#from analysisOnlyProSynASVs.txt (checked)
proSynASVTree <- read.tree("exported-rooted-tree_rep-seqs-onlyProSynSequences_noContamination_noMock_no1223/tree.nwk")

plot(proSynASVTree)
proSynASVTree$tip.label

seqG_withAllSamples_withProSyn %>% distinct(sequence) %>% nrow()
seqG_withAllSamples_withProSyn %>% distinct(sample) %>% nrow()

#exclues JW3 and 1223
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% filter(sample != "JW3")
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% filter(sample != "1223")

#includes Pro, Syn sequences
seqG_withAllSamples_withProSyn %>% distinct(sequence) %>% nrow()
seqG_withAllSamples_withProSyn %>% distinct(sample) %>% nrow()

head(seqG_withAllSamples_withProSyn)

#gets rid of unnecessary variables
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% select(-c(V1, featureContig))

#loads rep seqs from analysisOnlyProSynASVs.txt (checked)
#this corresponds to "exported-rooted-tree_rep-seqs-onlyProSynSequences_noContamination_noMock_no1223/tree.nwk
repSeqs_onlyProSyn <- read.table("rep-seqs-onlyProSynSequences_noContamination_noMock_no1223.fasta", fill = TRUE)

head(repSeqs_onlyProSyn)

#gets just rows with IDs
repSeqs_onlyProSyn <- repSeqs_onlyProSyn %>% filter(str_detect(V1, ">"))

#only pro and syn features
#not just ones that are dominant in a culture
nrow(repSeqs_onlyProSyn)

head(repSeqs_onlyProSyn)

#all of the IDs are distinct
repSeqs_onlyProSyn %>% distinct(V1) %>% nrow()
repSeqs_onlyProSyn %>% distinct(V2) %>% nrow()

head(seqG_withAllSamples_withProSyn)

#adds rep seq IDs to seqG_withAllSamples_withProSyn)
nrow(seqG_withAllSamples_withProSyn)
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% left_join(repSeqs_onlyProSyn, by = c("sequence" = "V2"))
nrow(seqG_withAllSamples_withProSyn)

head(seqG_withAllSamples_withProSyn)

#the non-Pro, non-Syn sequences
seqG_withAllSamples_withProSyn %>% filter(is.na(V1))

#there are 33 Pro, Syn sequences
#not just ones that are dominant in a culture
seqG_withAllSamples_withProSyn %>% filter(!is.na(V1)) %>% distinct(sequence) %>% nrow()
seqG_withAllSamples_withProSyn %>% filter(!is.na(V1)) %>% distinct(V1) %>% nrow()

#gets rid of non-Pro, non-Syn sequences
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% filter(!is.na(V1))
seqG_withAllSamples_withProSyn %>% distinct(V1) %>% nrow()

seqG_withAllSamples_withProSyn %>% distinct(sample) %>% nrow()

#there are more than 74 rows because a couple of the cultures have more than one 
#Pro, Syn sequence
seqG_withAllSamples_withProSyn %>% nrow()

#gets the cultures have more than one Pro, Syn sequence
mult <- seqG_withAllSamples_withProSyn %>% group_by(sample) %>% summarize(n = n()) %>% filter(n > 1)

#there is a dominant Pro or Syn ASV in each of the samples
#seqG_withAllSamples_withProSyn %>% semi_join(mult, by = c("sample")) %>% arrange(sample) %>% View()

#seqG_withAllSamples_withProSyn %>% group_by(sample) %>% arrange(desc(propSample)) %>% View()

seqG_withAllSamples_withProSyn %>% distinct(sequence) %>% nrow()

#gets the dominant Pro or Syn ASV in each of the samples
#seqG_withAllSamples_withProSyn %>% group_by(sample) %>% arrange(desc(propSample)) %>% View()
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% group_by(sample) %>% arrange(desc(propSample)) %>% slice(1)

seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% ungroup()

#seqG_withAllSamples_withProSyn %>% semi_join(mult, by = c("sample")) %>% arrange(sample) %>% View()

#now there is only one Pro, Syn sequence per sample
seqG_withAllSamples_withProSyn %>% nrow()

seqG_withAllSamples_withProSyn %>% 
  distinct(sequence) %>% nrow()

seqG_withAllSamples_withProSyn %>% 
  distinct(sequence) %>% 
  write.table("onlyDominantProSynSequences_noContamination_noMock_no1223.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


seqG_withAllSamples_withProSyn %>% nrow()

seqG_withAllSamples_withProSyn %>% 
  select(sequence) %>% 
  write.table("onlyDominantProSynSequences_noContamination_noMock_no1223_redundant.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#from analysisOnlyProSynASVs.txt (checked)
proSynASVTree <- read.tree("exported-rooted-tree_rep-seqs-onlyDominantProSynSequences_noContamination_noMock_no1223_redundant/tree.nwk")

plot(proSynASVTree)

#gets order of tree leaves
order <- proSynASVTree$tip.label
order <- order %>% as.data.frame()
colnames(order) <- "featureSequence"

nrow(order)
head(order)

seqG_withAllSamples_withProSyn %>% distinct(sequence) %>% nrow()
seqG_withAllSamples_withProSyn %>% nrow()

#arrange sequences in seqG_withAllSamples_withProSyn so that I can rename them 
#to fix tree labels
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% arrange(sequence)

#fix tree labels in seqG_withAllSamples_withProSyn to match tree labels
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% mutate(seq2 = ifelse(lag(sequence) == sequence, str_c(sequence, "_2"), sequence))
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% mutate(seq2 = ifelse(lag(seq2) == seq2, str_replace(seq2, "_2$", "_3"), seq2))
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% mutate(seq2 = ifelse(lag(seq2) == seq2, str_replace(seq2, "_3$", "_4"), seq2))
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% mutate(seq2 = ifelse(lag(seq2) == seq2, str_replace(seq2, "_4$", "_5"), seq2))
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% mutate(seq2 = ifelse(lag(seq2) == seq2, str_replace(seq2, "_5$", "_6"), seq2))
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% mutate(seq2 = ifelse(lag(seq2) == seq2, str_replace(seq2, "_6$", "_7"), seq2))
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% mutate(seq2 = ifelse(lag(seq2) == seq2, str_replace(seq2, "_7$", "_8"), seq2))
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% mutate(seq2 = ifelse(lag(seq2) == seq2, str_replace(seq2, "_8$", "_9"), seq2))
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% mutate(seq2 = ifelse(lag(seq2) == seq2, str_replace(seq2, "_9$", "_10"), seq2))
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% mutate(seq2 = ifelse(lag(seq2) == seq2, str_replace(seq2, "_10$", "_11"), seq2))

seqG_withAllSamples_withProSyn %>% filter(is.na(seq2))

#fills in missing fixed sequence IDs
seqG_withAllSamples_withProSyn <- seqG_withAllSamples_withProSyn %>% 
  mutate(seq2 = ifelse(sequence == "contig_1", c("contig_1", "contig_1_2", "contig_1_3", "contig_1_4", "contig_1_5", "contig_1_6",
                  "contig_1_7", "contig_1_8", "contig_1_9", "contig_1_10", "contig_1_11"), seq2))

#all of the fixed sequences are unique
seqG_withAllSamples_withProSyn %>% group_by(seq2) %>% summarize(n = n()) %>% filter(n > 1)

seqG_withAllSamples_withProSyn %>% filter(is.na(seq2))

head(order) 

#adds samples Pro, Syn sequences are in to the dataframe of the order of the 
#sequences in the tree
nrow(order)
order <- order %>% left_join(seqG_withAllSamples_withProSyn, by = c("featureSequence" = "seq2"))
nrow(order)

order %>% filter(is.na(sample))

str(toPlot$colsample)

#makes sample variable a character
toPlot$colsample <- as.character(toPlot$colsample)

str(toPlot$colsample)

#these are the sample IDs that don't match between toPlot and order
toPlot %>% anti_join(order, by = c("colsample" = "sample"))

#fixes sample IDs
order <- order %>% mutate(sample = ifelse(sample == "1213", "1313", sample))
order <- order %>% mutate(sample = ifelse(sample == "9201.2", "9311", sample))
order <- order %>% mutate(sample = ifelse(sample == "ASN9601", "AS9601", sample))
order <- order %>% mutate(sample = ifelse(sample == "9201.1", "9201", sample))
order <- order %>% mutate(sample = ifelse(sample == "S9503", "SYN9503", sample))
order <- order %>% mutate(sample = ifelse(sample == "SYN1320", "SYN1220", sample))

#seqG$sample <- ifelse(seqG$sample == "8102", "SYNCLADEXVI", seqG$sample)
order <- order %>% mutate(sample = ifelse(sample == "8102", "SYNCLADEXVI", sample))

#now all of the sample IDs match
toPlot %>% anti_join(order, by = c("colsample" = "sample"))

#makes a list of the samples in the order they appear in the phylogeny
order$sample
sampleOrder <- order$sample
length(sampleOrder)
sampleOrder

#puts samples in toPlot in the same order as in the tree
toPlot$colsample <- factor(toPlot$colsample, levels = sampleOrder)

unique(order$featureSequence == proSynASVTree$tip.label)

proSynASVTree$tip.label
order %>% select(featureSequence, sample) %>% head()

levels(toPlot$colsample)

toPlot %>% ungroup() %>% filter(is.na(colsample))

#plots the number of nonPro, nonSyn sequences in each sample, color coding by whether it is a Pro or Syn culture  
toPlot %>% ggplot(aes(x = colsample, y = numSeqs, fill = colculture)) + geom_bar(stat = 'identity', width = .5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "Number of ASVs, excluding Pro. and Syn.", x = "Sample", 
                                                                  fill = "")    

ggsave("numNonProNonSynSeqsBySample_inOrderOfProSynASVPhylogeny.png", dpi = 600, height = 6, width = 10)
ggsave("numNonProNonSynSeqsBySample_inOrderOfProSynASVPhylogeny.svg", dpi = 600, height = 6, width = 10)

plot(proSynASVTree)

##fixes tree leaves to be sample IDs
head(order)
order <- order %>% select(featureSequence, sample)

unique(proSynASVTree$tip.label == order$featureSequence)

proSynASVTree_edited <- proSynASVTree

proSynASVTree_edited$tip.label <- order$sample

plot(proSynASVTree_edited)

write.tree(proSynASVTree_edited, "rooted-tree_onlyDominantProSynSequences_redundant.nwk")



