library(tidyverse)
library(readr)
library(readxl)
library(picante)
library(phyloseq)
library(randomForest)
library(ggfortify); library(ggplot2)
library(lubridate)
library(phyloseq)
library(ggdendro)
library(ComplexHeatmap)
library(gridExtra)
library(cowplot)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads corrected counts of each taxa in each sample 
#this excludes sequences that were determined to be contamination from nearby wells
seq <- read.table("seqtab_corrected.txt")

numCols <- ncol(seq)

cols <- str_c("col", 1:numCols)
cols

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
head(seqG)

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

str(seqG)

#gets rid of the rows for which a sequence is not present in a sample
#this excludes the sequences that are just present in the unrelated samples (Con|Epi|Pcn|Lin)
seqG <- seqG %>% filter(abundance > 0)

seqG %>% distinct(sequence) %>% nrow()

#gets the proportion that each sequence is in each sample
seqG <- seqG %>% mutate(propSample = abundance/totalSeqs)

seqG %>% distinct(sequence) %>% nrow() 

str(seqG)

#excludes the sequences in samples that have a low relative abundances in the sample
seqG <- seqG %>% filter(propSample >= 0.002)

seqG %>% distinct(sequence) %>% nrow()

#gets rid of unnecessary variables
seqG <- seqG %>% select(-c(abundance, totalSeqs))

head(seqG)

#makes sequence names in seqG match sequence names in repSeqs
seqG$sequence <- str_replace(seqG$sequence, "col", "contig_")

head(seqG)

#loads the representative sequences fasta file so that I can get the 
#contig ID assigned to each feature ID
#this is with excluding Pro and Syn
repSeqs <- read.table("representative_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.fasta", fill = TRUE)

head(repSeqs)

#gets just the rows in repSeqs that correspond to IDs
repSeqs <- repSeqs %>% filter(str_detect(V1, ">"))

head(repSeqs)
repSeqs %>% distinct(V2) %>% nrow()

head(seqG)
seqG %>% distinct(sequence) %>% nrow()

nrow(repSeqs)

nrow(seqG)

#adds contig IDs to feature IDs in seqG
seqG <- seqG %>% left_join(repSeqs, by = c("sequence" = "V2"))

nrow(seqG)

head(seqG)

seqG %>% filter(!is.na(V1)) %>% distinct(sequence) %>% nrow()
seqG %>% filter(is.na(V1))

seqG %>% filter(!is.na(V1)) %>% head()

nrow(seqG)

#gets rid of ">" in feature IDs
seqG <- seqG %>% mutate(V1 = str_replace(V1, ">", ""))

nrow(seqG)

seqG %>% filter(!is.na(V1)) %>% head()

#makes feature ID, contig ID combined variable that matches 
#the IDs used in the tree file
seqG <- seqG %>% mutate(featureContig = str_c(V1, sequence, sep = ""))

seqG %>% filter(!is.na(V1)) %>% head()

nrow(seqG)


#loads cleaned up taxonomy of the same length features
#corresponds to representative_sameLengthSeqs.fasta
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv")

head(tax)

tax %>% distinct(`Feature ID`) %>% nrow()

seqG %>% filter(!is.na(V1)) %>% head()

nrow(seqG)

#adds taxonomy to seqG
seqG <- seqG %>% left_join(tax, by = c("V1" = "Feature ID"))

nrow(seqG)

seqG %>% filter(!is.na(V1)) %>% head()

#gets rid of pro and syn sequences
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

#correct number of sequences and samples after excluding Pro and Syn sequences
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

seqG %>% distinct(sample) %>% nrow()

#makes a variable for whether sample is Pro or Syn
seqG <- seqG %>% mutate(culture = ifelse(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916"), "Synechococcus", "Prochlorococcus"))

seqG %>% distinct(sample) %>% nrow()

#correct
seqG %>% filter(culture == "Synechococcus") %>% distinct(sample) %>% nrow()

#correct
seqG %>% filter(culture == "Prochlorococcus") %>% distinct(sample) %>% nrow()

seqG %>% filter(propSample == 0)
seqG %>% filter(propSample < .002)

head(seqG)

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

head(env)

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

#the entries that I don't know the month, day, or both still have a dayMonth value

#makes a variable for the month that cultures were isolated
env <- env %>% mutate(month = str_extract(`DATE ISOLATED_Edited`, "(?<=-)[0-9]*(?=-)"))

env$month

env %>% select(`DATE ISOLATED_Edited`, dayMonth, month)

env %>% distinct(month)
env %>% filter(is.na(month))

#all I know is that WH6501 is from 1965
#the culture website and excel sheet say MED4 and MED1 were isolated Jan. 1989
env %>% filter(month == "01")

#makes the month NA for WH6501
env %>% filter(is.na(month)) %>% nrow()
env <- env %>% mutate(month = ifelse(CULTURE == "WH6501", NA, month))
env %>% filter(is.na(month)) %>% nrow()

env %>% filter(month == "01")
env %>% filter(CULTURE == "WH6501")

env %>% distinct(month)

#makes a variable with the month name
env <- env %>% mutate(monthNew = ifelse(month == "01", "Jan.", month))
env <- env %>% mutate(monthNew = ifelse(month == "02", "Feb.", monthNew))
env <- env %>% mutate(monthNew = ifelse(month == "03", "Mar.", monthNew))
env <- env %>% mutate(monthNew = ifelse(month == "04", "Apr.", monthNew))
env <- env %>% mutate(monthNew = ifelse(month == "05", "May.", monthNew))
env <- env %>% mutate(monthNew = ifelse(month == "06", "Jun.", monthNew))
env <- env %>% mutate(monthNew = ifelse(month == "07", "Jul.", monthNew))
env <- env %>% mutate(monthNew = ifelse(month == "08", "Aug.", monthNew))
env <- env %>% mutate(monthNew = ifelse(month == "09", "Sep.", monthNew))
env <- env %>% mutate(monthNew = ifelse(month == "10", "Oct.", monthNew))
env <- env %>% mutate(monthNew = ifelse(month == "11", "Nov.", monthNew))
env <- env %>% mutate(monthNew = ifelse(month == "12", "Dec.", monthNew))

env %>% filter(is.na(month)) %>% nrow()
env %>% distinct(month, monthNew) %>% arrange(month)
env %>% filter(is.na(monthNew)) %>% nrow()

env %>% anti_join(seqSpread, by = c("CULTURE_originalName" = "colsample")) %>% nrow()
seqSpread %>% anti_join(env, by = c("colsample" = "CULTURE_originalName")) %>% nrow()

nrow(seqSpread)

#adds environmental variables to seqSpread
seqSpread <- seqSpread %>% left_join(env, by = c("colsample" = "CULTURE_originalName"))

nrow(seqSpread)

colnames(seqSpread)

env %>% distinct(`METHOD (for isolation)`)
seqSpread %>% distinct(`METHOD (for isolation)`)

dim(seqSpread)

###phylosift

#https://joey711.github.io/phyloseq/import-data.html#phyloseq-ize_data_already_in_r

#phyloseq - Takes as argument an otu_table and any unordered list of valid phyloseq components: 
#sample_data, tax_table, phylo, or XStringSet. The tip labels of a phylo-object (tree) must match 
#the OTU names of the otu_table, and similarly, the sequence names of an XStringSet object must match 
#the OTU names of the otu_table

colnames(seqSpread)

phylo_seqSpread <- seqSpread[c(1:237)]

dim(phylo_seqSpread)
colnames(phylo_seqSpread)

#makes the rownames in seqSpread the sample
rownames(phylo_seqSpread) <- phylo_seqSpread$colsample

#gets rid of the sample variable in seqSpread
phylo_seqSpread <- phylo_seqSpread %>% select(-c(colsample, colculture))

colnames(phylo_seqSpread)
colnames(phylo_seqSpread) <- str_replace(colnames(phylo_seqSpread), "col", "")
colnames(phylo_seqSpread)

#doesn't include Pro or Syn
#correct dimensions
dim(phylo_seqSpread)


#makes phyloseq object out of otu table
otu <- otu_table(phylo_seqSpread, taxa_are_rows = FALSE)



#loads cleaned up taxonomy of the same length features, no pro, no syn, no contamination, no mock samples, no 1223 or JW3
#this is from makeTreeNodes_noContamination_noMock_noProSyn_without1223.R
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv")

head(tax)

tax %>% distinct(`Feature ID`) %>% nrow()
tax %>% nrow()

tax <- tax %>% mutate(featureContig = str_c(`Feature ID`, V1))

head(tax)
tax$featureContig

tax %>% filter(is.na(Taxon))

#there are no Pro or Syn features in tax
tax %>% filter(str_detect(Taxon, "Prochlorococcus")) %>% distinct(Taxon)
tax %>% filter(str_detect(Taxon, "Syn")) %>% distinct(Taxon)

#makes dataframe with the taxonomy for each feature
taxPhyloseq <- tax %>% distinct(featureContig, shortTaxa)

nrow(tax)
nrow(taxPhyloseq)

head(taxPhyloseq)

#makes rownames of taxPhyloseq the feature ID
rownames(taxPhyloseq) <- taxPhyloseq$featureContig

#gets rid of feature ID variable
taxPhyloseq <- taxPhyloseq %>% select(-featureContig)

nrow(taxPhyloseq)

#makes phyloseq object out of taxPhyloseq
TAX <- tax_table(as.matrix(taxPhyloseq))

#loads rooted tree
#for this, I excluded sample 1223
tree <- read_tree("~/Dropbox (MIT)/Sean16SSep2019/exported-rooted-treeSameLength_noContamination_noMock_noProSyn_withoutSample1223/tree.nwk")


otuDF <- colnames(otu) %>% as.data.frame()

treeDF <- tree$tip.label %>% as.data.frame()

TAXDF <- TAX %>% rownames() %>% as.data.frame()


otuDF %>% nrow()

treeDF %>% nrow()

TAXDF %>% nrow()


colnames(otuDF) <- "V1"

colnames(treeDF) <- "V1"

colnames(TAXDF) <- "V1"



treeDF %>% anti_join(otuDF, by = c("V1"))
TAXDF %>% anti_join(otuDF, by = c("V1"))

otuDF %>% anti_join(treeDF, by = c("V1"))
TAXDF %>% anti_join(treeDF, by = c("V1"))

otuDF %>% anti_join(TAXDF, by = c("V1"))
treeDF %>% anti_join(TAXDF, by = c("V1"))

taxa_names(otu)
taxa_names(TAX)
taxa_names(tree)

#makes combined phyloseq object out of otu, TAX, tree
physeq <- phyloseq(otu, TAX, tree)

#https://groups.google.com/forum/#!topic/qiime-forum/7O8Bpz6T0VU
ordu <- ordinate(physeq, "PCoA", "unifrac", weighted = FALSE)
plot_ordination(physeq, ordu)
#plot_ordination(physeq, ordu, color = "SampleType", shape = "human")

set.seed(7)

#there is a warning
unifracRes <- UniFrac(physeq, weighted = FALSE)

#makes tree
sampleTree <- unifracRes %>% hclust()

#plot_heatmap(physeq, taxa.label="shortTaxa")

png("unifracTreeClusteringForMetadataHeatmap_nodesDoNotLineUp.png", 800, 800)
#par(mfrow=c(1,2))
plot(sampleTree)
dev.off()

#from https://stackoverflow.com/questions/14118033/horizontal-dendrogram-in-r-with-labels
#convert cluster object to use with ggplot
dendr <- dendro_data(sampleTree, type="rectangle") 

#from https://stackoverflow.com/questions/14118033/horizontal-dendrogram-in-r-with-labels
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank()) + theme(legend.position = "none")

ggsave("unifracTreeClusteringForMetadataHeatmap.png", height = 14, width = 14, dpi = 600)


sampleTree$labels[sampleTree$order]

#puts env sample variable in the same order as unifrac clustering of samples
env$CULTURE_originalName <- factor(env$CULTURE_originalName, levels = sampleTree$labels[sampleTree$order])
levels(env$CULTURE_originalName)

length(levels(env$CULTURE_originalName))
env %>% head()
env %>% select(CULTURE) %>% head()
env %>% select(`ALTERNATIVE NAME`) %>% head()

str(env)
env <- env %>% as.data.frame()
str(env)

#cultureType
colnames(env)
env %>% distinct(cultureType)
cultureType <- env %>% select(CULTURE_originalName, cultureType) %>% 
  gather(cultureType, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "cultureType", "Culture")) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Culture") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())



#ecotype
env %>% distinct(ECOTYPE)
env %>% arrange(ECOTYPE) %>% distinct(ECOTYPE)
#env <- env %>% mutate(ECOTYPE = ifelse(is.na(ECOTYPE), "NA", ECOTYPE))
#env %>% distinct(ECOTYPE)

order <- env %>% arrange(ECOTYPE) %>% distinct(ECOTYPE)
order <- order %>% select(ECOTYPE)
order <- order %>% as.list()
order <- order$ECOTYPE
order
env$ECOTYPE <- factor(env$ECOTYPE, levels = order, exclude = NULL)
levels(env$ECOTYPE)

ecotype <- env %>% select(CULTURE_originalName, ECOTYPE) %>% 
  gather(ECOTYPE, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "ECOTYPE", "Ecotype")) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Ecotype") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())


#clade
colnames(env)
env %>% distinct(CLADE)
#env <- env %>% mutate(CLADE = ifelse(is.na(CLADE), "NA", CLADE))
#env %>% distinct(CLADE)

order <- env %>% arrange(CLADE) %>% distinct(CLADE)
order <- order %>% select(CLADE)
order <- order %>% as.list()
order <- order$CLADE
order

#env <- env %>% mutate(CLADE = ifelse(is.na(CLADE), "NA", CLADE))
#env %>% distinct(CLADE)
env$CLADE <- factor(env$CLADE, levels = order, exclude = NULL)
levels(env$CLADE)

clade <- env %>% select(CULTURE_originalName, CLADE) %>% 
  gather(CLADE, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "CLADE", "Clade")) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Clade") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())


#place of origin
colnames(env)
env %>% distinct(`PLACE OF ORIGIN`)
#env <- env %>% mutate(`PLACE OF ORIGIN` = ifelse(is.na(`PLACE OF ORIGIN`), "NA", `PLACE OF ORIGIN`))
#env %>% distinct(`PLACE OF ORIGIN`)

order <- env %>% arrange(`PLACE OF ORIGIN`) %>% distinct(`PLACE OF ORIGIN`)
order <- order %>% select(`PLACE OF ORIGIN`)
order <- order %>% as.list()
order <- order$`PLACE OF ORIGIN`
order

#env <- env %>% mutate(`PLACE OF ORIGIN` = ifelse(is.na(`PLACE OF ORIGIN`), "NA", `PLACE OF ORIGIN`))
env$`PLACE OF ORIGIN` <- factor(env$`PLACE OF ORIGIN`, levels = order, exclude = NULL)
levels(env$`PLACE OF ORIGIN`)

place <- env %>% select(CULTURE_originalName, `PLACE OF ORIGIN`) %>% 
  gather(`PLACE OF ORIGIN`, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "PLACE OF ORIGIN", "Isolation location")) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Isolation location") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())


#cruise
colnames(env)
env %>% distinct(CRUISE)
#env <- env %>% mutate(CRUISE = ifelse(is.na(CRUISE), "NA", CRUISE))
#env %>% distinct(CRUISE)
env %>% filter(is.na(CRUISE)) %>% nrow()
env <- env %>% mutate(CRUISE = ifelse(str_detect(CRUISE, "EqPac/IRONEX\\?"), "Possibly EqPac/IRONEX", CRUISE))
env %>% filter(is.na(CRUISE)) %>% nrow()
env %>% distinct(CRUISE)

order <- env %>% arrange(CRUISE) %>% distinct(CRUISE)
order <- order %>% select(CRUISE)
order <- order %>% as.list()
order <- order$CRUISE
order

#env <- env %>% mutate(CRUISE = ifelse(is.na(CRUISE), "NA", CRUISE))
env$CRUISE <- factor(env$CRUISE, levels = order, exclude = NULL)
levels(env$CRUISE)

cruise <- env %>% select(CULTURE_originalName, CRUISE) %>% 
  gather(CRUISE, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "CRUISE", "Cruise")) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Cruise") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())



###deal with depth
colnames(env)
str(env$DEPTH)
env %>% distinct(DEPTH) %>% arrange(DEPTH)

env <- env %>% mutate(depthRange = ifelse(DEPTH <= 50, "0-50m", NA))
env <- env %>% mutate(depthRange = ifelse(DEPTH > 50 & DEPTH <= 100, "50-100m", depthRange))
env <- env %>% mutate(depthRange = ifelse(DEPTH > 100 & DEPTH <= 150, "100-150m", depthRange))
env <- env %>% mutate(depthRange = ifelse(DEPTH > 150 & DEPTH <= 200, "150-200m", depthRange))

env %>% filter(is.na(depthRange)) %>% distinct(DEPTH)
env %>% filter(is.na(depthRange)) %>% select(DEPTH, depthRange)

#env <- env %>% mutate(depthRange = ifelse(is.na(depthRange), "NA", depthRange))

env %>% distinct(DEPTH, depthRange) %>% arrange(DEPTH)
env %>% distinct(depthRange)

#is this okay to do???
env$depthRange <- factor(env$depthRange, levels = c("0-50m", "50-100m", "100-150m", "150-200m", NA), exclude = NULL)
levels(env$depthRange)

depth <- env %>% select(CULTURE_originalName, depthRange) %>% 
  gather(depthRange, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "depthRange", "Depth")) %>% 
  mutate(value = factor(value, levels = c("0-50m", "50-100m", "100-150m", "150-200m", NA), exclude = NULL)) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Depth") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())




###isolator
colnames(env)
env %>% distinct(ISOLATOR)
#env <- env %>% mutate(ISOLATOR = ifelse(is.na(ISOLATOR), "NA", ISOLATOR))
#env %>% distinct(ISOLATOR)

order <- env %>% arrange(ISOLATOR) %>% distinct(ISOLATOR)
order <- order %>% select(ISOLATOR)
order <- order %>% as.list()
order <- order$ISOLATOR
order

#env <- env %>% mutate(ISOLATOR = ifelse(is.na(ISOLATOR), "NA", ISOLATOR))
env %>% distinct(ISOLATOR)
env$ISOLATOR <- factor(env$ISOLATOR, levels = order, exclude = NULL)
levels(env$ISOLATOR)

isolator <- env %>% select(CULTURE_originalName, ISOLATOR) %>% 
  gather(ISOLATOR, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "ISOLATOR", "Isolator")) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Isolator") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())


###method of isolation
colnames(env)
env %>% distinct(`METHOD (for isolation)`)
#env <- env %>% mutate(`METHOD (for isolation)` = ifelse(is.na(`METHOD (for isolation)`), "NA", `METHOD (for isolation)`))
#env %>% distinct(`METHOD (for isolation)`)

order <- env %>% arrange(`METHOD (for isolation)`) %>% distinct(`METHOD (for isolation)`)
order <- order %>% select(`METHOD (for isolation)`)
order <- order %>% as.list()
order <- order$`METHOD (for isolation)`
order

#env <- env %>% mutate(`METHOD (for isolation)` = ifelse(is.na(`METHOD (for isolation)`), "NA", `METHOD (for isolation)`))
env %>% distinct(`METHOD (for isolation)`)
env$`METHOD (for isolation)` <- factor(env$`METHOD (for isolation)`, levels = order, exclude = NULL)
levels(env$`METHOD (for isolation)`)

method <- env %>% select(CULTURE_originalName, `METHOD (for isolation)`) %>% 
  gather(`METHOD (for isolation)`, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "METHOD \\(for isolation\\)", "Isolation method")) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Isolation method") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())



###year
env <- env %>% mutate(year = str_extract(`DATE ISOLATED_Edited`, "[0-9]*"))
env %>% select(`DATE ISOLATED_Edited`, year)
env %>% filter(is.na(year)) %>% distinct(`DATE ISOLATED_Edited`)
env %>% distinct(year) %>% arrange(year)
#env <- env %>% mutate(year = ifelse(is.na(year), "NA", year))
#env %>% distinct(year) %>% arrange(year)

order <- env %>% arrange(year) %>% distinct(year)
order <- order %>% select(year)
order <- order %>% as.list()
order <- order$year
order
env$year <- factor(env$year, levels = order, exclude = NULL)
levels(env$year)

year <- env %>% select(CULTURE_originalName, year) %>% 
  gather(year, key = "var", value = "value") %>% 
  mutate(var = str_replace(var, "year", "Year")) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Year") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())




###month
colnames(env)
env %>% distinct(monthNew)
#env <- env %>% mutate(`METHOD (for isolation)` = ifelse(is.na(`METHOD (for isolation)`), "NA", `METHOD (for isolation)`))
#env %>% distinct(`METHOD (for isolation)`)

order <- env %>% arrange(month) %>% distinct(monthNew)
order <- order %>% select(monthNew)
order <- order %>% as.list()
order <- order$monthNew
order

#env <- env %>% mutate(`METHOD (for isolation)` = ifelse(is.na(`METHOD (for isolation)`), "NA", `METHOD (for isolation)`))
#env %>% distinct(monthNew)
env$monthNew <- factor(env$monthNew, levels = order, exclude = NULL)
levels(env$monthNew)

month <- env %>% select(CULTURE_originalName, monthNew) %>% 
  gather(monthNew, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "monthNew", "Month")) %>% 
  mutate(value = factor(value, levels = order, exclude = NULL)) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Month") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())




###ocean type
colnames(env)
env %>% distinct(oceanType)
#env <- env %>% mutate(`METHOD (for isolation)` = ifelse(is.na(`METHOD (for isolation)`), "NA", `METHOD (for isolation)`))
#env %>% distinct(`METHOD (for isolation)`)

order <- env %>% arrange(oceanType) %>% distinct(oceanType)
order <- order %>% select(oceanType)
order <- order %>% as.list()
order <- order$oceanType
order

#env <- env %>% mutate(`METHOD (for isolation)` = ifelse(is.na(`METHOD (for isolation)`), "NA", `METHOD (for isolation)`))
#env %>% distinct(monthNew)
env$oceanType <- factor(env$oceanType, levels = order, exclude = NULL)
levels(env$oceanType)

oceanType <- env %>% select(CULTURE_originalName, oceanType) %>% 
  gather(oceanType, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "oceanType", "Isolation location type")) %>% 
  mutate(value = factor(value, levels = order, exclude = NULL)) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Isolation location type") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())

oceanType



###month, year

env %>% select(monthNew, year) %>% head()
env %>% select(monthNew, year) %>% str()
as.character(env$monthNew)
as.character(env$year)

#the culture that I am not sure about the month, day doesn't have a month value which is good
env %>% filter(is.na(month))

#makes a variable that has both the month and year of isolation
env <- env %>% mutate(monthNewYear = str_c(as.character(env$monthNew), as.character(env$year), sep = " "))

env %>% filter(is.na(monthNewYear)) %>% nrow()

#env %>% filter(is.na(monthNewYear)) %>% View()

env %>% distinct(monthNewYear) %>% nrow()

env %>% arrange(year, month) %>% distinct(monthNewYear)
env %>% arrange(year, month) %>% distinct(monthNewYear) %>% nrow()

order <- env %>% arrange(year, month) %>% distinct(monthNewYear)
order
order <- order %>% select(monthNewYear)
order
order <- order %>% as.list()
order <- order$monthNewYear
order

#env <- env %>% mutate(`METHOD (for isolation)` = ifelse(is.na(`METHOD (for isolation)`), "NA", `METHOD (for isolation)`))
#env %>% distinct(monthNew)

#puts the monthNewYear variable values in chronological order
env$monthNewYear <- factor(env$monthNewYear, levels = order, exclude = NULL)
levels(env$monthNewYear)
levels(env$monthNewYear) %>% length()

monthNewYear <- env %>% select(CULTURE_originalName, monthNewYear) %>% 
  gather(monthNewYear, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "monthNewYear", "Month, year")) %>% 
  mutate(value = factor(value, levels = order, exclude = NULL)) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_discrete(na.value="white") + labs(x = "", y = "", fill = "Month, year") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())

monthNewYear


###age of cultures

str(env$`DATE ISOLATED_Edited`)
str(ymd(env$`DATE ISOLATED_Edited`))

env$`DATE ISOLATED_Edited`
ymd(env$`DATE ISOLATED_Edited`)
ymd("2019-06-01 UTC")

#makes a variable for the age of the culture at the time they were sequenced
env <- env %>% mutate(yearsInAge = interval(ymd(env$`DATE ISOLATED_Edited`), ymd("2019-06-01 UTC"))/years(1))

env %>% filter(CULTURE_originalName == "WH6501")
env %>% filter(is.na(`DATE ISOLATED_Edited`)) %>% nrow()
env %>% filter(is.na(yearsInAge)) %>% nrow()
env %>% filter(is.na(as.character(monthNewYear))) %>% nrow()

str(env$yearsInAge)

yearsInAge <- env %>% select(CULTURE_originalName, yearsInAge) %>% 
  gather(yearsInAge, key = "var", value = "value") %>%
  mutate(var = str_replace(var, "yearsInAge", "Age of culture (years)")) %>% 
  ggplot() + geom_tile(aes(y = CULTURE_originalName, x = var, fill = value), width=0.7, height=0.7) + 
  theme(legend.position = 'none') + theme_classic() + 
  scale_fill_continuous(na.value="white", low = "#FFF7EC", high = "#7F0000") + 
  labs(x = "", y = "", fill = "Age of culture (years)") + 
  theme(axis.line=element_blank(),axis.ticks=element_blank())

yearsInAge


colnames(env)

htm <- grid.arrange(cultureType + theme(legend.position = 'none'), 
             ecotype + theme(legend.position = 'none'), 
             clade + theme(legend.position = 'none'), 
             cruise + theme(legend.position = 'none'), 
             place + theme(legend.position = 'none'),
             oceanType + theme(legend.position = 'none'),
             depth + theme(legend.position = 'none'), 
             year + theme(legend.position = 'none'), 
             month + theme(legend.position = 'none'), 
             isolator + theme(legend.position = 'none'), 
             method + theme(legend.position = 'none'), 
             monthNewYear + theme(legend.position = 'none'), 
             yearsInAge + theme(legend.position = 'none'), ncol=13)

ggsave("metadataHeatmap.svg", htm, height = 10, width = 30, dpi = 600)
ggsave("metadataHeatmap.png", htm, height = 30, width = 30, dpi = 600)


leg <- grid.arrange(cultureType %>% get_legend(), 
             ecotype %>% get_legend(), 
             clade %>% get_legend(), 
             cruise %>% get_legend(), 
             place %>% get_legend(), 
             oceanType %>% get_legend(), 
             depth %>% get_legend(), 
             year %>% get_legend(), 
             month %>% get_legend(),
             isolator %>% get_legend(), 
             method %>% get_legend(), 
             monthNewYear %>% get_legend(), 
             yearsInAge %>% get_legend(), ncol=13)

ggsave("metadataLegend.svg", leg, height = 14, width = 30, dpi = 600)

rev(sampleTree$labels[sampleTree$order])
write.table(rev(sampleTree$labels[sampleTree$order]), "unifracTreeClusteringForMetadataHeatmapOrder.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)




###run statistical test of pca 
###to see clustering is significantly correlated with any of the metadata

#sean says that * means interaction in formula while + is not for interactions 
#between variables

#https://groups.google.com/forum/#!topic/qiime-forum/7O8Bpz6T0VU
#https://rpubs.com/collnell/manova
#https://ichthyology.usm.edu/courses/multivariate/feb_7.pdf
#https://gist.github.com/claczny/3415270a6c6919969bff79d2c246a527

#unifracRes is a distance object
str(unifracRes)
unifracRes

env$CULTURE_originalName

levels(env$CULTURE_originalName)
sampleTree$labels[sampleTree$order]

labels(unifracRes)
length(labels(unifracRes))

##puts the sample IDs in env in the same order as they appear in unifracRes

order <- labels(unifracRes) %>% as.list()
order
env$CULTURE_originalName <- factor(env$CULTURE_originalName, levels = order)

levels(env$CULTURE_originalName) %>% length()
labels(unifracRes) %>% length()

levels(env$CULTURE_originalName)
labels(unifracRes)

unique(levels(env$CULTURE_originalName) == labels(unifracRes))

env <- env %>% arrange(CULTURE_originalName)

env$CULTURE_originalName
labels(unifracRes)
unique(env$CULTURE_originalName == labels(unifracRes))

colnames(env)

nrow(env)
labels(unifracRes) %>% length()

set.seed(7)

#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
adonis2(unifracRes ~ cultureType, data = env, method="unifrac")
#output tells us that our adonis test is significant so we can reject the null 
#hypothesis that our culture types have the same centroid

set.seed(7)

beta <- betadisper(unifracRes, env$cultureType)
permutest(beta)
#betadisper results are significant, meaning we can reject the null hypothesis that 
#our groups have the same dispersions, which means we are not confident that our 
#adonis result is a real result, and not due to differences in group dispersions

colnames(env)

#cultureType, clade, cruise, place of origin, depthRange, month are the variables that sean
#wants me to include in the formula
#env %>% select(cultureType, CLADE, CRUISE, `PLACE OF ORIGIN`, depthRange, monthNew) %>% View()

set.seed(7)

adonis2(unifracRes~cultureType+CLADE+CRUISE+`PLACE OF ORIGIN`+depthRange+monthNew, data=env, permutations = 999, method="unifrac")


##clade is nested within cultureType so I want to make the grouping variable cultureType

set.seed(7)

adonis2(unifracRes~CLADE+CRUISE+`PLACE OF ORIGIN`+depthRange+monthNew, data=env, permutations = 999, method="unifrac", strata="cultureType")



##cruise is nested within place of origin so I tried one at a time but the 
##the results don't differ much so I am going to include them both

set.seed(7)

adonis2(unifracRes~CLADE+`PLACE OF ORIGIN`+depthRange+monthNew, data=env, permutations = 999, method="unifrac", strata="cultureType")

set.seed(7)

adonis2(unifracRes~CLADE+CRUISE+depthRange+monthNew, data=env, permutations = 999, method="unifrac", strata="cultureType")

colnames(env)
str(env)
str(env$yearsInAge)

##makes yearsInAge variable a factor so I can run adonis2 formula with it
env %>% distinct(yearsInAge) %>% nrow()
env %>% distinct(yearsInAge) %>% filter(is.na(yearsInAge))

ageOrder <- env %>% distinct(yearsInAge)
nrow(ageOrder)
ageOrder %>% filter(is.na(yearsInAge))
ageOrder
ageOrder <- ageOrder %>% as.list()
ageOrder
ageOrder <- ageOrder$yearsInAge
ageOrder
length(ageOrder)

str(env$yearsInAge)
env$yearsInAge <- factor(env$yearsInAge, levels = ageOrder, exclude = NULL)
str(env$yearsInAge)
levels(env$yearsInAge)
length(levels(env$yearsInAge))

colnames(env)

set.seed(7)

div <- adonis2(unifracRes~CLADE+CRUISE+`PLACE OF ORIGIN`+depthRange+monthNew+`METHOD (for isolation)`+yearsInAge, data=env, permutations = 999, method="unifrac", strata="cultureType")
div

write.table(div, "permutationTestUnifrac.txt", quote = FALSE)


set.seed(7)
adonis2(unifracRes~CLADE, data=env, permutations = 999, method="unifrac", strata="cultureType")

set.seed(7)
adonis2(unifracRes~CRUISE, data=env, permutations = 999, method="unifrac", strata="cultureType")

set.seed(7)
adonis2(unifracRes~`PLACE OF ORIGIN`, data=env, permutations = 999, method="unifrac", strata="cultureType")

set.seed(7)
adonis2(unifracRes~depthRange, data=env, permutations = 999, method="unifrac", strata="cultureType")

set.seed(7)
adonis2(unifracRes~monthNew, data=env, permutations = 999, method="unifrac", strata="cultureType")

set.seed(7)
adonis2(unifracRes~`METHOD (for isolation)`, data=env, permutations = 999, method="unifrac", strata="cultureType")

set.seed(7)
adonis2(unifracRes~yearsInAge, data=env, permutations = 999, method="unifrac", strata="cultureType")











