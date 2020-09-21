library(tidyverse)
library(phyloseq)
library(ggplot2)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads taxonomy of sean and diatom combined
#from analysisSeanAndDiatom_V2.txt (checked)
tax <- read_tsv("diatomNoMitochondriaChloroplast_seanNoContamination_noMock_withProSyn_no1223Taxonomy.tsv")

#gets rid of first row because it just has variable names
tax <- tax[-1,]

#correct number of clusters
nrow(tax)

#there are no eukaryotic, mitochondria, or chloroplast sequences which is good
tax %>% filter(str_detect(Taxon, "eukary|Eukary|chloroplast|Chloroplast|mitochondria|Mitochondria"))

nrow(tax)

#makes a variable for the class
tax <- tax %>% mutate(class = str_extract(Taxon, "c__[a-z,A-Z,0-9,\\[,\\],\\-]*;{0,1}"))
tax %>% distinct(class)
tax <- tax %>% mutate(class = str_replace_all(class, ";", ""))
tax %>% distinct(class)
tax <- tax %>% mutate(class = str_replace_all(class, " ", ""))
tax %>% distinct(class)
tax %>% filter(is.na(class))
tax$class <- ifelse(is.na(tax$class), "Unclassified at class level", tax$class)
tax %>% filter(class == "c__")
tax %>% distinct(class)

#gets rid of "[" and "]" in class names
tax$class <- str_replace(tax$class, "\\[", "")
tax$class <- str_replace(tax$class, "\\]", "")

#gets rid of "c__" in class names
tax$class <- str_replace(tax$class, "c__", "")

tax %>% distinct(class)

#checked
#View(tax %>% distinct(Taxon, class))

#correct number of clusters
nrow(tax)


###diatom

#loads the abundance of sequences in each sample for the diatom dataset
#Sean sent this to me over Slack on 11/26/2019
seq <- read.table("Diatoms16S.txt")

numCols <- ncol(seq)

#makes the column IDs match what I used for qiime2
cols <- str_c(1:numCols, "diatom")
colnames(seq) <- cols

nrow(seq)

#gets list of samples
samples <- rownames(seq)

#makes a variable for the sample in seq instead of the samples just being
#the rownames
seq <- seq %>% mutate(sample = samples)

#in the diatom dataset, there are 404 sequences (without clustering) and without 
#excluding mitochondria, chloroplast, or low relative abundance sequences
numCols
ncol(seq)
nrow(seq)

colnames(seq)[400:405]

#gathers the diatom abundance data into long form
seqG <- seq %>% gather(1:404, key = "sequence", value = "abundance")

str(seqG)

#in the diatom dataset, gets rid of the rows for which a sequence is 
#not present in a sample
seqG <- seqG %>% filter(abundance > 0)

#there are 404 diatom sequences
seqG %>% distinct(sequence) %>% nrow()

str(seqG)

#calculates the total number of sequence in each sample
summ <- seqG %>% group_by(sample) %>% summarize(total = sum(abundance))

nrow(seqG)

#adds total number of reads in each sample to seqG
seqG <- seqG %>% left_join(summ, by = c("sample"))

nrow(seqG)

seqG %>% filter(is.na(total))

str(seqG)

#calculates the relative abundance (prop) of each sequence in each sample
seqG <- seqG %>% mutate(prop = abundance/total)

seqG %>% group_by(sample) %>% summarize(totalProp = sum(prop)) %>% filter(totalProp != 1)

#excludes the sequences with low relative abundance 
seqG <- seqG %>% filter(prop >= 0.002)

seqG %>% distinct(sequence) %>% nrow()

#gets rid of unnecessary variables
seqG <- seqG %>% select(-c(total, prop))

#calculates the total number of sequence in each sample after excluding the low abundance sequences
summ <- seqG %>% group_by(sample) %>% summarize(total = sum(abundance))

nrow(seqG)

#adds total number of reads in each sample after excluding the low abundance sequences to seqG
seqG <- seqG %>% left_join(summ, by = c("sample"))

nrow(seqG)

seqG %>% filter(is.na(total))

str(seqG)

#calculates the relative abundance (prop) of each sequence in each sample after excluding the low abundance sequences
seqG <- seqG %>% mutate(prop = abundance/total)

seqG %>% group_by(sample) %>% summarize(totalProp = sum(prop)) %>% filter(totalProp != 1)

write_csv(seqG, "seqG_diatom_WithMitochondriaChloroplast_noLowProp.csv")



###sean

#for Sean dataset, loads abundance of sequences in each sample
#this came from taxonomyStackedBarPlot_noContamination_noProSyn_withoutSample1223.R (checked)
#Pro and Syn sequences have not been excluded
#samples 1223 and JW3 have not been excluded
#prop < .002 have been excluded 
#Mock, Con, Epi, Pcn, Lin samples have been excluded
seanSeq <- read_csv("seqG_forDiatomAndSeanCombinedAnalysis.csv")

#excludes this sample because it only has one sequence
seanSeq %>% filter(sample == "JW3")
seanSeq <- seanSeq %>% filter(sample != "JW3")

#excludes 1223 because it doesn't have Pro or Syn
seanSeq %>% filter(sample == "1223")
seanSeq <- seanSeq %>% filter(sample != "1223")

#correct number of samples
seanSeq %>% distinct(sample) %>% nrow()

#this includes Pro and Syn sequences
seanSeq %>% distinct(sequence) %>% nrow()

#low prop sequences in samples have already been filtered out
seanSeq %>% filter(propSample < .002)

#gets rid of unnecessary variables
seanSeq <- seanSeq %>% select(-c(totalSeqs, propSample))

#calculates the total number of sequences in each sample after excluding
#low prop sequences
summ <- seanSeq %>% group_by(sample) %>% summarize(totalWithoutLowProp = sum(abundance))

nrow(summ)

nrow(seanSeq)

#adds total number of sequences in each sample after excluding
#low prop sequences to seanSeq
seanSeq <- seanSeq %>% left_join(summ, by = c("sample"))

nrow(seanSeq)

seanSeq %>% filter(is.na(totalWithoutLowProp))

#calculates the relative abundance of each sequence, after excluding low prop sequences
seanSeq <- seanSeq %>% mutate(propWithoutLowProp = abundance/totalWithoutLowProp)

str(seanSeq)

#makes the sequence IDs match what I used for qiime2
seanSeq <- seanSeq %>% mutate(sequence = str_replace(sequence, "col", ""))
seanSeq <- seanSeq %>% mutate(sequence = str_c(sequence, "sean"))

#renames variables
colnames(seanSeq)[4:5] <- c("total", "prop")

seanSeq %>% group_by(sample) %>% summarize(totalProp = sum(prop)) %>% filter(totalProp != 1)


###combined diatom and sean sequence abundance data

#makes a variable for the experiment/origin of the data
seqG <- seqG %>% mutate(sampleSource = "diatom")
seanSeq <- seanSeq %>% mutate(sampleSource = "sean")

#combines the diatom and sean data
seqComb <- bind_rows(seqG, seanSeq)

seqComb %>% distinct(sequence) %>% nrow()

seqComb %>% group_by(sampleSource) %>% distinct(sequence) %>% summarize(n = n())

###fix sample IDs

#1213 is 1313
#WH7803 is WH7803
#SYN1320 is SYN1220
#9201.2 is 9311
#9201.1 is 9201

#AS9601 is correct ID rather than ASN9601
#SYN9503 is correct ID rather than S9503

seqComb$sample <- ifelse(seqComb$sample == "1213", "1313", seqComb$sample)
#seqComb$sample <- ifelse(seqComb$sample == "WH7803", "WH7801", seqComb$sample)
seqComb$sample <- ifelse(seqComb$sample == "SYN1320", "SYN1220", seqComb$sample)
seqComb$sample <- ifelse(seqComb$sample == "9201.2", "9311", seqComb$sample)
seqComb$sample <- ifelse(seqComb$sample == "9201.1", "9201", seqComb$sample)

seqComb$sample <- ifelse(seqComb$sample == "ASN9601", "AS9601", seqComb$sample)
seqComb$sample <- ifelse(seqComb$sample == "S9503", "SYN9503", seqComb$sample)

seqComb$sample <- ifelse(seqComb$sample == "8102", "SYNCLADEXVI", seqComb$sample)


seqComb %>% distinct(sample) %>% nrow()

#loads IDs for diatom samples
#from extractDiatomMetadata.R
ID <- read_csv("diatomIDs.csv")

#adds diatom sample IDs to seqComb
nrow(seqComb)
seqComb <- seqComb %>% left_join(ID, by = c("sample" = "Run"))
nrow(seqComb)

seqComb %>% filter(sampleSource == "diatom") %>% filter(is.na(ID))

#replaces diatom sample IDs with the newly added ones
seqComb <- seqComb %>% mutate(sample = ifelse(sampleSource == "diatom", ID, sample))

nrow(seqComb)

seqComb <- seqComb %>% select(-ID)

seqComb %>% distinct(sample) %>% nrow()
seqComb %>% distinct(sample)

###sean and diatom combined biom

#checked in checkFeatureTableFromBiomSeanDiatomDN1Edited.R
#loads biom data which came from running the combined diatom and sean data through qiime2
#in analysisSeanAndDiatom_V2.txt
biom <- read.table("feature-table-fromBiom-diatomNoMitochondriaChloroplast_seanNoContamination_noMock_withProSyn_no1223_dn-1Edited.tsv", colClasses = "character")

str(biom)

#gets the input sequence/sample IDs
#they are the first row of biom
namesNew <- biom[1,]
namesNew

ncol(biom)
length(namesNew)

#makes the column names the sequence/sample IDs
colnames(biom) <- namesNew

#gets rid of the first row because it just has variable names
biom <- biom[-1,]

#there are 133 clustered features which is correct
nrow(biom)

#there are 327 sequence/samples, which is the correct number of sequences
ncol(biom)

colnames(biom)[1:5]
colnames(biom)[322:328]

#gathers the biom dataset into long format
biomG <- biom %>% gather(2:328, key = "sample", value = "present")

#there are 133 clustered features and 327 sequences/samples
biomG %>% distinct(OTU_ID) %>% nrow()
biomG %>% distinct(sample) %>% nrow()

str(biomG)

#makes variable indicating whether a feature is present in a sample numeric
biomG$present <- as.numeric(biomG$present)

biomG %>% filter(present > 0) %>% distinct(OTU_ID) %>% nrow()
biomG %>% filter(present > 0) %>% distinct(sample) %>% nrow()

biomG %>% filter(present > 1)

#gets just the rows for which a feature is present in a sample
biomG <- biomG %>% filter(present > 0)

biomG %>% group_by(OTU_ID) %>% summarize(n = n()) %>% filter(n > 1)
biomG %>% group_by(OTU_ID) %>% distinct(sample) %>% summarize(n = n()) %>% filter(n > 1)


#seqComb has more sequences than biomG because I excluded 
#the mitochondria and chloroplast diatom sequences from biomG 
#but not from seqComb
biomG %>% distinct(sample) %>% nrow()
seqComb %>% distinct(sequence) %>% nrow()

biomG %>% anti_join(seqComb, by = c("sample" = "sequence")) 
seqComb %>% anti_join(biomG, by = c("sequence" = "sample")) %>% distinct(sequence)

#excludes the sequences in seqComb that were identified as mitochondria or chloroplast
seqComb <- seqComb %>% semi_join(biomG, by = c("sequence" = "sample"))

#now there are the same number of sequences in biomG and seqComb
biomG %>% distinct(sample) %>% nrow()
seqComb %>% distinct(sequence) %>% nrow()

biomG %>% anti_join(seqComb, by = c("sample" = "sequence")) 
seqComb %>% anti_join(biomG, by = c("sequence" = "sample")) %>% distinct(sequence)

nrow(biomG)

#adds the abundance of each sequence in each sample data to biomG
biomG <- biomG %>% left_join(seqComb, by = c("sample" = "sequence"))

nrow(biomG)

biomG %>% filter(is.na(sample.y))

#there are 327 sequences
biomG %>% distinct(sample) %>% nrow()

#there are 93 samples
biomG %>% distinct(sample.y) %>% nrow()

#theere are 133 clusters
biomG %>% distinct(OTU_ID) %>% nrow()
tax %>% distinct(`Feature ID`) %>% nrow()

nrow(biomG)

#adds taxonomy to biomG
biomG <- biomG %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))

nrow(biomG)

biomG %>% filter(is.na(class))


###calculate prop abundance of each diatom sequence in each sample

#gets the total number of reads in each diatom sample after excluding low relative abundance sequences and 
#sequences that were identified as mitochondria or chloroplast
summ <- biomG %>% filter(sampleSource == "diatom") %>% group_by(sample.y) %>% summarize(total = sum(abundance))

nrow(biomG)

#adds total number of reads in each diatom sample after excluding low relative abundance sequences and 
#sequences that were identified as mitochondria or chloroplast to seqG
biomG <- biomG %>% left_join(summ, by = c("sample.y"))

nrow(biomG)

biomG %>% filter(sampleSource == "diatom") %>% filter(is.na(total.y))

#makes one variable for the total number of sequences in a sample
biomG <- biomG %>% mutate(total.x = ifelse(sampleSource == "diatom", total.y, total.x))

biomG %>% filter(is.na(total.x))

biomG %>% filter(abundance == 0)

#recalculates in the relative abundance for each diatom sequence in each sample
biomG <- biomG %>% mutate(prop = ifelse(sampleSource == "diatom", abundance/total.x, prop))

biomG %>% filter(is.na(prop))

#gets rid of unnecessary variable
biomG <- biomG %>% select(-total.y)

colnames(biomG)[6] <- "total"


biomG %>% group_by(sample.y) %>% summarize(totalProp = sum(prop)) %>% filter(totalProp != 1)

#these are the clusters arranged by their number of sequences
biomG %>% group_by(OTU_ID) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(desc(n))

#these are the clusters and their clusters that are found in both the sean and diatom samples
biomG %>% group_by(OTU_ID, class) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n > 1)

#there are 133 clusters
biomG %>% distinct(OTU_ID) %>% summarize(n = n())
#12 clusters are found in both the sean and diatom sequences
biomG %>% group_by(OTU_ID) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n > 1) %>% ungroup() %>% distinct(OTU_ID) %>% nrow()
#these are the distinct taxonomies of the clusters that were found in both the sean and diatom sequences
biomG %>% group_by(OTU_ID, class, Taxon) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n > 1) %>% ungroup() %>% distinct(Taxon) %>% arrange(Taxon)

#gets the clusters that were found in just one of the sample sources: either sean or diatom but not both
foundInOneSource <- biomG %>% group_by(OTU_ID) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n == 1) %>% distinct(OTU_ID)

nrow(foundInOneSource)

#104 clusters were just found in sean samples
biomG %>% semi_join(foundInOneSource, by = c("OTU_ID")) %>% filter(sampleSource == "sean") %>% distinct(OTU_ID) %>% summarize(n = n())

#17 clusters were just found in diatom samples
biomG %>% semi_join(foundInOneSource, by = c("OTU_ID")) %>% filter(sampleSource == "diatom") %>% distinct(OTU_ID) %>% summarize(n = n())

12+104+17
104+17

#Synechococcophycideae clusters were only found in sean samples not diatom samples
biomG %>% filter(class == "Synechococcophycideae") %>% group_by(OTU_ID, class) %>% distinct(sampleSource)

#these are the Synechococcophycideae clusters arranged by the number of sequences they have
biomG %>% filter(class == "Synechococcophycideae") %>% group_by(OTU_ID) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(desc(n))

#these are the number of sequences in each sample source
biomG %>% group_by(sampleSource) %>% distinct(sample) %>% summarize(n = n())

#makes a variable for whether sample is a Pro or Syn culture
biomG <- biomG %>% mutate(culture = ifelse(str_detect(sample.y, "SYN|WH|9220|S9503|8102|RS9916") & sampleSource == "sean", "Synechococcus", ifelse(sampleSource == "sean", "Prochlorococcus", NA)))

nrow(biomG)

#for sean seqs, if the class is Synechococcophycideae, adds whether the sequence is Pro or Syn
biomG <- biomG %>% mutate(class = ifelse(class == "Synechococcophycideae" & sampleSource == "sean", str_c(class, culture, sep = " "), class))

nrow(biomG)

#there are no Synechococcophycideae clusters in the diatom samples
biomG %>% filter(sampleSource == "diatom") %>% filter(class == "Synechococcophycideae") %>% distinct(Taxon)
biomG %>% filter(sampleSource == "diatom") %>% filter(class == "Synechococcophycideae") %>% distinct(sample)

biomG %>% filter(sampleSource == "diatom" & class == "Synechococcophycideae") %>% distinct(sample) %>% 
  mutate(sample = str_replace(sample, "diatom", "")) %>% 
  mutate(sample = str_c(sample, "_ending")) %>% 
  write.table("diatomSamples_SynechococcophycideaeSeqs.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

#for diatom seqs, if the class is Synechococcophycideae, adds whether the sequence is Pro or Syn
biomG <- biomG %>% mutate(class = ifelse(class == "Synechococcophycideae" & sampleSource == "diatom", str_c(class, "Synechococcus", sep = " "), class))


#loads rooted tree
#this came from analysisSeanAndDiatom_V2.txt (checked)
tree <- read_tree("~/Dropbox (MIT)/Sean16SSep2019/exported-rooted-rep-seqs-diatomNoMitochondriaChloroplast_seanNoContamination_noMock_withProSyn_no1223-dn-1_tree/tree.nwk")

##to make otu phyloseq object

head(biomG)

#gets necessary variables to spread biomG
seqSpread <- biomG %>% select(OTU_ID, sample.y, abundance)

head(seqSpread)

#these are the clusters that have multiple sequences 
seqSpread %>% group_by(OTU_ID, sample.y) %>% summarize(n = n()) %>% filter(n > 1)

#gets the total abundance of each cluster in each sample
seqSpreadSumm <- seqSpread %>% group_by(OTU_ID, sample.y) %>% summarize(sumAbundance = sum(abundance))

seqSpreadSumm %>% group_by(OTU_ID, sample.y) %>% summarize(n = n()) %>% filter(n > 1)

head(seqSpreadSumm)

#makes total abundance of each cluster in each sample into a dataframe
seqSpreadSumm <- seqSpreadSumm %>% ungroup() %>% as.data.frame()

#puts into wide form
seqSpreadSumm <- seqSpreadSumm %>% spread(key = OTU_ID, value = sumAbundance, fill = 0)

#makes the rownames the sample IDs
rownames(seqSpreadSumm) <- seqSpreadSumm$sample.y
seqSpreadSumm <- seqSpreadSumm %>% select(-c(sample.y))

dim(seqSpreadSumm)

#makes phyloseq object out of otu table
otu <- otu_table(seqSpreadSumm, taxa_are_rows = FALSE)


otuDF <- colnames(otu) %>% as.data.frame()
treeDF <- tree$tip.label %>% as.data.frame()

otuDF %>% nrow()
treeDF %>% nrow()

colnames(otuDF) <- "V1"
colnames(treeDF) <- "V1"

otuDF %>% anti_join(treeDF, by = c("V1"))
treeDF %>% anti_join(otuDF, by = c("V1"))

taxa_names(otu)
taxa_names(tree)

#makes combined phyloseq object out of otu, TAX, tree
#physeq <- phyloseq(otu, TAX, tree)
physeq <- phyloseq(otu, tree)

set.seed(7)

#there is a warning
unifracRes <- UniFrac(physeq, weighted = FALSE)

#makes tree
sampleTree <- unifracRes %>% hclust()

#plot_heatmap(physeq, taxa.label="shortTaxa")

plot(sampleTree)

sampleTree$labels

sampleTree$order

sampleTree$labels[sampleTree$order]

#puts seqG sample variable in the same order as unifrac clustering of samples
biomG$sample.y <- factor(biomG$sample.y, levels = sampleTree$labels[sampleTree$order])

levels(biomG$sample.y)

#colors <- c("cyan3", "sienna4", "olivedrab3", "mediumvioletred", "gray58", "lightsteelblue1", "springgreen4", "hotpink1", "red", "tan3", "blueviolet", "plum", "goldenrod1", "black", "blue1", "lightcoral", "khaki2", "gray92", "springgreen", "pink", "darkorange1", "yellow")

colors <- c("cyan3", "mediumvioletred", "gray58", "springgreen4", "hotpink1", "red", "tan3", "blueviolet", "black", "blue1", "lightcoral", "springgreen", "pink", "darkorange1")

biomG %>% ggplot(aes(x = sample.y, y = prop, fill = class)) + geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y = "Relative abundance", fill = "") + scale_fill_manual(values = colors) 

ggsave("diatomAndSeanTaxonomyStackedBarPlot_withProSyn.png", dpi = 600, height = 6, width = 15)


png("diatomAndSeanUnifracTree_withProSyn.png")

plot(sampleTree)

dev.off()






###without Pro and Syn

biomG %>% distinct(class)

#makes biomG dataframe without Pro or Syn
noProSyn <- biomG %>% filter(class != "Synechococcophycideae Prochlorococcus") %>% filter(class != "Synechococcophycideae Synechococcus")

noProSyn %>% distinct(class)

#gets rid of unnecessary variables
noProSyn <- noProSyn %>% select(-c(total, prop))

#gets the total number of reads in each sample after excluding the Pro and Syn sequences
summ <- noProSyn %>% group_by(sample.y) %>% summarize(totalWithoutProSyn = sum(abundance))

nrow(summ)

nrow(noProSyn)

#adds total number of reads in each sample after excluding 
#Pro and Syn sequences to noProSyn
noProSyn <- noProSyn %>% left_join(summ, by = c("sample.y"))

nrow(noProSyn)

noProSyn %>% filter(is.na(totalWithoutProSyn))

#calculates the relative abundance of each sequence in each sample
#after excluding Pro and Syn sequences
noProSyn <- noProSyn %>% mutate(propWithoutProSyn = abundance/totalWithoutProSyn)

noProSyn %>% group_by(sample.y) %>% summarize(totalProp = sum(propWithoutProSyn)) %>% filter(totalProp != 1)

noProSyn %>% ggplot(aes(x = sample.y, y = propWithoutProSyn, fill = class)) + geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y = "Relative abundance", fill = "")


##to make otu phyloseq object

#gets just the necessary variables to spread the abundance data
seqSpread <- noProSyn %>% select(OTU_ID, sample.y, abundance)

seqSpread %>% group_by(OTU_ID, sample.y) %>% summarize(n = n()) %>% filter(n > 1)

#gets the total abundance of each cluster in each sample
seqSpreadSumm <- seqSpread %>% group_by(OTU_ID, sample.y) %>% summarize(sumAbundance = sum(abundance))

seqSpreadSumm %>% group_by(OTU_ID, sample.y) %>% summarize(n = n()) %>% filter(n > 1)

#makes seqSpreadSumm into a dataframe
seqSpreadSumm <- seqSpreadSumm %>% ungroup() %>% as.data.frame()

#puts seqSpreadSumm into wide form
seqSpreadSumm <- seqSpreadSumm %>% spread(key = OTU_ID, value = sumAbundance, fill = 0)

#makes the rownames the samples
rownames(seqSpreadSumm) <- seqSpreadSumm$sample.y
seqSpreadSumm <- seqSpreadSumm %>% select(-c(sample.y))

dim(seqSpreadSumm)

#makes phyloseq object out of otu table
otu <- otu_table(seqSpreadSumm, taxa_are_rows = FALSE)


###is it okay to do it this way--not make a new tree without Pro and Syn??

otuDF <- colnames(otu) %>% as.data.frame()
treeDF <- tree$tip.label %>% as.data.frame()

otuDF %>% nrow()
treeDF %>% nrow()

colnames(otuDF) <- "V1"
colnames(treeDF) <- "V1"

otuDF %>% anti_join(treeDF, by = c("V1"))
treeDF %>% anti_join(otuDF, by = c("V1"))

taxa_names(otu)
taxa_names(tree)

#makes combined phyloseq object out of otu, TAX, tree
#physeq <- phyloseq(otu, TAX, tree)
physeq <- phyloseq(otu, tree)

set.seed(7)

#there is a warning
unifracRes <- UniFrac(physeq, weighted = FALSE)

#makes tree
sampleTree <- unifracRes %>% hclust()

#plot_heatmap(physeq, taxa.label="shortTaxa")

plot(sampleTree)

sampleTree$labels

sampleTree$order


sampleTree$labels[sampleTree$order]

#puts seqG sample variable in the same order as unifrac clustering of samples
noProSyn$sample.y <- factor(noProSyn$sample.y, levels = sampleTree$labels[sampleTree$order])

levels(noProSyn$sample.y)

#colors <- c("cyan3", "sienna4", "olivedrab3", "mediumvioletred", "gray58", "lightsteelblue1", "springgreen4", "hotpink1", "red", "tan3", "blueviolet", "plum", "goldenrod1", "black", "blue1", "lightcoral", "khaki2", "gray92", "darkorange1", "yellow")

colors <- c("cyan3", "mediumvioletred", "gray58", "springgreen4", "hotpink1", "red", "tan3", "blueviolet", "black", "blue1", "lightcoral", "darkorange1")

noProSyn %>% ggplot(aes(x = sample.y, y = propWithoutProSyn, fill = class)) + geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y = "Relative abundance", fill = "") + scale_fill_manual(values = colors)

ggsave("diatomAndSeanTaxonomyStackedBarPlot_noProSyn.png", dpi = 600, height = 6, width = 15)


png("diatomAndSeanUnifracTree_noProSyn.png")

plot(sampleTree)

dev.off()

