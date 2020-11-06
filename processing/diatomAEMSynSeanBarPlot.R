library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggdendro)

setwd("~/Dropbox (MIT)/Sean16SSep2019/diatomAEMSynSean/")

#loads taxonomy
#diatomAEMSynSean-dn-1_taxonomy.tsv (checked) is from analysisDiatomAEMSynSean.txt (checked)
tax <- read_tsv("diatomAEMSynSean-dn-1_taxonomy.tsv")

#gets rid of first row because it just has variable names
tax <- tax[-1,]

#this is correct because there are 175 clustered features
nrow(tax)

#makes a variable for the class
tax <- tax %>% mutate(class = str_extract(Taxon, "c__[a-z,A-Z,0-9,\\[,\\],\\-]*;{0,1}"))
tax %>% distinct(class)
tax <- tax %>% mutate(class = str_replace_all(class, ";", ""))
tax %>% distinct(class)
tax <- tax %>% mutate(class = str_replace_all(class, " ", ""))
tax %>% distinct(class)
tax %>% filter(is.na(class)) %>% distinct(Taxon)
tax$class <- ifelse(is.na(tax$class), "Unclassified at class level", tax$class)
tax %>% filter(class == "c__")
tax %>% distinct(class)

#gets rid of "[" and "]" in class names
#tax$class <- str_replace(tax$class, "\\[", "")
#tax$class <- str_replace(tax$class, "\\]", "")

tax %>% distinct(class)

#gets rid of "c__" in class names
tax$class <- str_replace(tax$class, "c__", "")

tax %>% distinct(class)

#still has the correct number of rows
nrow(tax)

#checked
#View(tax %>% distinct(Taxon, class))


#there are no chloroplast or mitochondria sequences
tax %>% filter(str_detect(Taxon, "chloroplast|Chloroplast|mitochondria|Mitochondria"))


###syn 

#from AEMSynAndSeanBarPlot.R
syn <- read_csv("~/Dropbox (MIT)/Sean16SSep2019/comparisonToAEMPaperSynCultures/seqG_syn_noLowProp.csv")

###diatom 

#from seanDiatomTaxonomyBarPlot.R
diatom <- read_csv("~/Dropbox (MIT)/Sean16SSep2019/seqG_diatom_WithMitochondriaChloroplast_noLowProp.csv")

nrow(diatom)

#loads IDs for diatom samples
#from extractDiatomMetadata.R
ID <- read_csv("~/Dropbox (MIT)/Sean16SSep2019/diatomIDs.csv")

#adds diatom sample IDs to seqComb
nrow(diatom)
diatom <- diatom %>% left_join(ID, by = c("sample" = "Run"))
nrow(diatom)

diatom %>% filter(is.na(ID))

#replaces diatom sample IDs with the newly added ones
diatom <- diatom %>% mutate(sample = ID)

nrow(diatom)

diatom <- diatom %>% select(-ID)

diatom %>% distinct(sample) %>% nrow()
diatom %>% distinct(sample)
diatom %>% distinct(sequence) %>% nrow()

###sean 

#for Sean dataset, loads abundance of sequences in each sample
#this came from taxonomyStackedBarPlot_noContamination_noProSyn_withoutSample1223.R (checked)
#Pro and Syn sequences have not been excluded
#samples 1223 and JW3 have not been excluded
#prop < .002 have been excluded 
#Mock, Con, Epi, Pcn, Lin samples have been excluded
seanSeq <- read_csv("~/Dropbox (MIT)/Sean16SSep2019/seqG_forDiatomAndSeanCombinedAnalysis.csv")

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

seanSeq %>% filter(is.na(propWithoutLowProp))

str(seanSeq)

#makes the sequence IDs match what I used for qiime2
seanSeq <- seanSeq %>% mutate(sequence = str_replace(sequence, "col", ""))
seanSeq <- seanSeq %>% mutate(sequence = str_c(sequence, "sean"))

colnames(seanSeq)[4:5] <- c("total", "prop")




###combined diatom, AEM syn, and sean abundance data

syn %>% head()
diatom %>% head()
seanSeq %>% head()

syn %>% distinct(sequence) %>% nrow()
diatom %>% distinct(sequence) %>% nrow()
seanSeq %>% distinct(sequence) %>% nrow()

#makes a variable for the experiment/origin of the data
syn <- syn %>% mutate(sampleSource = "syn")
diatom <- diatom %>% mutate(sampleSource = "diatom")
seanSeq <- seanSeq %>% mutate(sampleSource = "sean")

syn %>% head()
diatom %>% head()
seanSeq %>% head()

#combines the diatom, AEM syn, and sean data
seqComb <- bind_rows(syn, diatom, seanSeq)

seqComb %>% distinct(sequence) %>% nrow()

head(seqComb)

###fix sample IDs

#1213 is 1313
#WH7803 is WH7801
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




###sean and syn combined biom

#checked in checkFeatureTableFromBiomDN1Edited.R
#loads biom data which came from running the combined diatom and sean data through qiime2
#in analysisDiatomAEMSynSean.txt
biom <- read.table("exported-feature-tableBiom-diatomAEMSynSean-dn-1Edited.tsv", colClasses = "character")

str(biom)

#gets the input sequence/sample IDs
#they are the first row of biom
namesNew <- biom[1,]

ncol(biom)
length(namesNew)

#makes the column names the sequence/sample IDs
colnames(biom) <- namesNew

#gets rid of the first row because it just has variable names
biom <- biom[-1,]

#there are 175 clustered features
nrow(biom)

#there are 426 sequence/samples; this is the correct number of sequences/samples
ncol(biom)

colnames(biom)[1:5]
colnames(biom)[422:427]

#gathers the biom dataset into long format
biomG <- biom %>% gather(2:427, key = "sample", value = "present")

#there are 175 clustered features and 426 sequences/samples
biomG %>% distinct(OTU_ID) %>% nrow()
biomG %>% distinct(sample) %>% nrow()

str(biomG)

#makes variable indicating whether a feature is present in a sample 
#numeric
biomG$present <- as.numeric(biomG$present)

str(biomG)

biomG %>% filter(present > 0) %>% distinct(OTU_ID) %>% nrow()
biomG %>% filter(present > 0) %>% distinct(sample) %>% nrow()

str(biomG)

biomG %>% filter(present > 1)

#gets just the rows for which a feature is present in a sample
biomG <- biomG %>% filter(present > 0)

#these are the clusters that have multiple sequences
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

nrow(biomG)

#adds the abundance of each sequence in each sample data to biomG
biomG <- biomG %>% left_join(seqComb, by = c("sample" = "sequence"))

nrow(biomG)

biomG %>% filter(is.na(sample.y))

#this is the correct number of sequences
biomG %>% distinct(sample) %>% nrow()

#this is the correct number of samples
biomG %>% distinct(sample.y) %>% nrow()

seanSeq %>% distinct(sample) %>% nrow()
diatom %>% distinct(sample) %>% nrow()
syn %>% distinct(sample) %>% nrow()
74+19+9

#this is the correct number of clusters
biomG %>% distinct(OTU_ID) %>% nrow()
tax %>% distinct(`Feature ID`) %>% nrow()

nrow(biomG)

#adds taxonomy to biomG
biomG <- biomG %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))

nrow(biomG)

biomG %>% filter(is.na(class))

biomG %>% filter(is.na(total))
biomG %>% filter(is.na(prop))

#these diatom samples don't have props that add up to 1 because the prop of each sequence in each 
#diatom sample was calculated before I excluded the mitochondria and chloroplast sequences
biomG %>% group_by(sample.y) %>% summarize(totalProp = sum(prop)) %>% filter(totalProp != 1)


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


#these are the clusters arranged by the number of sequences that they have
biomG %>% group_by(OTU_ID) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(desc(n))

#these are the clusters and their taxonomies that are found in multiple datasets
biomG %>% group_by(OTU_ID, class) %>% distinct(sampleSource) %>% summarize(n = n()) %>% arrange(desc(n)) %>% filter(n > 1)

#there are 175 clusters
biomG %>% distinct(OTU_ID) %>% summarize(n = n())
biomG %>% group_by(sampleSource) %>% distinct(OTU_ID) %>% summarize(n = n())
#11 clusters are found in the sean, syn and diatom datasets
biomG %>% group_by(OTU_ID) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n == 3) %>% ungroup() %>% distinct(OTU_ID) %>% nrow()
#31 clusters were found in both sean and syn datasets
biomG %>% filter(sampleSource != "diatom") %>% group_by(OTU_ID) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n == 2) %>% ungroup() %>% distinct(OTU_ID) %>% nrow()
#12 clusters were found in both diatom and syn datasets
biomG %>% filter(sampleSource != "sean") %>% group_by(OTU_ID) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n == 2) %>% ungroup() %>% distinct(OTU_ID) %>% nrow()
#13 clusters were found in both diatom and sean datasets
biomG %>% filter(sampleSource != "syn") %>% group_by(OTU_ID) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n == 2) %>% ungroup() %>% distinct(OTU_ID) %>% nrow()


#these are the distinct taxonomies of the clusters that were found in the sean, syn and diatom datasets
biomG %>% group_by(OTU_ID, class, Taxon) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n == 3) %>% ungroup() %>% distinct(Taxon) %>% arrange(Taxon)

#gets the clusters that were found in just one of the sample sources
foundInOneSource <- biomG %>% group_by(OTU_ID) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n == 1) %>% distinct(OTU_ID)

#141 clusters were found in just one data source
nrow(foundInOneSource)

#85 clusters were just found in sean samples
biomG %>% semi_join(foundInOneSource, by = c("OTU_ID")) %>% filter(sampleSource == "sean") %>% distinct(OTU_ID) %>% summarize(n = n())

#41 clusters were just found in syn samples
biomG %>% semi_join(foundInOneSource, by = c("OTU_ID")) %>% filter(sampleSource == "syn") %>% distinct(OTU_ID) %>% summarize(n = n())

#15 clusters were just found in diatom samples
biomG %>% semi_join(foundInOneSource, by = c("OTU_ID")) %>% filter(sampleSource == "diatom") %>% distinct(OTU_ID) %>% summarize(n = n())

85+41+15
141+13+12+31-11-11

biomG %>% filter(str_detect(class, "Synech|synech")) %>% distinct(class)

#these are the Synechococcophycidea clusters arranged by how many sequences they had
biomG %>% filter(class == "Synechococcophycideae") %>% group_by(OTU_ID) %>% distinct(sample) %>% summarize(n = n()) %>% arrange(desc(n))

#these are the Synechococcophycidea clusters that were found in multiple datasets
biomG %>% filter(class == "Synechococcophycideae") %>% group_by(OTU_ID, class) %>% distinct(sampleSource) %>% summarize(n = n()) %>% filter(n > 1)

#this is the number of sequences in each dataset
biomG %>% group_by(sampleSource) %>% distinct(sample) %>% summarize(n = n())

#makes a variable for whether sample is a Pro or Syn culture
biomG <- biomG %>% mutate(culture = ifelse(str_detect(sample.y, "SYN|WH|9220|S9503|8102|RS9916") & sampleSource == "sean", "Synechococcus", ifelse(sampleSource == "sean", "Prochlorococcus", NA)))

biomG %>% distinct(sampleSource, culture)

nrow(biomG)

#for sean seqs, if the class is Synechococcophycideae, adds whether the sequence is Pro or Syn
biomG <- biomG %>% mutate(class = ifelse(class == "Synechococcophycideae" & sampleSource == "sean", str_c(class, culture, sep = " "), class))

nrow(biomG)

biomG %>% filter(sampleSource == "syn") %>% filter(class == "Synechococcophycideae") %>% distinct(Taxon)
biomG %>% filter(sampleSource == "syn") %>% filter(class == "Synechococcophycideae") %>% distinct(sample)

#biomG %>% filter(sampleSource == "syn" & class == "Synechococcophycideae") %>% View()
biomG %>% filter(sampleSource == "syn" & class == "Synechococcophycideae") %>% distinct(Taxon)
biomG %>% filter(sampleSource == "syn" & class == "Synechococcophycideae" & !str_detect(Taxon, "g__Synechococcus")) %>% distinct(Taxon)
#biomG %>% filter(sampleSource == "syn" & class == "Synechococcophycideae" & !str_detect(Taxon, "g__Synechococcus")) %>% View()

#gets the distinct syn culture sequences that were identified as Synechococcophycideae
#so I can blast them against other sequences in analyzeSeanAndAEMSynSequences.txt
#biomG %>% filter(sampleSource == "syn" & class == "Synechococcophycideae") %>% distinct(sample) %>% 
  #mutate(sample = str_c(sample, "_ending")) %>% 
  #write.table("synSamples_SynechococcophycideaeSeqs.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



###they are all Syn cultures in the AEM paper so I am going to assume all of the Synechococcophycideae sequences in the AEM paper 
###are Synechococcus sequences

biomG %>% filter(sampleSource == "syn") %>% filter(class == "Synechococcophycideae") %>% distinct(Taxon)

#for syn culture seqs, if the class is Synechococcophycideae, adds whether the sequence is Pro or Syn
biomG <- biomG %>% mutate(class = ifelse(class == "Synechococcophycideae" & sampleSource == "syn", str_c(class, "Synechococcus", sep = " "), class))

biomG %>% filter(str_detect(class, "Synechococcophycideae")) %>% distinct(sampleSource, culture, class)

#loads rooted tree
#this came from analysisDiatomAEMSynSean.txt (checked)
tree <- read_tree("exported-rooted-rep-seqs-diatomAEMSynSean-dn-1_tree/tree.nwk")

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

#makes total abundance of each cluster in each sample into a dataframe
seqSpreadSumm <- seqSpreadSumm %>% ungroup() %>% as.data.frame()

head(seqSpreadSumm)

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
sampleTree <- unifracRes %>% hclust(method="ward.D")

#plot_heatmap(physeq, taxa.label="shortTaxa")

plot(sampleTree)

sampleTree$labels

sampleTree$order

sampleTree$labels[sampleTree$order]

#puts seqG sample variable in the same order as unifrac clustering of samples
biomG$sample.y <- factor(biomG$sample.y, levels = sampleTree$labels[sampleTree$order])

levels(biomG$sample.y)

biomG %>% filter(is.na(class))
biomG %>% filter(class == "")

###why do these syn culture samples have multiple syn clusters???
###it seems that there are multiple syn strains coexisting in cultures
biomG %>% filter(str_detect(class, "Synechococcophycideae")) %>% filter(abundance > 0) %>% group_by(sample.y) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% filter(n > 1)
#biomG %>% filter(str_detect(class, "Synechococcophycideae")) %>% filter(prop > .002) %>% group_by(sample.y) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% filter(n > 1)
multipleSyn <- biomG %>% filter(str_detect(class, "Synechococcophycideae")) %>% group_by(sample.y) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% filter(n > 1)
multipleSyn
biomG %>% semi_join(multipleSyn, by = c("sample.y")) %>% filter(str_detect(class, "Synechococcophycideae")) %>% group_by(OTU_ID, sample.y) %>% summarize(abundance = sum(abundance)) %>% arrange(sample.y, desc(abundance)) %>% View()
#biomG %>% semi_join(multipleSyn, by = c("sample.y")) %>% filter(str_detect(class, "Synechococcophycideae")) %>% filter(prop > .002) %>% group_by(OTU_ID, sample.y) %>% summarize(abundance = sum(abundance)) %>% arrange(sample.y, desc(abundance)) %>% View()


#there aren't chloroplast or mitochondria sequences after excluding prop < .002 sequences
biomG %>% filter(str_detect(Taxon, "chloroplast|Chloroplast|mitochondria|Mitochondria"))

#samplesWithChloro <- biomG %>% filter(str_detect(Taxon, "chloroplast|Chloroplast|mitochondria|Mitochondria")) %>% distinct(sample.y)

#all of the samples that have chloroplast sequences also have Syn sequences
#samplesWithChloro
#biomG %>% semi_join(samplesWithChloro, by = c("sample.y")) %>% filter(abundance > 0) %>% filter(str_detect(class, "Synechococcophycideae")) %>% distinct(sample.y, class)

biomG %>% group_by(sample.y) %>% summarize(totalProp = sum(prop)) %>% filter(totalProp != 1)
biomG %>% distinct(sample) %>% nrow()
biomG %>% distinct(OTU_ID) %>% nrow()
biomG %>% distinct(sample.y) %>% nrow()

#colors <- c("cyan3", "sienna4", "olivedrab3", "mediumvioletred", "gray58", "lightsteelblue1", "springgreen4", "hotpink1", "red", "tan3", "blueviolet", "plum", "goldenrod1", "black", "blue1", "lightcoral", "khaki2", "gray92", "springgreen", "pink", "darkorange1", "yellow")

order <- c("A712011", "Acidimicrobiia", 
           "Actinobacteria", "Alphaproteobacteria", "Betaproteobacteria", "Cytophagia", "Deltaproteobacteria", 
           "Flavobacteriia", "Gammaproteobacteria", "[Leptospirae]", "OM190", "Phycisphaerae", "[Rhodothermi]", 
           "[Saprospirae]", "Synechococcophycideae Prochlorococcus", 
           "Synechococcophycideae Synechococcus", "Unclassified at class level")          

length(order)
biomG %>% distinct(class) %>% nrow()

biomG$class <- factor(biomG$class, levels = order)
levels(biomG$class)

colors <- c("cyan3", "khaki2", "olivedrab3", "mediumvioletred", "gray58", "springgreen4", "hotpink1", "red", "tan3", "blueviolet", "lightsteelblue1", "black", "blue1", "lightcoral", "springgreen", "pink", "darkorange1")

length(colors)

biomG %>% ggplot(aes(x = sample.y, y = prop, fill = class)) + geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y = "Relative abundance", fill = "") + scale_fill_manual(values = colors) 

ggsave("diatomAEMSynAndSeanTaxonomyStackedBarPlot_withProSyn.png", dpi = 600, height = 6, width = 16)
ggsave("diatomAEMSynAndSeanTaxonomyStackedBarPlot_withProSyn.svg", dpi = 600, height = 6, width = 16)



###checked everything above


#from https://stackoverflow.com/questions/14118033/horizontal-dendrogram-in-r-with-labels
#convert cluster object to use with ggplot
dendr <- dendro_data(sampleTree, type="rectangle") 

#gets the ordering of the samples on the tree
colors <- label(dendr)
head(colors)
colors <- colors$label
colors <- colors %>% as.data.frame()
colnames(colors) <- "label"
head(colors)

#gets the dataset that each sample is from
sampleSource <- biomG %>% distinct(sample.y, sampleSource)
nrow(sampleSource)
head(sampleSource)

str(colors)
colors$label <- as.character(colors$label)
str(colors)

#adds dataset that samples came from to the ordering of the samples on the tree
colors <- colors %>% left_join(sampleSource, by = c("label" = "sample.y"))

nrow(colors)
colors %>% filter(is.na(sampleSource))

#assigns a color to each distinct data source
colors <- colors %>% mutate(color = ifelse(sampleSource == "sean", "red", ifelse(sampleSource == "syn", "blue", "black")))

colors %>% distinct(sampleSource, color)

head(colors)
tail(colors)
plot(sampleTree)

#gets a list of the colors so that I can plot the labels with the color corresponding to each sample's datasource
toPlot <- colors %>% select(color)
toPlot <- toPlot$color
c(toPlot)
head(colors)
tail(colors)

#from https://stackoverflow.com/questions/14118033/horizontal-dendrogram-in-r-with-labels
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0, color = toPlot), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank()) + theme(legend.position = "none")

ggsave("diatomAEMSynAndSeanUnifracTree_withProSyn.png", dpi = 600, width = 6, height = 12)


plot(as.dendrogram(sampleTree), horiz = T)
plot(sampleTree)






###without Pro and Syn

biomG %>% distinct(class) %>% arrange(class)

#makes biomG dataframe without Pro or Syn
noProSyn <- biomG %>% filter(class != "Synechococcophycideae Prochlorococcus") %>% filter(class != "Synechococcophycideae Synechococcus")

noProSyn %>% distinct(class)

head(noProSyn)

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

head(noProSyn)

#gets just the necessary variables to spread the abundance data
seqSpread <- noProSyn %>% select(OTU_ID, sample.y, abundance)

head(seqSpread)

seqSpread %>% group_by(OTU_ID, sample.y) %>% summarize(n = n()) %>% filter(n > 1)

#gets the total abundance of each cluster in each sample
seqSpreadSumm <- seqSpread %>% group_by(OTU_ID, sample.y) %>% summarize(sumAbundance = sum(abundance))

head(seqSpreadSumm)

seqSpreadSumm %>% group_by(OTU_ID, sample.y) %>% summarize(n = n()) %>% filter(n > 1)

#makes seqSpreadSumm into a dataframe
seqSpreadSumm <- seqSpreadSumm %>% ungroup() %>% as.data.frame()

head(seqSpreadSumm)

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
sampleTree <- unifracRes %>% hclust(method="ward.D")

#plot_heatmap(physeq, taxa.label="shortTaxa")

plot(sampleTree)

sampleTree$labels

sampleTree$order


sampleTree$labels[sampleTree$order]

#puts seqG sample variable in the same order as unifrac clustering of samples
noProSyn$sample.y <- factor(noProSyn$sample.y, levels = sampleTree$labels[sampleTree$order])

levels(noProSyn$sample.y)

noProSyn %>% group_by(sample.y) %>% summarize(totalProp = sum(propWithoutProSyn)) %>% filter(totalProp != 1)

noProSyn %>% distinct(class)

order <- c("A712011", "Acidimicrobiia", 
           "Actinobacteria", "Alphaproteobacteria", "Betaproteobacteria", "Cytophagia", "Deltaproteobacteria", 
           "Flavobacteriia", "Gammaproteobacteria", "[Leptospirae]", "OM190", "Phycisphaerae", "[Rhodothermi]", 
           "[Saprospirae]", "Unclassified at class level")  

order %>% length()

noProSyn %>% distinct(class) %>% nrow()

#colors <- c("cyan3", "sienna4", "olivedrab3", "mediumvioletred", "gray58", "lightsteelblue1", "springgreen4", "hotpink1", "red", "tan3", "blueviolet", "plum", "goldenrod1", "black", "blue1", "lightcoral", "khaki2", "gray92", "darkorange1", "yellow")

colors <- c("cyan3", "khaki2", "olivedrab3", "mediumvioletred", "gray58", "springgreen4", "hotpink1", "red", "tan3", "blueviolet", "lightsteelblue1", "black", "blue1", "lightcoral", "darkorange1")

length(colors)

noProSyn$class <- factor(noProSyn$class, levels = order)

levels(noProSyn$class)

noProSyn %>% ggplot(aes(x = sample.y, y = propWithoutProSyn, fill = class)) + geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y = "Relative abundance", fill = "") + scale_fill_manual(values = colors)

ggsave("diatomAEMSynAndSeanTaxonomyStackedBarPlot_noProSyn.png", dpi = 600, height = 6, width = 16)
ggsave("diatomAEMSynAndSeanTaxonomyStackedBarPlot_noProSyn.svg", dpi = 600, height = 6, width = 16)


noProSyn %>% group_by(sampleSource) %>% distinct(OTU_ID) %>% summarize(n = n())


#from https://stackoverflow.com/questions/14118033/horizontal-dendrogram-in-r-with-labels
#convert cluster object to use with ggplot
dendr <- dendro_data(sampleTree, type="rectangle") 

#gets the ordering of the samples on the tree
colors <- label(dendr)
head(colors)
colors <- colors$label
colors <- colors %>% as.data.frame()
colnames(colors) <- "label"
head(colors)

#gets the dataset that each sample is from
sampleSource <- biomG %>% distinct(sample.y, sampleSource)
nrow(sampleSource)
head(sampleSource)

str(colors)
colors$label <- as.character(colors$label)
str(colors)

#adds dataset that samples came from to the ordering of the samples on the tree
colors <- colors %>% left_join(sampleSource, by = c("label" = "sample.y"))

nrow(colors)
colors %>% filter(is.na(sampleSource))

#assigns a color to each distinct data source
colors <- colors %>% mutate(color = ifelse(sampleSource == "sean", "red", ifelse(sampleSource == "syn", "blue", "black")))

colors %>% distinct(sampleSource, color)

head(colors)
tail(colors)
plot(sampleTree)

#gets a list of the colors so that I can plot the labels with the color corresponding to each sample's datasource
toPlot <- colors %>% select(color)
toPlot <- toPlot$color
c(toPlot)
head(colors)
tail(colors)

#from https://stackoverflow.com/questions/14118033/horizontal-dendrogram-in-r-with-labels
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0, color = toPlot), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank()) + theme(legend.position = "none")

ggsave("diatomAEMSynAndSeanUnifracTree_noProSyn.png", dpi = 600, width = 6, height = 12)


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
ggsave("diatomAEMSynAndSeanUnifracTree_noProSyn.svg", dpi = 600, width = 6, height = 12)



plot(as.dendrogram(sampleTree), horiz = T)

plot(sampleTree)

png("diatomAEMSynAndSeanUnifracTree_noProSyn_simple.png")
plot(sampleTree, labels = FALSE)
dev.off()

svg("diatomAEMSynAndSeanUnifracTree_noProSyn_simple.svg")
plot(sampleTree, labels = FALSE)
dev.off()




