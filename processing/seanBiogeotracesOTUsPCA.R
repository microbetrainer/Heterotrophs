library(tidyverse)
library(lubridate)
library(phyloseq)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads the sean and biogeotraces sequences present in each cluster
#from clusterSeanAndBiogeotracesSeqs.txt (checked)
biomG <- read_csv("combinedBiogeotracesBiom.csv")

#28022 clusters (correct)
#1184489 sequences (correct)
biomG %>% distinct(OTU_ID) %>% nrow()
biomG %>% distinct(sample) %>% nrow()
biomG %>% nrow()

str(biomG)

biomG %>% distinct(present)

#each sequence is only in one cluster
biomG %>% filter(present == 1) %>% group_by(sample) %>% summarize(n = n()) %>% filter(n != 1)
biomG %>% group_by(sample) %>% summarize(n = n()) %>% filter(n != 1)

#loads the abundance of each sean sequence in each culture
#this includes pro and syn 
#sequences with prop < .002 have been excluded
#from heatmap_noContamination_noLowAbundance_noMock_noProSyn_without1223_propReads.R
seqG <- read_csv("forSeanComparisonToBiogeotraces.csv")

#this includes Pro and Syn sequences and the JW3 and 1223 samples
seqG %>% distinct(sequence) %>% nrow()

seqG %>% filter(sample == "JW3")
seqG %>% filter(sample == "1223")

#excludes JW3 and 1223
seqG <- seqG %>% filter(sample != "JW3")
seqG <- seqG %>% filter(sample != "1223")

#includes Pro and Syn sequences
seqG %>% distinct(sequence) %>% nrow()
head(seqG)

#there is only one row for each sequence, sample pair
seqG %>% group_by(sample, sequence) %>% summarize(n = n()) %>% filter(n > 1)

#1213 is 1313
#WH7803 is WH7801
#SYN1320 is SYN1220
#9201.2 is 9311
#9201.1 is 9201

#AS9601 is correct ID rather than ASN9601
#SYN9503 is correct ID rather than S9503

#fixes sample IDs that were mismatched
#Sean and I determined which sample IDs were mismatched over slack
seqG <- seqG %>% mutate(sample = ifelse(sample == "1213", "1313", sample))
#seqG <- seqG %>% mutate(sample = ifelse(sample == "WH7803", "WH7801", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "SYN1320", "SYN1220", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "9201.2", "9311", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "9201.1", "9201", sample))

seqG <- seqG %>% mutate(sample = ifelse(sample == "ASN9601", "AS9601", sample))
seqG <- seqG %>% mutate(sample = ifelse(sample == "S9503", "SYN9503", sample))

seqG$sample <- ifelse(seqG$sample == "8102", "SYNCLADEXVI", seqG$sample)

seqG %>% head()

#correct number of samples
seqG %>% distinct(sample) %>% nrow()

#loads the sequences after excluding pro and syn sequences
#sequences with prop < .002 have been excluded
#from getSequencesForTree_noContamination_noProSyn_withoutSample1223.R
seqsNoProSyn <- read.table("sequencesForTree_noContamination_noMock_noProSyn_no1223.txt")

#includes Pro and Syn
seqG %>% distinct(sequence) %>% nrow()

#doesn't include Pro and Syn 
head(seqsNoProSyn)
seqsNoProSyn %>% distinct(V1) %>% nrow()

head(seqG)
head(seqsNoProSyn)

#excludes Pro and Syn sequences from seqG
seqG <- seqG %>% semi_join(seqsNoProSyn, by = c("sequence" = "V1"))

#correct number of sequences and samples
seqG %>% distinct(sequence) %>% nrow()
seqG %>% distinct(sample) %>% nrow()

#biomG includes Pro and Syn sequences in sean samples
biomG %>% filter(str_detect(sample, "sean")) %>% nrow()

biomG %>% filter(str_detect(sample, "sean")) %>% head()

##makes the sean sequence IDs in biomG match seqG
biomG <- biomG %>% mutate(sample = ifelse(str_detect(sample, "sean"), str_replace(sample, "X", ""), sample))
biomG %>% filter(str_detect(sample, "sean")) %>% head()

biomG <- biomG %>% mutate(sample = ifelse(str_detect(sample, "sean"), str_c("contig_", sample), sample))
biomG %>% filter(str_detect(sample, "sean")) %>% head()

biomG <- biomG %>% mutate(sample = ifelse(str_detect(sample, "sean"), str_replace(sample, "sean", ""), sample))

biomG %>% filter(str_detect(sample, "contig")) %>% head()

#these are the sean sequences
#this includes Pro and Syn
biomG %>% filter(str_detect(sample, "contig")) %>% nrow()

biomG %>% filter(str_detect(sample, "contig")) %>% head()
seqG %>% head()

#finds the sean sequences in biomG that do not have matches in seqG
#these are the sean pro and syn sequences
biomG %>% filter(str_detect(sample, "contig")) %>% left_join(seqG, by = c("sample" = "sequence")) %>% filter(is.na(abundance)) %>% nrow()
268-33

#these are the clusters arranged by how many sequences they have across 
#biogeotraces and sean samples
biomG %>% group_by(OTU_ID) %>% summarize(n = n()) %>% arrange(desc(n))

#these are the clusters in sean samples arranged by how many sequences they have 
#across sean samples
biomG %>% filter(str_detect(sample, "contig")) %>% group_by(OTU_ID) %>% summarize(n = n()) %>% arrange(desc(n))

head(seqG)

#excludes unnecessary variable
seqG <- seqG %>% select(-abundance)

head(biomG)
head(seqG)

#adds the cluster each sean sequence is assigned to to the abundance of
#sean sequences across sean samples
#with Pro and Syn
nrow(biomG)
biomG <- biomG %>% left_join(seqG, by = c("sample" = "sequence"))
nrow(biomG)

head(biomG)

colnames(biomG)
colnames(biomG)[2]
colnames(biomG)[2] <- "sequence"
colnames(biomG)[4]
colnames(biomG)[4] <- "sample"

#these are the sean samples and sequences
head(biomG)
biomG %>% filter(!is.na(sample)) %>% distinct(sample) %>% nrow()
biomG %>% filter(!is.na(sample)) %>% distinct(sequence) %>% nrow()

#makes a variable for whether the sequence is a sean culture or biogeotraces sequence
biomG <- biomG %>% mutate(source = ifelse(str_detect(sequence, "biogeotraces"), "biogeotraces", "sean"))

biomG %>% filter(is.na(source))

biomG %>% head()

#extracts the biogeotraces sample from the biogeotraces sequence ID
biomG <- biomG %>% mutate(sample = ifelse(source == "biogeotraces", str_extract(sequence, "[a-zA-Z0-9]*"), sample))

biomG %>% filter(source == "sean") %>% filter(!is.na(sample)) %>% distinct(sample) %>% nrow()
biomG %>% filter(source == "sean") %>% filter(!is.na(sample)) %>% distinct(sequence) %>% nrow()

biomG %>% filter(is.na(sample)) %>% distinct(source)

head(biomG)

#these are sean culture pro and syn sequences
biomG %>% filter(is.na(sample)) %>% distinct(sequence) %>% nrow()

#excludes pro and syn sequences that were in sean samples
biomG <- biomG %>% filter(!is.na(sample))

biomG %>% filter(is.na(sample))

#correct
biomG %>% group_by(source) %>% distinct(sequence) %>% summarize(n = n())
biomG %>% group_by(source) %>% distinct(sample) %>% summarize(n = n())

head(biomG)

#there are multiiple rows for OTU, sample pairs if there multiple sequences in the OTU 
#that are present in the sample
biomG %>% group_by(OTU_ID, sample) %>% summarize(n = n()) %>% filter(n > 1)

head(biomG)

#makes it so there is only one row for each OTU, sample pair
biomG <- biomG %>% distinct(OTU_ID, present, sample, source)

#there is just one row for each OTU, sample pair
biomG %>% group_by(OTU_ID, sample) %>% summarize(n = n()) %>% filter(n > 1)

#loads the taxonomy of the biogeotraces and sean OTUs
#from clusterSeanAndBiogeotracesSeqs.txt (checked)
tax <- read_tsv("sequencesForTree_noContamination_noMock_withProSyn_no1223_biogeotracesTaxonomy.tsv")

#excludes first row because it just has comments
tax[1,]
tax <- tax[-1,]
head(tax)

#gets rid of unnecessary variable
tax <- tax %>% select(1,2)

#correct because there are 28022 clusters
nrow(tax)

head(biomG)
head(tax)

#adds taxonomies of clusters to biomG
nrow(biomG)
biomG <- biomG %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(biomG)

head(biomG)
biomG %>% filter(is.na(Taxon))

biomG %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)
biomG %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(Taxon)

#this makes sense because I already excluded the pro and syn sequences in sean samples
#there are still Pro and Syn OTUs in biomG biogeotraces samples
biomG %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(source)

#excludes pro and syn OTUs from biomG
biomG <- biomG %>% filter(!str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc"))

biomG %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)

#gets just the OTU, sample pairs
head(biomG)
phylo_seqSpread <- biomG %>% select(OTU_ID, sample, present)
head(phylo_seqSpread)

phylo_seqSpread %>% distinct(present)

phylo_seqSpread %>% distinct() %>% nrow()
phylo_seqSpread %>% nrow()

#spread phylo_seqSpread so samples are rows and OTUs are columns
head(phylo_seqSpread)
phylo_seqSpread <- phylo_seqSpread %>% spread(key = OTU_ID, value = present, fill = 0)
phylo_seqSpread[1:5,1:5]

dim(phylo_seqSpread)

phylo_seqSpread %>% filter(is.na(sample)) %>% nrow()

#makes rownames the sample IDs
rownames(phylo_seqSpread) <- phylo_seqSpread$sample
phylo_seqSpread <- phylo_seqSpread %>% select(-sample)

phylo_seqSpread[1:5,1:5]
dim(phylo_seqSpread)

#makes phyloseq object out of otu table
otu <- otu_table(phylo_seqSpread, taxa_are_rows = FALSE)

#makes rownames of taxonomy dataframe the OTU IDs
head(tax)
rownames(tax) <- tax$`Feature ID`
head(tax)

#makes phyloseq object out of taxPhyloseq
tax <- tax_table(as.matrix(tax))

#loads rooted tree
#from clusterSeanAndBiogeotracesSeqs.txt (checked)
tree <- read_tree("exported-rooted-seanBiogeotraces-dn-1_tree/tree.nwk")

##makes dataframe of OTUs in each phyloseq object

colnames(otu)
otuDF <- colnames(otu) %>% as.data.frame()
head(otuDF)

tree$tip.label
treeDF <- tree$tip.label %>% as.data.frame()
head(treeDF)

tax %>% rownames()
TAXDF <- tax %>% rownames() %>% as.data.frame()
head(TAXDF)

#doesn't include pro and syn
otuDF %>% nrow()

#includes pro and syn
treeDF %>% nrow()
TAXDF %>% nrow()

#gives column a variable name
colnames(otuDF) <- "V1"
colnames(treeDF) <- "V1"
colnames(TAXDF) <- "V1"


treeDF %>% anti_join(otuDF, by = c("V1")) %>% nrow()
TAXDF %>% anti_join(otuDF, by = c("V1")) %>% nrow()

otuDF %>% anti_join(treeDF, by = c("V1"))
TAXDF %>% anti_join(treeDF, by = c("V1"))

otuDF %>% anti_join(TAXDF, by = c("V1"))
treeDF %>% anti_join(TAXDF, by = c("V1"))

taxa_names(otu)
taxa_names(tax)
taxa_names(tree)


#gets the source of each sample: sean or biogeotraces
head(biomG)
sam_data <- biomG %>% distinct(sample, source)
head(sam_data)
#correct number of samples
sam_data %>% group_by(source) %>% summarize(n = n())

#fixes source values to make them more explanatory in plots
head(sam_data)
sam_data <- sam_data %>% mutate(source = str_replace(source, "biogeotraces", "BioGeoTraces sample"))
sam_data <- sam_data %>% mutate(source = str_replace(source, "sean", "Culture sample"))
sam_data %>% group_by(source) %>% summarize(n = n())

#makes rownames the sample IDs
head(sam_data)
rownames(sam_data) <- sam_data$sample
head(sam_data)
#sam_data <- sam_data %>% select(-sample)
#head(sam_data)

#makes dataframe of the sample IDs from sam_data
rownames(sam_data)
samDF <- rownames(sam_data)
head(samDF)
samDF <- samDF %>% as.data.frame()
head(samDF)
colnames(samDF) <- "V1"
head(samDF)
nrow(samDF)

#makes dataframe of the sample IDs from phylo_seqSpread
rownames(phylo_seqSpread)
otuDF <- rownames(phylo_seqSpread)
head(otuDF)
otuDF <- otuDF %>% as.data.frame()
head(otuDF)
colnames(otuDF) <- "V1"
head(otuDF)

#the sample IDs match between sam_data and phylo_seqSpread
samDF %>% anti_join(otuDF, by = c("V1"))
otuDF %>% anti_join(samDF, by = c("V1"))

#makes phyloseq sample data object
sam <- sample_data(sam_data)

#makes combined phyloseq object out of otu, TAX, tree, and sam
physeq <- phyloseq(otu, tax, tree, sam)

saveRDS(physeq, "physeq_seanBiogeotraces")

#https://groups.google.com/forum/#!topic/qiime-forum/7O8Bpz6T0VU
ordu <- readRDS("ordu_seanBiogeotraces")
plot_ordination(physeq, ordu, color = "source")

ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_pca_OTU.png", dpi = 600, height = 12, width = 12)
ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_pca_OTU.svg", dpi = 600, height = 12, width = 12)

