library(tidyverse)

setwd("~/Dropbox (MIT)/clusterAllisonSeanIllumina16S/")

#this came from clusterSeanAllisonIllumina16SSeqs_excludingLowAbundanceAllison.txt (checked)
biom <- read.table("feature-table-fromBiom_ExcludingLowCountsAllison.tsv") 

#these are feature/sample names
#these came from feature-table-fromBiom_ExcludingLowCountsAllison.tsv (checked)
names <- "V1	07d0bb5d023b91f79a21fa1f08de76da9475afa7	16e13db71aa93d0c849a0a4c8a5c14104ef149ea	1acb4b1dcfdf5be27c536949eb15515631711d41	27fedabca0a1e1f5d8d3b279d42c25113299de01	3033e7739060a50de88c2b06168e26d681fca9c0	39c7cac93b7e1913f21cbcf959ec738710dc3381	45ac47fbf72fbdad0846bcf5f5b804f65f32fd54	47d3b61a3ec945aecde3475e6d98986b35e8c6f3	48d7e3444f85af8c2bd1f10624778437d82e583d	4a54d389b1a158f478d861637299bf6ae04aae80	4d42d41088f08ee2db9f70375560715a507c91c5	53e028c6f56c4210d9312ed07831258ad003a9a2	5715644b7f4888973e7f8c1456391a084587d576	57872c999308bfec89b949c43c4a4b51aafdb0f3	5865070a00b3df85a584690732ff897c19063217	5c8141fea9c85ab22abad95fbf8c1b66873d18de	6b28f251bd68700bcbe4a7e4986b19d7c825bcf6	7176f58de210a7cea377d555b1d247c67b69855f	7f131455a0d7d583a795db8bf875920e1138077f	8128f23d8559bee53a0a97f540b4b47c77775982	8785e4a7f77eacd1d9863a2e9db9099f7e2314bf	9773ae42f9c1076c7bb03948d9dbc4ab89461f62	9b5f9ebe90d1eaaa37237e4a51152aa79c966e84	b1d24a3d6bc8f67e2e20782a80fc9a1981a6506b	b28441d052be31e9e23600490c26e647125695f3	b92cf2c1956ebbe6b0e7a4b7edd646d4dd84f9dc	c4f4e9f84d490e39681ebe5ee2b11ad0a3db4dfd	ca0cace32a7c6c5813227f0a06c9f227ca845671	d089cd89ef62f29a2a3399c9fe6344baa260aa45"

#makes sample names into list
namesNew <- str_split(names, pattern = "\t")
namesNew <- namesNew[[1]]

#same length which is good
ncol(biom)
length(namesNew)

#makes the column names the sample names
colnames(biom) <- namesNew

#there are 29 clustered features
nrow(biom)

#there are 29 sample/sequence columns
ncol(biom)


#gathers into long format
biomG <- biom %>% gather(2:30, key = "sample", value = "present")

#there are 29 clustered features and 29 samples
biomG %>% distinct(V1) %>% nrow()
biomG %>% distinct(sample) %>% nrow()

biomG %>% filter(present > 0) %>% distinct(V1) %>% nrow()
biomG %>% filter(present > 0) %>% distinct(sample) %>% nrow()

#there are 29 times a sequence is found in a sample 
#because this is just presence/absence because I imported features to qiime2
#I need to connect the features to their abundances in each sample
biomG %>% filter(present > 0)

#old <- read.table("feature-table-fromBiom.tsv")

#oldG <- old %>% gather(2:60, key = "sample", value = "present")

#this is the same as the nonclustered results which is good!
#oldG %>% filter(present > 0)

#loads taxonomy of clustered features
taxa <- read.table("taxonomySeanAllisonFeatures_ExcludingLowCountsAllison.tsv", fill = TRUE)

#replaces NAs with empty strings so that I can join all of the taxonomy variables together
taxa[is.na(taxa)]
taxa[is.na(taxa)] <- ""

#combines all of the taxonomy variables together into one variable
taxa <- taxa %>% mutate(taxon = str_c(V2, V3, V4, V5, V6, V7, V8, V9, sep = " "))

#gets rid of unnecessary taxonomy variables
taxa <- taxa %>% select(V1, taxon)

#gets rid of row with variable names
taxa <- taxa[-1,]

#there are 29 features
nrow(taxa)

#fixes column names
colnames(taxa) <- c("feature", "taxon")

#gets rid of semicolons
taxa$taxon <- str_replace_all(taxa$taxon, ";", "")

str_extract_all(taxa$taxon, "[a-z]__") %>% unique()

#gets rid of letters with underscores
taxa$taxon <- str_replace_all(taxa$taxon, "[a-z]__", "")

str_extract_all(taxa$taxon, "[0-9]*\\.[0-9]*")

#gets rid of probabilities
taxa$taxon <- str_replace_all(taxa$taxon, "[0-9]*\\.[0-9]*", "")

str_extract_all(taxa$taxon, "\\s*$") %>% unique()

#gets rid of extra spaces at the end of taxonomies
taxa$taxon <- str_replace_all(taxa$taxon, "\\s*$", "")

#makes a variable for the shortened taxonomy
taxa <- taxa %>% mutate(shortTaxa = str_extract(taxon, "[A-Z][a-z,\\s]*$"))

taxa

#makes shortTaxa that didn't get extracted
taxa <- taxa %>% mutate(shortTaxa = ifelse(is.na(shortTaxa), str_extract(taxon, "Alteromonadaceae ZD0117"), shortTaxa))

taxa

nrow(biomG)

#adds taxonomy to biomG
biomG <- biomG %>% left_join(taxa, by = c("V1" = "feature"))

nrow(biomG)

biomG %>% filter(is.na(shortTaxa))

#from getSeanAllisonFeaturesAfterExcludingAllisonLowAbundance.R (checked)
summarySpread <- read_csv("summarySpread_ExcludingLowCountsAllison.csv")

#gets rid of unnecessary variables
summarySpread <- summarySpread %>% select(-c(propSean, totalAllisonSean, foundInSeanAllison, sumPro, logNormSeanAllisonByPro))

#gets the number of Pro sequences in each Sean and Allison sample
proAbundance <- summarySpread %>% distinct(sample.y, proAbundanceSean, proAbundanceAllison)

#puts proAbundance into long form
proAbundance <- proAbundance %>% gather(2:3, key = "source", value = "proAbundance")

#makes source match summary
proAbundance$source <- str_extract(proAbundance$source, "Sean|Allison")
proAbundance$source <- str_to_lower(proAbundance$source)

#gets rid of unnecessary variables from summarySpread so I can gather 
#into long form
summarySpread <- summarySpread %>% select(-c(proAbundanceAllison, proAbundanceSean))

#there are 29 features which is correct
summarySpread %>% distinct(V1) %>% nrow()

#gathers sequence data into long form
summary <- summarySpread %>% gather(allison:sean, key = "source", value = "sumAbundance")

#there are still 29 features which is good
summary %>% distinct(V1) %>% nrow()

nrow(summary)

#adds Sean and Allison Pro read abundance to summary by sample and source
summary <- summary %>% left_join(proAbundance, by = c("sample.y", "source"))

nrow(summary)

summary %>% filter(is.na(proAbundance))

#makes a variable for the log of the total number of reads of each feature 
#divided by the number of Pro in each sample by source
#takes log10 because this is what I did for Allison seqs previously and Sean said to take 
#log base 10 or base 2
summary <- summary %>% mutate(logNormTotalByPro = log10(sumAbundance/proAbundance))

#summarySpread$foundInSeanAllison <- factor(summarySpread$foundInSeanAllison, levels = c(TRUE, FALSE))

order <- rev(c("48d7e3444f85af8c2bd1f10624778437d82e583d",
               "57872c999308bfec89b949c43c4a4b51aafdb0f3",
               "9773ae42f9c1076c7bb03948d9dbc4ab89461f62",
               "ca0cace32a7c6c5813227f0a06c9f227ca845671",
               "47d3b61a3ec945aecde3475e6d98986b35e8c6f3",
               "3033e7739060a50de88c2b06168e26d681fca9c0",
               "4d42d41088f08ee2db9f70375560715a507c91c5",
               "8785e4a7f77eacd1d9863a2e9db9099f7e2314bf",
               "d089cd89ef62f29a2a3399c9fe6344baa260aa45",
               "b92cf2c1956ebbe6b0e7a4b7edd646d4dd84f9dc",
               "c4f4e9f84d490e39681ebe5ee2b11ad0a3db4dfd",
               "07d0bb5d023b91f79a21fa1f08de76da9475afa7",
               "5715644b7f4888973e7f8c1456391a084587d576",
               "53e028c6f56c4210d9312ed07831258ad003a9a2",
               "39c7cac93b7e1913f21cbcf959ec738710dc3381",
               "5c8141fea9c85ab22abad95fbf8c1b66873d18de",
               "5865070a00b3df85a584690732ff897c19063217",
               "7176f58de210a7cea377d555b1d247c67b69855f",
               "7f131455a0d7d583a795db8bf875920e1138077f",
               "16e13db71aa93d0c849a0a4c8a5c14104ef149ea",
               "1acb4b1dcfdf5be27c536949eb15515631711d41",
               "45ac47fbf72fbdad0846bcf5f5b804f65f32fd54",
               "b1d24a3d6bc8f67e2e20782a80fc9a1981a6506b",
               "9b5f9ebe90d1eaaa37237e4a51152aa79c966e84",
               "27fedabca0a1e1f5d8d3b279d42c25113299de01",
               "6b28f251bd68700bcbe4a7e4986b19d7c825bcf6",
               "8128f23d8559bee53a0a97f540b4b47c77775982",
               "4a54d389b1a158f478d861637299bf6ae04aae80",
               "b28441d052be31e9e23600490c26e647125695f3"))

#same as number of features which is good
length(order)

#puts features in the same order on the heatmap as on the tree
summary$V1 <- factor(summary$V1, levels = order)

summary %>% filter(is.na(logNormTotalByPro))

summary %>% filter(sumAbundance == 0)

#excludes rows for which the abundance is 0 because you can't log(0) is -Inf
#and if I exclude abundance 0 rows, the cell for these rows on the heatmaps will be white 
#which is what I want
summary <- summary %>% filter(sumAbundance != 0)

summary %>% ggplot(aes(sample.y, V1)) + geom_tile(aes(fill = logNormTotalByPro)) + facet_wrap(~source) +
  labs(x = "", fill = "", y = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_gradient(low = "white", high = "red") + scale_color_manual(values = c("red", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.ticks.x = element_blank(), 
        axis.text=element_text(size=12, color = "black"))

ggsave("heatmapWithoutLowAbundanceAllison_sideBySide.png", dpi = 600, height = 6, width = 8)

#ungroup dataframe
summary <- summary %>% ungroup() 

#pro and syn features have been excluded
summary %>% filter(str_detect(taxon, "Proc|proc|Syn|syn"))

#calculates the number of each feature divided by the number of Pro sequences in each sample
summary <- summary %>% mutate(sumAbundanceNormByProAbundance = sumAbundance/proAbundance)

summary %>% filter(is.na(sumAbundanceNormByProAbundance))

head(summary)

#gets rid of unnecessary variables
forCorr <- summary %>% select(V1,source,sample.y,sumAbundanceNormByProAbundance)

head(forCorr)

#makes a variable that includes both the source and the sample
forCorr <- forCorr %>% mutate(sourceSample = str_c(source, sample.y, sep = " "))

head(forCorr)

#gets rid of unnecessary variables
forCorr <- forCorr %>% select(-c(source, sample.y))

head(forCorr)

#spreads the data that each feature has a sumAbundanceNormByProAbundance value for each Sean 
#and Allison sample
forCorr %>% group_by(V1) %>% summarize(n = n()) %>% filter(n != 2)
head(forCorr)
forCorr <- forCorr %>% spread(key = V1, value = sumAbundanceNormByProAbundance, fill = 0)
head(forCorr)

#makes the rownames the sample IDs
rownames(forCorr) <- forCorr$sourceSample
forCorr <- forCorr %>% select(-sourceSample)

head(forCorr)

#gets the transpose of forCorr so that the sample IDs are column names because that 
#is what I want to calculate correlations between
forCorr <- t(forCorr)

dim(forCorr)

library(psych)

#gets the spearman correlation between each pair of non-Pro, non-Syn features
#this is based on each feature's normalized abundance across the samples

#also, gets the pvalues for each correlation so I can calculate qvalues and then 
#get pairs with FDR < .1
#does not adjust pvalue for multiple tests because you calculate q values based 
#on unadjusted pvalues 
#from https://stackoverflow.com/questions/13112238/a-matrix-version-of-cor-test
test <- corr.test(forCorr, method = "spearman", adjust = "none")

library(qvalue)

###I can't calculate the qvalue because there are two few values to estimate correctly 
###which I learned from https://github.com/StoreyLab/edge/issues/13

#makes a dataframe with the correlations between each pair of features
corDF <- test$r %>% as.data.frame()

#correct dimensions
dim(corDF)

corDF

#clusters the samples based on their correlations
#from https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/hclust
hc <- hclust(dist(corDF), "ward.D2")

plot(hc)

hc$labels

hc$order

hc$labels[hc$order]

corDF

#makes a variable for one set of the samples
corDF <- corDF %>% mutate(sample1 = row.names(corDF))

corDF

#gathers the correlations into long form
corDF <- corDF %>% gather(1:6, key = "sample2", value = "correlation")

head(corDF)

#there 6 samples which is correct
corDF %>% distinct(sample1) %>% nrow()
corDF %>% distinct(sample2) %>% nrow()

#puts the samples in the order that they appear on the dendrogram
corDF$sample1 <- factor(corDF$sample1, levels = hc$labels[hc$order])
corDF$sample2 <- factor(corDF$sample2, levels = hc$labels[hc$order])

levels(corDF$sample1)
levels(corDF$sample2)

corDF %>% arrange(sample1) %>% distinct(sample1)

#plots correlations for each pair of samples
corDF %>% ggplot(aes(x = sample1, y = sample2, fill = correlation)) + geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits=c(-1,1)) + 
  labs(fill = "Spearman correlation based on\nabundance of non-Pro, non-Syn\nASVs normalized by Pro. abundance", x = "", y = "") + 
  theme(legend.title=element_text(size=6), legend.text=element_text(size=6)) + 
  theme(legend.key.size = unit(.4, "cm")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#ggsave("correlationHeatmapWithoutLowAbundanceAllison_combinedSeanAllisonIllumina.png", dpi = 600, height = 6, width = 8)
#ggsave("correlationHeatmapWithoutLowAbundanceAllison_combinedSeanAllisonIllumina.svg", dpi = 600, height = 6, width = 8)


levels(corDF$sample1)
levels(corDF$sample2)
corDF$sample1 <- factor(corDF$sample1, levels = c("allison MED4", "sean MED4", "allison MIT9313", "sean MIT9313", "allison NATL2A", "sean NATL2A"))
corDF$sample2 <- factor(corDF$sample2, levels = c("sean NATL2A", "allison NATL2A", "sean MIT9313", "allison MIT9313", "sean MED4", "allison MED4"))

levels(corDF$sample1)
levels(corDF$sample2)

corDF %>% ggplot(aes(x = sample1, y = sample2, fill = correlation)) + geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits=c(-1,1)) + 
  labs(fill = "Spearman correlation based on\nabundance of non-Pro, non-Syn\nASVs normalized by Pro. abundance", x = "", y = "") + 
  theme(legend.title=element_text(size=6), legend.text=element_text(size=6)) + 
  theme(legend.key.size = unit(.4, "cm")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("correlationHeatmapWithoutLowAbundanceAllison_combinedSeanAllisonIllumina.png", dpi = 600, height = 6, width = 8)
ggsave("correlationHeatmapWithoutLowAbundanceAllison_combinedSeanAllisonIllumina.svg", dpi = 600, height = 6, width = 8)





