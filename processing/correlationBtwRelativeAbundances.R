library(tidyverse)
library(readr)
library(readxl)
library(picante)
library(phyloseq)


setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads abundances of each feature in each sample
#this came from heatmap_noContamination_noLowAbundance_noMock_noProSyn_without1223.R
seqG <- read_csv("seqG_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.csv")

#there are 74 samples and 235 features 
#this is after excluding contamination, and Pro and Syn features, and 1223 and JW3
seqG %>% distinct(sample) %>% nrow()
seqG %>% distinct(V1) %>% nrow()

seqG %>% filter(sample == "JW3")
seqG %>% filter(sample == "1223")

#each sample has an abundance for each of the 235 features and eacj 
#feature has an abundance for each of the 74 samples
seqG %>% group_by(sample) %>% summarize(n = n()) %>% filter(n != 235)
seqG %>% group_by(V1) %>% summarize(n = n()) %>% filter(n != 74)

#spread data so I can run cumulative sum normalization and then 
#calculate spearman correlation between each pair of features
seqGS <- seqG %>% spread(key = V1, value = abundance)

#makes seqGS into a dataframe
seqGS <- seqGS %>% as.data.frame()

#makes rownames the sample
rownames(seqGS) <- seqGS$sample
seqGS <- seqGS %>% select(-sample)

#correct number of dimensions
dim(seqGS)


library(metagenomeSeq)

#from https://rdrr.io/bioc/metagenomeSeq/man/newMRexperiment.html:
#count data is representative of the number of reads annotated for a feature
#rows should correspond to features and columns to samples

seqGS[1:5, 1:4]

#makes dataframe for which rows are features and columns are samples
counts <- t(seqGS)

counts[1:5, 1:4]

obj <- newMRexperiment(counts)

#from https://rdrr.io/bioc/metagenomeSeq/man/cumNorm.html and 
#https://support.bioconductor.org/p/64061/
#res <- cumNorm(obj, p = cumNormStat(obj))

#these are the results from cumNorm
#the number for each sample tells me what to normalize the count data by 
#for that sample
#head(normFactors(res))

#normFactors(res) %>% str()
#normFactors(res) %>% length()

#normalizes the number of sequences in each sample by the normalization
#factor outputted by cumNorm
#2 indicates to apply over columns
#from https://www.r-bloggers.com/using-apply-sapply-lapply-in-r/
#norm <- apply(seqGS, 2, function(x) x/normFactors(res))

#seqGS[,1]
#norm[,1]
#seqGS[,1]/norm[,1]
#normFactors(res)

#seqGS[,2]
#norm[,2]
#seqGS[,2]/norm[,2]
#normFactors(res)

#seqGS[3,]
#norm[3,]
#seqGS[3,]/norm[3,]
#normFactors(res)



#from https://www.rdocumentation.org/packages/metagenomeSeq/versions/1.14.0/topics/cumNormMat
#this is the same as how I previously normalized by dividing the number of sequences by the 
#normalization factors, except this then multiplies by 1000
norm <- cumNormMat(obj, p = cumNormStat(obj), sl = 1000)

#none of the resulting values are more than 0 but less than 1
norm[norm < 1 & norm > 0]

str(seqGS)
str(norm)

#makes normalized counts into a dataframe
norm <- norm %>% as.data.frame()

dim(seqGS)
dim(norm)

#gets transpose of norm
norm <- t(norm)

str(norm)

#makes norm into dataframe
norm <- norm %>% as.data.frame()

#same dimensions which is good
dim(seqGS)
dim(norm)

#gets the spearman correlation between each pair of non-Pro, non-Syn features
#this is based on each feature's normalized abundance across the samples
#corDF <- cor(norm, method = "spearman")

library(psych)

#gets the spearman correlation between each pair of non-Pro, non-Syn features
#this is based on each feature's normalized abundance across the samples

#also, gets the pvalues for each correlation so I can calculate qvalues and then 
#get pairs with FDR < .1
#does not adjust pvalue for multiple tests because you calculate q values based 
#on unadjusted pvalues 
#from https://stackoverflow.com/questions/13112238/a-matrix-version-of-cor-test
test <- corr.test(norm, method = "spearman", adjust = "none")

#r is the matrix of correlations
#test$r %>% View()

#p is the unadjusted two tailed probability (pvalue) for each correlation
#test$p %>% View()

#the rownames and colnames in test$p are the same
unique(colnames(test$p) == rownames(test$p))

#from this: http://www.nonlinear.com/support/progenesis/comet/faq/v2.0/pq-values.aspx, 
#makes top half of dataframe and the diagonal NA because
#I want to calculate qvalues based on just unique pairs of sequences, not including 
#sequences compared to themselves (the diagonal)
test$p[upper.tri(test$p, diag = TRUE)] <- NA

#makes dataframe of unadjusted pvalues
unadjustedP <- test$p
str(unadjustedP)
unadjustedP <- unadjustedP %>% as.data.frame()

#correct number of dimensions
dim(unadjustedP)


#makes the rownames (one set of features) into a variable
unadjustedP_g <- unadjustedP %>% mutate(feature1 = row.names(unadjustedP))

#gathers the unadjusted pvalues into long form
unadjustedP_g <- unadjustedP_g %>% gather(1:235, key = "feature2", value = "unadjustedPValue")

#correct number of features
unadjustedP_g %>% distinct(feature1) %>% nrow()


library(qvalue)

#calculates the qvalue of each pvalue
#the q-value of a test (pvalue) measures the proportion of false positives incurred 
#(called the false discovery rate) when that particular test is called significant
#https://www.rdocumentation.org/packages/qvalue/versions/2.4.2/topics/qvalue
qvalue(unadjustedP)

#gets the qvalues
qValue <- qvalue(unadjustedP)$qvalue

#makes the rownames (one set of features) into a variable
qValue_g <- qValue %>% mutate(feature1 = row.names(qValue))

#gathers the qvalues into long form
qValue_g <- qValue_g %>% gather(1:235, key = "feature2", value = "qValue")

#correct number of features
qValue_g %>% distinct(feature1) %>% nrow() 

nrow(unadjustedP_g)

#adds qvalues to unadjusted pvalues
unadjustedP_g <- unadjustedP_g %>% left_join(qValue_g, by = c("feature1", "feature2"))

nrow(unadjustedP_g)

#correct number of features
unadjustedP_g %>% distinct(feature1) %>% nrow()

#the qvalues joined correctly to the unadjusted pvalues
unadjustedP_g %>% filter(is.na(unadjustedPValue)) %>% filter(!is.na(qValue))
unadjustedP_g %>% filter(is.na(qValue)) %>% filter(!is.na(unadjustedPValue))

#when I exclude rows with NA pvalues and qvalues, there are 234 feature1s and 234 feature2s
#because 02e4bec9ed20bee6ece6416ddbbadf3c30c7cb64contig_159 is missing from feature1 and 
#ff5d4f76f5d2b4390102b832486c425b662b37ebcontig_808 is missing from feature2
unadjustedP_g %>% distinct(feature1) %>% anti_join(unadjustedP_g %>% filter(!is.na(qValue)), by = c("feature1"))
unadjustedP_g %>% distinct(feature2) %>% anti_join(unadjustedP_g %>% filter(!is.na(qValue)), by = c("feature2"))

#when I exclude rows with NA pvalues and qvalues, feature1 is missing 02e4bec9ed20bee6ece6416ddbbadf3c30c7cb64contig_159
#but feature2 has 02e4bec9ed20bee6ece6416ddbbadf3c30c7cb64contig_159
unadjustedP_g %>% filter(feature1 == "02e4bec9ed20bee6ece6416ddbbadf3c30c7cb64contig_159") %>% filter(!is.na(qValue))
unadjustedP_g %>% filter(feature2 == "02e4bec9ed20bee6ece6416ddbbadf3c30c7cb64contig_159") %>% filter(!is.na(qValue))

#when I exclude rows with NA pvalues and qvalues, feature1 has ff5d4f76f5d2b4390102b832486c425b662b37ebcontig_808 but 
#feature2 is missing ff5d4f76f5d2b4390102b832486c425b662b37ebcontig_808
unadjustedP_g %>% filter(feature1 == "ff5d4f76f5d2b4390102b832486c425b662b37ebcontig_808") %>% filter(!is.na(qValue))
unadjustedP_g %>% filter(feature2 == "ff5d4f76f5d2b4390102b832486c425b662b37ebcontig_808") %>% filter(!is.na(qValue))

unadjustedP_g %>% filter(!is.na(qValue)) %>% distinct(feature1) %>% nrow()
unadjustedP_g %>% filter(!is.na(qValue)) %>% distinct(feature2) %>% nrow()

#excludes the rows with NA qValues and pValues
unadjustedP_g <- unadjustedP_g %>% filter(!is.na(qValue))

str(unadjustedP_g)

#Sean told me to only include correlations with FDR/qvalue <.1 
#when I do this, there are still a large number of correlations and the 
#mean unadjusted pvalue increases which makes sense
unadjustedP_g %>% filter(qValue < .1) %>% summarize(meanPValue = mean(unadjustedPValue), n = n())
unadjustedP_g %>% summarize(meanPValue = mean(unadjustedPValue), n = n())

unadjustedP_g %>% distinct(feature1) %>% nrow()
unadjustedP_g %>% distinct(feature2) %>% nrow()

unadjustedP_g %>% filter(is.na(unadjustedPValue))
unadjustedP_g %>% filter(is.na(qValue))

#makes qvalue NA if qvalue is >= .1
unadjustedP_g <- unadjustedP_g %>% mutate(qValue = ifelse(qValue < .1, qValue, NA))

unadjustedP_g %>% filter(qValue < .1)
unadjustedP_g %>% filter(qValue >= .1)

#gets rid of unadjusted pvalue variable
unadjustedP_g <- unadjustedP_g %>% select(-unadjustedPValue)

##makes a copy of unadjustedP_g but with feature1 and feature2 switched
copy <- unadjustedP_g
copy$feature1New <- copy$feature2
copy$feature2New <- copy$feature1
copy <- copy %>% select(-c(feature1, feature2))

colnames(copy)[2:3] <- c("feature1", "feature2")

copy <- copy %>% select(feature1, feature2, qValue)

#makes a dataframe with the original unadjustedP_g and the unadjustedP_g but 
#with the feature1 and feature2 variables switched
cellsToGet <- bind_rows(unadjustedP_g, copy)

#now feature1 and feature2 each have 235 distinct IDs
cellsToGet %>% distinct(feature1) %>% nrow()
cellsToGet %>% distinct(feature2) %>% nrow()

#makes a dataframe with the correlations between each pair of features
corDF <- test$r %>% as.data.frame()

#makes a variable for one set of the features
corDF <- corDF %>% mutate(feature1 = row.names(corDF))

#gathers the correlations into long form
corDF <- corDF %>% gather(1:235, key = "feature2", value = "correlation")

#there are 235 distinct feature1s and feature2 which is correct
corDF %>% distinct(feature1) %>% nrow()
corDF %>% distinct(feature2) %>% nrow()

#the only pairs of features that are in corDF but not cellsToGet are for pairs in which
#a feature is compared to itself
corDF %>% anti_join(cellsToGet, by = c("feature1", "feature2"))
corDF %>% anti_join(cellsToGet, by = c("feature1", "feature2")) %>% filter(feature1 != feature2)
cellsToGet %>% anti_join(corDF, by = c("feature1", "feature2")) 

nrow(corDF)

#adds cellsToGet (which has the info on whether a correlation was significant) to corDF
corDF <- corDF %>% left_join(cellsToGet, by = c("feature1", "feature2"))

nrow(corDF)

#there are no correlations for which the qvalue is >= .1 because I made these qvalues NAf
corDF %>% filter(qValue >= .1)

#make the correlation value NA if the qvalue is NA 
#this makes nonsignificant correlations and correlations for features compared to themselves NA
corDF <- corDF %>% mutate(correlation = ifelse(is.na(qValue), NA, correlation))

corDF %>% filter(is.na(correlation)) %>% filter(!is.na(qValue))
corDF %>% filter(is.na(qValue)) %>% filter(!is.na(correlation))

#str(corDF)

#makes correlation results into a dataframe
#corDF <- corDF %>% as.data.frame()

#makes a variable for one feature in the pair
#corDF <- corDF %>% mutate(feature1 = rownames(corDF))

#gather correlation dataframe
#corDF <- corDF %>% gather(1:235, key = "feature2", value = "correlation")

###checking correlation

#norm[1:5, 1:4]

#head(corDF)

#cor(norm[,1], norm[,2], method = "spearman")

#cor(norm[,2], norm[,3], method = "spearman")

#corDF %>% filter(feature1 == "031a23c740bcf047c85fd29fb9b75822ea473a09contig_82") %>% filter(feature2 == "0468fcf94fc57999b4d6518905ea6408230ae54fcontig_404")

#gets rid of qvalue variable because I am not going to plot it
corDF <- corDF %>% select(-qValue)

#spreads correlation data
corDFS <- corDF %>% spread(key = feature2, value = correlation)

#makes rownames of spread correlation data the feature1 variable
rownames(corDFS) <- corDFS$feature1
corDFS <- corDFS %>% select(-feature1)

#correct dimensions
dim(corDFS)

#I can't cluster the features when they have NA values
#hc <- hclust(dist(corDFS), "ward.D2")

#so I replaced the NAs with 0
#https://support.bioconductor.org/p/45575/
corDFS_forClustering <- corDFS
corDFS_forClustering[is.na(corDFS_forClustering)] <- 0

#clusters the features based on their correlations after replacing 
#the NAs with 0
#from https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/hclust
hc <- hclust(dist(corDFS_forClustering), "ward.D2")

png("correlationBasedOnHeterotrophAbundances_noContamination_noMock_noProSyn_withoutSample1223_featureClustering.png")
#plots the clustering
plot(hc, labels = FALSE)
dev.off()

plot(hc)

hc$labels

hc$order

hc$labels[hc$order]

corDF %>% distinct(feature1) %>% nrow()
corDF %>% distinct(feature2) %>% nrow()

#puts the features in the order that they appear on the dendrogram
corDF$feature1 <- factor(corDF$feature1, levels = hc$labels[hc$order])
corDF$feature2 <- factor(corDF$feature2, levels = hc$labels[hc$order])

levels(corDF$feature1)

corDF %>% arrange(feature1) %>% distinct(feature1)

#plots correlations for each pair of features, leaving the nonsignificant correlations 
#and correlations between a feature compared to itself NA
corDF %>% ggplot(aes(x = feature1, y = feature2, fill = correlation)) + geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  labs(fill = "Significant Spearman correlation across samples") + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.title=element_text(size=6), legend.text=element_text(size=6)) + 
  theme(legend.key.size = unit(.4, "cm"))

ggsave("correlationBasedOnHeterotrophAbundances_noContamination_noMock_noProSyn_withoutSample1223.png", dpi = 300, height = 8, width = 8)

corDF %>% ggplot(aes(x = feature1, y = feature2, fill = correlation)) + geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  labs(fill = "Significant Spearman correlation across samples") + 
  theme(legend.title=element_text(size=6), legend.text=element_text(size=6)) + 
  theme(legend.key.size = unit(.4, "cm")) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

str(corDF)

#corDF %>% arrange(desc(feature1)) %>% distinct(feature1, shortTaxa.x) %>% View()
#corDF %>% arrange(feature1) %>% distinct(feature1, shortTaxa.x) %>% View()

#loads cleaned up taxonomy of the same length features, no pro, no syn, no contamination, no mock samples, no 1223 or JW3
#this is from makeTreeNodes_noContamination_noMock_noProSyn_without1223.R
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv")

tax %>% distinct(`Feature ID`) %>% nrow()
tax %>% nrow()

#makes a variable for the combined different feature IDs
tax <- tax %>% mutate(feature = str_c(`Feature ID`, V1))

#gets rid of unnecessary variables
tax <- tax %>% select(feature, shortTaxa)

nrow(corDF)

#adds taxonomy to the correlation data
#this messes up the ordering of the features
corDF <- corDF %>% left_join(tax, by = c("feature1" = "feature")) %>% left_join(tax, by = c("feature2" = "feature"))

nrow(corDF)

corDF %>% filter(is.na(shortTaxa.x))

#puts the features in the order that they appear on the dendrogram
corDF$feature1 <- factor(corDF$feature1, levels = hc$labels[hc$order])
corDF$feature2 <- factor(corDF$feature2, levels = hc$labels[hc$order])

str(corDF)
corDF %>% arrange(feature1) %>% distinct(feature1)

#writes the taxonomies in the correct order to text files
#writes in the forward and reverse order because this is necessary 
#to label both axes

corDF %>% arrange(feature1) %>% distinct(feature1, shortTaxa.x) %>% select(shortTaxa.x) %>% 
  write.table("treeOrderTaxonomyForCorrelationHeatmap.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

corDF %>% arrange(desc(feature1)) %>% distinct(feature1, shortTaxa.x) %>% select(shortTaxa.x) %>% 
  write.table("treeOrderTaxonomyForCorrelationHeatmap_reverse.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


#from getColorCodingForPhylogeneticTree.R
color <- read_csv("colorCodingPhylogeneticTree_notInOrder.csv")

#gets rid of unnecessary variables
color <- color %>% select(shortTaxa, class, featureContig)

#puts the features in the order that they appear on the dendrogram
color$featureContig <- factor(color$featureContig, levels = levels(corDF$feature1))

color %>% arrange(featureContig) %>% distinct(featureContig)
color %>% arrange(desc(featureContig)) %>% distinct(featureContig)

color %>% arrange(featureContig) %>% distinct(featureContig, shortTaxa, class) %>% select(featureContig, shortTaxa, class)
color %>% arrange(desc(featureContig)) %>% distinct(featureContig, shortTaxa, class) %>% select(featureContig, shortTaxa, class)

color %>% arrange(featureContig) %>% distinct(featureContig, shortTaxa, class) %>% select(shortTaxa, class) %>% 
  write.table("treeOrderTaxonomyForCorrelationHeatmap_colorCoding.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

color %>% arrange(desc(featureContig)) %>% distinct(featureContig, shortTaxa, class) %>% select(shortTaxa, class) %>% 
  write.table("treeOrderTaxonomyForCorrelationHeatmap_colorCoding_reverse.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


#contig_7 is the the sequence that matches to Alt Mac MIT1002 (the lab strain)
read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv") %>% filter(V1 == "contig_7")

#"f1e1d1742851920183773cc5710788dbe75225e8contig_7" is the 155th entry
levels(corDF$feature1)[155]

corDF %>% arrange(feature1) %>% distinct(feature1, shortTaxa.x)
corDF %>% arrange(feature1) %>% distinct(feature1, shortTaxa.x) %>% slice(155:160)



