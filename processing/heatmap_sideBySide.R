library(tidyverse)

setwd("~/Dropbox (MIT)/clusterAllisonSeanIllumina16S/")

#this came from clusterSeanAllisonIllumina16SSeqs.txt (checked)
biom <- read.table("feature-table-fromBiom.tsv") 

#these are feature/sample names
#these came from feature-table-fromBiom.tsv
names <- "V1	8349da6d92720a4837cd8527f0423894	contig3	contig37	contig7	10933022f9b86ef255f1b5b010b05854	contig53	17d786ce8f68c39745f6ba3530db951e	contig173	19fcd1e5141188d040a4efe2da213f8f	contig57	230a096a7635854364026cb6ca02206a	contig111	2e9599d4d171842cde3dfc750417c84c	contig213	3266fefce67a9d030cfce9f907203300	contig70	5198bb2290cfb9665678ecca1522e24f	contig334	63586b50310b55083318fa8d1ab1b37c	contig165	6a5409e7695ca0fef67c8846b4e8ea6d	contig149	7dc54e9f3b846ff5b9afa2e14b7f286a	contig172	8334bd9ea83ac028f4e11278852ba47f	contig11	83ae3a40a367793ac7aaf9c79d69df6a	contig305	981cec2043503f95e40975ffb4cfa4da	contig164	98b0b4fceed91f45a392fac33e8e5660	contig25	9d4acd652505aeb73c49c8f312279c74	contig12	9dd7c31f5138dffa8e77f2ebfb64b949	contig163	a3fcbae413e192f52aacfc6db72144af	contig261	b8cd19f852be08768528b4cc05f25652	contig108	bbe5a62f75965b6575b4ef1c5de9e73c	contig177	c00cad92f1901b38bdda414be9557ffa	contig160	c405ad29c5aed7e6eb4a194bce9f095a	contig459	c5b7f441b9603531adf4612503363cd0	contig76	d360005c2c6d97fc01a5894838840162	contig8	dfd81eaecbffe0d522375fcf8ee36a95	contig36	f4eea3fb91daf0103a84cec6bb25596a	contig109	f503c181a3ba69fddc62a7dd67b25af5	contig400	f574bfc1c168d5c704760ed661715512	contig61	07c8b17c30d6bced44e26dec897564f5	0997194aed46dc2a2ba3e386654157d8	618055a769277c14994a90850183555f	71cd617896f5096431158ecea9c01f47	b0071d99acd2fc62de30895fcc0b5ccc	contig444	efaa00ade414bb47f67daa3ae3fbd7d3"


#makes sample names into list
namesNew <- str_split(names, pattern = "\t")
namesNew <- namesNew[[1]]

#same length which is good
ncol(biom)
length(namesNew)

#makes the column names the sample names
colnames(biom) <- namesNew

#there are 35 clustered features
nrow(biom)

#there are 65 sample columns
ncol(biom)


#gathers into long format
biomG <- biom %>% gather(2:66, key = "sample", value = "present")

#there are 35 clustered features and 65 samples
biomG %>% distinct(V1) %>% nrow()
biomG %>% distinct(sample) %>% nrow()

biomG %>% filter(present > 0) %>% distinct(V1) %>% nrow()
biomG %>% filter(present > 0) %>% distinct(sample) %>% nrow()

#there are 65 times a sequence is found in a sample 
#because this is just presence/absence because I imported features to qiime2
#I need to connect the features to their abundances in each sample
biomG %>% filter(present > 0)

#old <- read.table("feature-table-fromBiom.tsv")

#oldG <- old %>% gather(2:60, key = "sample", value = "present")

#this is the same as the nonclustered results which is good!
#oldG %>% filter(present > 0)

#loads taxonomy of clustered features
taxa <- read.table("taxonomySeanAllisonFeatures.tsv", fill = TRUE)

#replaces NAs with empty strings so that I can join all of the taxonomy variables together
taxa[is.na(taxa)]
taxa[is.na(taxa)] <- ""

#combines all of the taxonomy variables together into one variable
taxa <- taxa %>% mutate(taxon = str_c(V2, V3, V4, V5, V6, V7, V8, V9, sep = " "))

#gets rid of unnecessary taxonomy variables
taxa <- taxa %>% select(V1, taxon)

#gets rid of row with variable names
taxa <- taxa[-1,]

#there are 35 features
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



###sean

#from heatmap_noContamination_noLowAbundance_noMock_noProSyn_without1223.R (checked) in Sean16SSep2019 directory
seqGSean <- read_csv("~/Dropbox (MIT)/Sean16SSep2019/seqG_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.csv")

seqGSean %>% distinct(V1) %>% nrow()
seqGSean %>% distinct(sample) %>% nrow()

#gets the cultures/samples that Allison sampled for Illumina 16S
seqGSean <- seqGSean %>% filter(sample %in% c("MED4", "9313", "NATL2A"))

seqGSean %>% distinct(sample)

#gets only sequences that are present in samples
seqGSean <- seqGSean %>% filter(abundance > 0)

seqGSean %>% distinct(V1) %>% nrow()
seqGSean %>% nrow()

#these are the features that have multiple sequences
seqGSean %>% group_by(V1) %>% summarize(n = n()) %>% filter(n > 1)

#makes variable for data source
seqGSean <- seqGSean %>% mutate(source = "sean")




###allison 

#from abundanceHeatmap.R (checked) in analysisOfAllison16S directory
seqGAllison <- read_csv("~/Dropbox (MIT)/analysisOfAllison16S/seqG.csv")

seqGAllisonSummary <- seqGAllison %>% group_by(sample) %>% summarize(total = sum(`relative abundance`))

nrow(seqGAllison)

#adds total number of Pro reads in each sample to seqGAllison
seqGAllison <- seqGAllison %>% left_join(seqGAllisonSummary, by = c("sample"))

nrow(seqGAllison)

#there are no NA proSample values in seqNamesG
seqGAllison %>% filter(is.na(total))

#makes variable for the number of reads divided by the total number of reads in the sample
seqGAllison <- seqGAllison %>% mutate(normAbun = `relative abundance`/total) 

#excludes sequences that are not in a sample
seqGAllison <- seqGAllison %>% filter(`relative abundance` > 0)

#seqGAllison %>% filter(normAbun < 0.002)

seqGAllison %>% filter(str_detect(taxa, "Prochlorococcus")) %>% distinct(seq)

#excludes Pro sequences
seqGAllison <- seqGAllison %>% filter(!str_detect(taxa, "Prochlorococcus"))


seqGAllison %>% distinct(seq) %>% nrow()
seqGAllison %>% nrow()

#excludes taxa variable because I am going to add the combined Sean and Allison
#taxonomy dataset
seqGAllison <- seqGAllison %>% select(-taxa)

#makes a variable for the data source
seqGAllison <- seqGAllison %>% mutate(source = "allison")


###both 

#makes the seqG dataset column names match
colnames(seqGSean)[2] <- "seq"
colnames(seqGAllison)[3] <- "abundance"

#replaces 9313 with MIT9313 in sample IDs for Sean seqG 
seqGSean %>% distinct(sample)
seqGSean$sample <- str_replace(seqGSean$sample, "9313", "MIT9313")

seqGSean %>% distinct(sample)
seqGAllison %>% distinct(sample)

#combines seqGSean and seqGAllison
seqG <- bind_rows(seqGSean, seqGAllison %>% select(-c(total, normAbun)))


#none of the seq IDs are used in both the Allison and Sean analysis
seqG %>% group_by(seq) %>% distinct(source) %>% summarize(n = n()) %>% filter(n > 1)


#all of Allison feature names match across biomG and the abundance data
#and none of Sean feature names match
biomG %>% semi_join(seqG, by = c("sample" = "seq"))
seqG %>% semi_join(biomG, by = c("seq" = "sample"))
seqG %>% semi_join(biomG, by = c("seq" = "sample")) %>% distinct(source)
seqGAllison %>% anti_join(biomG, by = c("seq" = "sample"))


seqG %>% mutate(contig = ifelse(source == "sean", str_extract(seq, "contig.*"), seq)) %>% 
  group_by(contig) %>% distinct(seq) %>% summarize(n = n()) %>% filter(n > 1)

#fixes sean seq IDs so that they match biomG
seqG$seq <- ifelse(seqG$source == "sean", str_extract(seqG$seq, "contig.*"), seqG$seq)
seqG$seq <- ifelse(seqG$source == "sean", str_replace(seqG$seq, "_", ""), seqG$seq)

#now all of the sequence IDs in seqG and biomG match
biomG %>% anti_join(seqG, by = c("sample" = "seq"))
seqG %>% anti_join(biomG, by = c("seq" = "sample"))


nrow(biomG)

#adds abundances to biomG
biomG <- biomG %>% left_join(seqG, by = c("sample" = "seq"))

nrow(biomG)

biomG %>% distinct(present)

#gets rid sequences that aren't in samples
biomG <- biomG %>% filter(present == 1)

#these are the features that have multiple sequences/features
biomG %>% group_by(V1) %>% summarize(n = n()) %>% arrange(desc(n))

biomG %>% distinct(sample.y)

head(biomG)

#gets sum of abundances for features in samples
summary <- biomG %>% group_by(V1, source, sample.y, taxon, shortTaxa) %>% summarize(sumAbundance = sum(abundance))

head(summary)

#these are the features that are found in multiple cultures and/or in both Allison and Sean's samples
summary %>% group_by(V1) %>% summarize(n = n()) %>% arrange(desc(n)) %>% filter(n > 1)

#there are 35 features
summary %>% ungroup() %>% distinct(V1) %>% nrow()

#spreads data by whether from Allison or Sean
#summarySpread <- summary %>% spread(key = "source", value = "sumAbundance", fill = 0)

#these are the features that are found in different types of cultures
#summarySpread %>% group_by(V1) %>% summarize(n = n()) %>% filter(n > 1)

#makes a variable for the proportion of the reads of each feature that are in each culture 
#that are from Sean's experiment
#summarySpread <- summarySpread %>% mutate(propSean = sean/sum(sean, allison))

#gets the total number of reads of each feature that are in each culture across Allison 
#and Sean's experiment
#summarySpread <- summarySpread %>% mutate(totalAllisonSean = sum(sean, allison))

#makes a variable for whether the feature in a culture is in both Allison and Sean's experiment
#summarySpread <- summarySpread %>% mutate(foundInSeanAllison = ifelse((sean > 0 & allison > 0), TRUE, FALSE))

#summarySpread <- summarySpread %>% select(-c(sean, allison))


#from abundanceHeatmap.R (checked) in analysisOfAllison16S
proAllison <- read_csv("~/Dropbox (MIT)/analysisOfAllison16S/seqG.csv")

#gets the Allison Pro. features
proAllison <- proAllison %>% filter(str_detect(taxa, "Prochlorococcus"))

#gets the Allison Pro. features that are in samples
proAllison <- proAllison %>% filter(`relative abundance` > 0)

#renames abundance column
#colnames(proAllison)[4] <- "proAbundanceAllison"

#gets just the number of Pro reads in each sample
proAllison <- proAllison %>% select(3:4)


#from getSequencesForTree_noContamination_noPro.R (checked) in Sean16SSep2019 directory
proSean <- read_csv("~/Dropbox (MIT)/Sean16SSep2019/seqG_noContamination_withProAndSyn.csv")

proSean %>% distinct(Feature.ID) %>% nrow()

#gets the sequences from the samples that were also in Allison's experiment
proSean <- proSean %>% filter(sample %in% c("MED4", "9313", "NATL2A"))

proSean %>% distinct(sample)

#gets just the Pro sequences
proSean <- proSean %>% filter(shortTaxa == "Prochlorococcus")

#replaces 9313 sample ID with MIT9313 
proSean <- proSean %>% mutate(sample = ifelse(sample == "9313", "MIT9313", sample))

proSean %>% distinct(sample)

#gets the number of Pro reads in each sample
proSean <- proSean %>% select(sample, abundance)

#renames Pro read abundance variable
#colnames(proSean)[2] <- "proAbundanceSean"

#combines Sean and Allison Pro read abundance data sets
#pro <- proSean %>% left_join(proAllison, by = c("sample"))

#makes variable for source of Pro abundance
proAllison <- proAllison %>% mutate(source = "allison")
proSean <- proSean %>% mutate(source = "sean")

#renames variables
colnames(proSean)[2] <- "proAbundance"
colnames(proAllison)[2] <- "proAbundance"

#binds pro abundance datasets together
pro <- bind_rows(proAllison, proSean)


#gets the total number of Pro reads in each culture across Sean and Allison's experiments
#summarySpread <- summarySpread %>% mutate(sumPro = sum(proAbundanceSean, proAbundanceAllison))

#summarySpread %>% filter(totalAllisonSean == 0)

#there are 35 features which is correct
summary %>% ungroup() %>% distinct(V1) %>% nrow()

#summaryS <- summary %>% spread(key = source, value = sumAbundance, fill = 0)

#summary <- summaryS %>% gather(allison:sean, key = "source", value = "sumAbundance")

#there are 35 features which is correct
summary %>% ungroup() %>% distinct(V1) %>% nrow()

#adds 1 to the relative abundances of 0 so that I can take the log transformation
#summary$sumAbundance <- ifelse(summary$sumAbundance == 0, 
                               #summary$sumAbundance + 1, summary$sumAbundance)

nrow(summary)

#adds Sean and Allison Pro read abundance to summary by sample and source
summary <- summary %>% left_join(pro, by = c("sample.y" = "sample", "source"))

nrow(summary)

summary %>% filter(is.na(proAbundance))


#makes a variable for the log of the total number of reads of each feature 
#divided by the number of Pro in each sample by source
#takes log10 because this is what I did for Allison seqs previously and Sean said to take 
#log base 10 or base 2
summary <- summary %>% mutate(logNormTotalByPro = log10(sumAbundance/proAbundance))

summary %>% filter(is.na(logNormTotalByPro))

#there are 35 features which is correct
summary %>% ungroup() %>% distinct(V1) %>% nrow()

#summarySpread$foundInSeanAllison <- factor(summarySpread$foundInSeanAllison, levels = c(TRUE, FALSE))


order <- rev(c("48d7e3444f85af8c2bd1f10624778437d82e583d",
               "57872c999308bfec89b949c43c4a4b51aafdb0f3",
               "0b51c0f2d5930c46f268b8592b9fca23fe00f8c1",
               "ca0cace32a7c6c5813227f0a06c9f227ca845671",
               "9773ae42f9c1076c7bb03948d9dbc4ab89461f62",
               "47d3b61a3ec945aecde3475e6d98986b35e8c6f3",
               "fe66dcc0a86631393dce511d9e2348888ef42c46",
               "53e028c6f56c4210d9312ed07831258ad003a9a2",
               "39c7cac93b7e1913f21cbcf959ec738710dc3381",
               "5c8141fea9c85ab22abad95fbf8c1b66873d18de",
               "c4f4e9f84d490e39681ebe5ee2b11ad0a3db4dfd",
               "5715644b7f4888973e7f8c1456391a084587d576",
               "07d0bb5d023b91f79a21fa1f08de76da9475afa7",
               "b92cf2c1956ebbe6b0e7a4b7edd646d4dd84f9dc",
               "3033e7739060a50de88c2b06168e26d681fca9c0",
               "4d42d41088f08ee2db9f70375560715a507c91c5",
               "8785e4a7f77eacd1d9863a2e9db9099f7e2314bf",
               "d089cd89ef62f29a2a3399c9fe6344baa260aa45",
               "ebf4dff2e8291d0807859a8289a864d05d52698c",
               "16e13db71aa93d0c849a0a4c8a5c14104ef149ea",
               "6ee71749cdfca7c10ec389ae08d312ecbe37d0aa",
               "7176f58de210a7cea377d555b1d247c67b69855f",
               "7f131455a0d7d583a795db8bf875920e1138077f",
               "b1d24a3d6bc8f67e2e20782a80fc9a1981a6506b",
               "f85b9acee94c13d2d091a47ba52907ad2327947b",
               "1acb4b1dcfdf5be27c536949eb15515631711d41",
               "45ac47fbf72fbdad0846bcf5f5b804f65f32fd54",
               "5865070a00b3df85a584690732ff897c19063217",
               "27fedabca0a1e1f5d8d3b279d42c25113299de01",
               "8128f23d8559bee53a0a97f540b4b47c77775982",
               "4a54d389b1a158f478d861637299bf6ae04aae80",
               "b28441d052be31e9e23600490c26e647125695f3",
               "6b28f251bd68700bcbe4a7e4986b19d7c825bcf6",
               "9b5f9ebe90d1eaaa37237e4a51152aa79c966e84",
               "8391c3f2f297927f0197174d956e327bd6ae696a"))


summary$V1 <- factor(summary$V1, levels = order)

#some of the squares which do not have any of the respective sequence were showing up as 
#light pink instead of white so this makes so they will be white
#summary <- summary %>% mutate(logNormTotalByPro = ifelse(sumAbundance == 1, -4.84761584, logNormTotalByPro))

summary %>% ggplot(aes(sample.y, V1)) + geom_tile(aes(fill = logNormTotalByPro)) + facet_wrap(~source) +
  labs(x = "", fill = "", y = "") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_gradient(low = "white", high = "red") + scale_color_manual(values = c("red", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.ticks.x = element_blank(), 
        axis.text=element_text(size=12, color = "black"))

ggsave("heatmap_sideBySide.png", dpi = 600, height = 6, width = 8)

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

ggsave("correlationHeatmap_combinedSeanAllisonIllumina.png", dpi = 600, height = 6, width = 8)
ggsave("correlationHeatmap_combinedSeanAllisonIllumina.svg", dpi = 600, height = 6, width = 8)







