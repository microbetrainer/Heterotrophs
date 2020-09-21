library(tidyverse)
library(readr)
library(readxl)
library(picante)
library(phyloseq)


setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads distances between nodes (hetetrophs) on phylogenetic tree
#this came from heatmap_noContamination_noLowAbundance_noMock_noProSyn_without1223.R
load("distancesBtwPhylogeneticTreeNodes_noMock_noProSyn_withoutSample1223.Rda")

dist <- distancesBtwNodes

#makes into dataframe
dist <- dist %>% as.data.frame()

#makes a variable for the rownames (one of the features)
dist$feature1 <- rownames(dist)

colnames(dist)[235:236]

#gathers into long form
dist <- dist %>% gather(1:235, key = "feature2", value = "dist")

#loads abundances of each feature in each sample
#this came from heatmap_noContamination_noLowAbundance_noMock_noProSyn_without1223.R
seqG <- read_csv("seqG_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.csv")

#gets just the sequences that are present in a sample
seqG <- seqG %>% filter(abundance > 0)

#correct numher of features
seqG %>% distinct(V1) %>% nrow()

#gets rid of unnecessary variable
seqG <- seqG %>% select(-abundance)

#correct number of features
dist %>% distinct(feature1) %>% nrow()
dist %>% distinct(feature2) %>% nrow()

#gets list of samples
samples <- seqG %>% distinct(sample)
samples <- samples$sample %>% as.list()

#correct number of samples
length(samples)

dfList <- list()
num = 1

#iterates through list of samples
for (i in samples) {
  print(i)
  #gets the distint features that are present in the sample
  features <- seqG %>% filter(sample == i) %>% distinct(V1)
  #makes a dataframe that only has the distances betweeen pairs of features, for which both features are present 
  #in the sample
  df <- dist %>% semi_join(features, by = c("feature1" = "V1")) %>% semi_join(features, by = c("feature2" = "V1"))
  #makes a variable for the sample ID
  df <- df %>% mutate(sample = i)
  #adds df to the list of dataframes
  dfList[[num]] <- df
  print(num)
  #adds 1 to the number so the next dataframe will get appended to the end of the list
  num <- num + 1
}

#binds the list of dataframes into one dataframe
dfComb <- bind_rows(dfList)

#correct number of samples
dfComb %>% distinct(sample) %>% nrow()

#excludes rows for which the two rows are identical
dfComb <- dfComb %>% filter(feature1 != feature2)

dfComb %>% nrow()

#this gets rid of the duplicate rows for each pairwise comparison
#so there are half the rows when you do this
dfComb %>% filter(feature1 < feature2) %>% nrow()

5061 * 2


#1213 is 1313
#WH7803 is WH7801
#SYN1320 is SYN1220
#9201.2 is 9311

#now there are 73 samples instead of 74 because
#1213 is no longer in dfComb because it only has one feature so it was only 
#compared to itself
dfComb %>% distinct(sample) %>% nrow()
dfComb %>% filter(sample == "1213")

dfComb %>% filter(sample == "WH7803")
dfComb %>% filter(sample == "SYN1320")
dfComb %>% filter(sample == "9201.2")
dfComb %>% filter(sample == "9201.1")

#fixes messed up sample IDs
#dfComb$sample <- ifelse(dfComb$sample == "WH7803", "WH7801", dfComb$sample)
dfComb$sample <- ifelse(dfComb$sample == "SYN1320", "SYN1220", dfComb$sample)
dfComb$sample <- ifelse(dfComb$sample == "9201.2", "9311", dfComb$sample)
dfComb$sample <- ifelse(dfComb$sample == "9201.1", "9201", dfComb$sample)


#ASN9601 corresponds to AS9601
dfComb$sample <- ifelse(dfComb$sample == "ASN9601", "AS9601", dfComb$sample)

#S9503 corresponds to SYN9503
dfComb$sample <- ifelse(dfComb$sample == "S9503", "SYN9503", dfComb$sample)

#seqComb$sample <- ifelse(seqComb$sample == "8102", "SYNCLADEXVI", seqComb$sample)
dfComb$sample <- ifelse(dfComb$sample == "8102", "SYNCLADEXVI", dfComb$sample)



#gets the samples ordered by descending median pairwise distance
#calculates the median pairwise distance for each sample
order <- dfComb %>% filter(feature1 < feature2) %>% group_by(sample) %>% summarize(med = median(dist))
order <- order %>% arrange(desc(med))
order <- order %>% distinct(sample)
order <- order$sample %>% as.list()

#correct number of samples
length(order)

#puts samples in order of median pairwise distance
dfComb$sample <- factor(dfComb$sample, levels = order)

dfComb %>% filter(feature1 < feature2) %>% ggplot(aes(x = sample, y = dist)) + geom_point(color = 'grey') + 
  geom_boxplot(outlier.shape = NA, fill = NA) +
  labs(y = "Pairwise phylogenetic distance", x = "Sample") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("phylogeneticDistancesBtwHeterotrophs_noContamination_noMock_noProSyn_withoutSample1223.png", dpi = 600, height = 10, width = 13)
ggsave("phylogeneticDistancesBtwHeterotrophs_noContamination_noMock_noProSyn_withoutSample1223.svg", dpi = 600, height = 10, width = 13)














