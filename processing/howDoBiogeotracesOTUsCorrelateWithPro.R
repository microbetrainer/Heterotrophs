#for the top 10 OTUs in cultures can you check how their abundance 
#varies with the abundance of things called Prochlorococcus? 
#basically just calculate/plot correlations 

library(tidyverse)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

read_csv("mostUbiquitiousNonProNonSynOTUsInSeanCultures_seanAndBiogeotracesClusteredTogether.csv") %>% nrow()
read_csv("occupancyAbundanceNonProNonSynOTUsInSeanCultures_seanAndBiogeotracesClusteredTogether.csv") %>% nrow()

read_csv("mostUbiquitiousNonProNonSynOTUsInSeanCultures_seanAndBiogeotracesClusteredTogether.csv") %>% 
  anti_join(read_csv("occupancyAbundanceNonProNonSynOTUsInSeanCultures_seanAndBiogeotracesClusteredTogether.csv"), by = c("OTU_ID"))

#loads the most ubiquitious non-pro, non-syn OTUs in sean samples 
#this is based on the sean and biogeotraces sequences being clustered together
#from ubiquityAnalysisSeanBiogeotraces_clustered.R
ubi <- read_csv("mostUbiquitiousNonProNonSynOTUsInSeanCultures_seanAndBiogeotracesClusteredTogether.csv")

ubi %>% head()
nrow(ubi)

#loads the sean and biogeotraces sequences present in each cluster
#from clusterSeanAndBiogeotracesSeqs.txt (checked)
biomG <- read_csv("combinedBiogeotracesBiom.csv")

#28022 clusters 
#1184489 sequences
head(biomG)
biomG %>% distinct(OTU_ID) %>% nrow()
biomG %>% distinct(sample) %>% nrow()
biomG %>% nrow()

str(biomG)

biomG %>% distinct(present)

#each sequence is only in one cluster
biomG %>% group_by(sample) %>% summarize(n = n()) %>% filter(n != 1)

#renames biomG column names
head(biomG)
colnames(biomG)
colnames(biomG)[2]
colnames(biomG)[2] <- "biogeotracesSequence" 
head(biomG)

#extracts the biogeotraces sample from the biogeotraces sequence ID
biomG <- biomG %>% mutate(biogeotracesSample = str_extract(biogeotracesSequence, "[a-zA-Z0-9]*"))
head(biomG)

biomG %>% filter(is.na(biogeotracesSample))

#gets rid of unnecessary variable
biomG <- biomG %>% select(-present)
head(biomG)

#excludes sean sequences
biomG %>% filter(!str_detect(biogeotracesSequence, "SRR"))
biomG <- biomG %>% filter(str_detect(biogeotracesSequence, "SRR"))

biomG %>% distinct(biogeotracesSample) %>% nrow()

#calculate total number of sequences in each biogeotraces samples
#including Pro, Syn sequences
head(biomG)
numSeqsBioGeoSample <- biomG %>% group_by(biogeotracesSample) %>% distinct(biogeotracesSequence) %>% 
  summarize(totalSeqsBioGeo = n())
nrow(numSeqsBioGeoSample)
head(numSeqsBioGeoSample)
str(numSeqsBioGeoSample)

biomG %>% distinct(OTU_ID) %>% nrow()
biomG %>% distinct(biogeotracesSequence) %>% nrow()

head(biomG)
head(ubi)

nrow(ubi)
ubi %>% semi_join(biomG, by = c("OTU_ID")) %>% nrow()
biomG %>% semi_join(ubi, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

#gets the rows in biomG that have the OTUs that are ubiquitious across 
#sean culture samples
#biomG_ubi and biomG just have the biogeotraces samples not the sean samples
biomG_ubi <- biomG %>% semi_join(ubi, by = c("OTU_ID"))

head(biomG_ubi)
biomG_ubi %>% distinct(OTU_ID) %>% nrow()

#loads taxonomies of OTUs
#based on clustering of sean and biogeotraces sequences
#from clusterSeanAndBiogeotracesSeqs.txt (checked)
tax <- read_tsv("sequencesForTree_noContamination_noMock_withProSyn_no1223_biogeotracesTaxonomy.tsv")

#gets rid of first row because it just has comments
tax[1,]
tax <- tax[-1,]

head(tax)

#gets rid of unnecessary variable
tax <- tax %>% select(1,2)

nrow(tax)

##add taxonomy of OTUs to biomG
#biomG has only the biogeotraces samples
#biomG was based on the clustering of sean and biogeotraces sequences
biomG %>% head()
biomG %>% filter(!str_detect(biogeotracesSequence, "SRR"))
biomG %>% distinct(OTU_ID) %>% nrow()

nrow(biomG)
biomG <- biomG %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(biomG)

biomG %>% filter(is.na(Taxon)) %>% nrow()

head(biomG)

biomG %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)
biomG %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(Taxon)

#gets just Pro and Syn OTUs
biomG_proSyn <- biomG %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc"))

biomG_proSyn %>% distinct(Taxon)

#gets the number of biogeotraces sequences of each ubiquitious OTU in each biogeotraces sample
biomG_ubi %>% head()
biomG_ubi <- biomG_ubi %>% group_by(OTU_ID, biogeotracesSample) %>% distinct(biogeotracesSequence) %>% 
  summarize(numSeqOfUbiOTU = n())
biomG_ubi %>% head()

biomG_ubi <- biomG_ubi %>% ungroup()

biomG_ubi %>% head()
biomG_ubi %>% nrow()
biomG_ubi %>% distinct(biogeotracesSample) %>% nrow()
biomG_ubi %>% distinct(OTU_ID) %>% nrow()

#gets the number of Pro, Syn sequences in each biogeotraces sample
biomG_proSyn
biomG_proSyn %>% head()
biomG_proSyn <- biomG_proSyn %>% group_by(biogeotracesSample) %>% distinct(biogeotracesSequence) %>% 
  summarize(numProSynSeq = n())
biomG_proSyn %>% head()

#ungroup dataframes
biomG_ubi <- biomG_ubi %>% ungroup()
biomG_proSyn <- biomG_proSyn %>% ungroup()

#not all of the biogeotraces samples have Pro, Syn sequences
biomG_proSyn %>% head()
biomG_proSyn %>% distinct(biogeotracesSample) %>% nrow()

#no sean sequences are in biomG
biomG %>% distinct(biogeotracesSample) %>% nrow()
biomG_proSyn %>% distinct(biogeotracesSample) %>% nrow()

#these are the biogeotraces samples that do not have Pro, Syn sequences
biomG %>% anti_join(biomG_proSyn, by = c("biogeotracesSample")) %>% distinct(biogeotracesSample)
biomG %>% anti_join(biomG_proSyn, by = c("biogeotracesSample")) %>% distinct(biogeotracesSample) %>%
  filter(!str_detect(biogeotracesSample, "SRR"))


#gets the biogeotraces samples that are missing from biomG_proSyn
#missingBiogeo <- biomG %>% anti_join(biomG_proSyn, by = c("biogeotracesSample")) %>% distinct(biogeotracesSample) 

#biomG_proSyn %>% head()
#missingBiogeo %>% head()

#makes numProSynSeq value 0 for the biogeotraces samples that are missing from biomG_proSyn
#missingBiogeo <- missingBiogeo %>% mutate(numProSynSeq = 0)
#missingBiogeo %>% head()

#adds missing biogeotraces samples to biomG_proSyn
#nrow(biomG_proSyn)
#biomG_proSyn <- bind_rows(biomG_proSyn, missingBiogeo)
#nrow(biomG_proSyn)

#now all of the biogeotraces samples are in biomG_proSyn
biomG_proSyn %>% head()
biomG_proSyn %>% nrow()
biomG_proSyn %>% distinct(biogeotracesSample) %>% nrow()

#these are the 10 ubiquitious OTUs
biomG_ubi %>% distinct(OTU_ID)

#makes a dataframe for each OTU in biomG_ubi
#ubi_02ba11823bfe93a0250096d95e07e745e40736bc <- biomG_ubi %>% filter(OTU_ID == "02ba11823bfe93a0250096d95e07e745e40736bc")
#ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 <- biomG_ubi %>% filter(OTU_ID == "08d46233bb26721e2336eb640c2dca0dd94497d3")
#ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 <- biomG_ubi %>% filter(OTU_ID == "1237c64eb05da167644a4bff9d681bde9b9daa54")
#ubi_1618797b285d4650872ef897004644b2990e6d73 <- biomG_ubi %>% filter(OTU_ID == "1618797b285d4650872ef897004644b2990e6d73")
#ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 <- biomG_ubi %>% filter(OTU_ID == "3dcd4e930322140ea6d4c075826c2dfbb7124553")
#ubi_6e447d49df572a0622adad2e68409bcf687108e8 <- biomG_ubi %>% filter(OTU_ID == "6e447d49df572a0622adad2e68409bcf687108e8")
#ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- biomG_ubi %>% filter(OTU_ID == "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d")
#ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- biomG_ubi %>% filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf")
#ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f <- biomG_ubi %>% filter(OTU_ID == "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f")
#ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d <- biomG_ubi %>% filter(OTU_ID == "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d")

#these are the 10 ubiquitious OTUs
biomG_ubi %>% distinct(OTU_ID)

#adds the number of Pro, Syn sequences in each biogeotraces sample to each dataframe
#so that there is a row for each biogeotraces sample for each ubiquitious OTU
#ubi_02ba11823bfe93a0250096d95e07e745e40736bc <- ubi_02ba11823bfe93a0250096d95e07e745e40736bc %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
#ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 <- ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
#ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 <- ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
#ubi_1618797b285d4650872ef897004644b2990e6d73 <- ubi_1618797b285d4650872ef897004644b2990e6d73 %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
#ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 <- ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
#ubi_6e447d49df572a0622adad2e68409bcf687108e8 <- ubi_6e447d49df572a0622adad2e68409bcf687108e8 %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
#ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
#ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
#ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f <- ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
#ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d <- ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))

#these are the 10 ubiquitious OTUs
biomG_ubi %>% distinct(OTU_ID)

#now, there is a row for each biogeotraces sample for each ubiquitious OTU
#ubi_02ba11823bfe93a0250096d95e07e745e40736bc %>% nrow()
#ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 %>% nrow() 
#ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 %>% nrow() 
#ubi_1618797b285d4650872ef897004644b2990e6d73 %>% nrow() 
#ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 %>% nrow() 
#ubi_6e447d49df572a0622adad2e68409bcf687108e8 %>% nrow() 
#ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d %>% nrow() 
#ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf %>% nrow() 
#ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f %>% nrow() 
#ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d %>% nrow() 

#these are the 10 ubiquitious OTUs
biomG_ubi %>% distinct(OTU_ID)

#makes a variable with the OTU ID in each dataframe
#ubi_02ba11823bfe93a0250096d95e07e745e40736bc <- ubi_02ba11823bfe93a0250096d95e07e745e40736bc %>% 
  #mutate(OTU_ID = "02ba11823bfe93a0250096d95e07e745e40736bc")
#ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 <- ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 %>% 
  #mutate(OTU_ID = "08d46233bb26721e2336eb640c2dca0dd94497d3")
#ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 <- ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 %>% 
  #mutate(OTU_ID = "1237c64eb05da167644a4bff9d681bde9b9daa54")
#ubi_1618797b285d4650872ef897004644b2990e6d73 <- ubi_1618797b285d4650872ef897004644b2990e6d73 %>% 
  #mutate(OTU_ID = "1618797b285d4650872ef897004644b2990e6d73")
#ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 <- ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 %>% 
  #mutate(OTU_ID = "3dcd4e930322140ea6d4c075826c2dfbb7124553")
#ubi_6e447d49df572a0622adad2e68409bcf687108e8 <- ubi_6e447d49df572a0622adad2e68409bcf687108e8 %>% 
  #mutate(OTU_ID = "6e447d49df572a0622adad2e68409bcf687108e8")
#ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d %>% 
  #mutate(OTU_ID = "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d")
#ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf %>% 
  #mutate(OTU_ID = "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf")
#ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f <- ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f %>% 
  #mutate(OTU_ID = "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f")
#ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d <- ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d %>% 
  #mutate(OTU_ID = "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d")

#these are the 10 ubiquitious OTUs
biomG_ubi %>% distinct(OTU_ID)

#binds dataframes together
#data <- bind_rows(ubi_02ba11823bfe93a0250096d95e07e745e40736bc,
  #ubi_08d46233bb26721e2336eb640c2dca0dd94497d3,
  #ubi_1237c64eb05da167644a4bff9d681bde9b9daa54, 
  #ubi_1618797b285d4650872ef897004644b2990e6d73, 
  #ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553, 
  #ubi_6e447d49df572a0622adad2e68409bcf687108e8, 
  #ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d, 
  #ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf, 
  #ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f, 
  #ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d) 

head(biomG_ubi)
head(biomG_proSyn)

biomG_ubi %>% distinct(OTU_ID) %>% nrow()

biomG_ubi %>% distinct(biogeotracesSample) %>% nrow()
biomG_proSyn %>% distinct(biogeotracesSample) %>% nrow()

#adds number of pro, syn sequences in each biogeotraces sample to the number
#of sequences of each OTU I am interested in
nrow(biomG_ubi)
data <- biomG_ubi %>% left_join(biomG_proSyn, by = c("biogeotracesSample"))
nrow(data)

head(data)

#these are the biogeotraces samples that don't have any pro, syn sequences
data %>% filter(is.na(numProSynSeq)) %>% distinct(biogeotracesSample) %>% nrow()
data %>% filter(!is.na(numProSynSeq)) %>% distinct(biogeotracesSample) %>% nrow()

data %>% distinct(OTU_ID) %>% nrow()

data %>% filter(is.na(numSeqOfUbiOTU))

data %>% head()
data %>% filter(is.na(OTU_ID))
data %>% filter(is.na(biogeotracesSample))
data %>% filter(is.na(numSeqOfUbiOTU))
data %>% filter(is.na(numProSynSeq))

#excludes biogeotraces samples that don't have any pro, syn sequences even if 
#the heterotroph OTU from culture is present in biogeotraces
#I only want to take correlation based on the instances where both the heterotroph 
#OTU and Pro/Syn are present in the biogeotraces sample, because this is what 
#Sean suggested
data %>% filter(is.na(numProSynSeq))
data <- data %>% filter(!is.na(numProSynSeq))

str(data)
data %>% head()

data %>% nrow()
data %>% distinct(OTU_ID) %>% nrow()

#the number of biogeotraces samples each OTU is found in, for which the biogeotraces 
#sample also has Pro/Syn sequences
data %>% group_by(OTU_ID) %>% summarize(n = n())
data %>% group_by(OTU_ID) %>% distinct(biogeotracesSample) %>% summarize(n = n())

data %>% head()

#adds number of sequences in each biogeotraces sample, including pro and syn 
#to the number of sequences of the ubiquitious heterotrophs in and number 
#of pro, syn sequences in each biogeotraces sample
head(data)
head(numSeqsBioGeoSample)
nrow(data)
data <- data %>% left_join(numSeqsBioGeoSample, by = c("biogeotracesSample"))
nrow(data)

head(data)
data %>% filter(is.na(totalSeqsBioGeo))

#calculates the proportion of sequences that are in each heterotroph OTU 
#and that are pro, syn sequences in each biogeotraces sample
data <- data %>% mutate(propSeqOfUbiOTU = numSeqOfUbiOTU/totalSeqsBioGeo)
data <- data %>% mutate(propProSynSeq = numProSynSeq/totalSeqsBioGeo)

head(data)
data %>% filter(is.na(propSeqOfUbiOTU))
data %>% filter(is.na(propProSynSeq))

#spearman correlation between relative abundance of OTUs of interest and relative 
#abundance of Pro, Syn sequences in each biogeotraces sample 
#the OTUs of interest are the OTUs that are ubiquitious in sean culture samples
df <- cor(data$propSeqOfUbiOTU, data$propProSynSeq, method = "spearman")

df <- df %>% as.data.frame()

df

colnames(df) <- "Correlation between 10 Ubiquitious OTUs and Pro., Syn. relative abundance"

df

write_csv(df, "10UbiquitiousOTUsAcrossCulturesGroupedTogether_howTheirAbundanceCorrelatesWithProSynAbundanceAcrossBiogeotraces.csv")


#these are the 10 ubiquitious OTUs
data %>% distinct(OTU_ID)

#calculates correlation for each OTU separately
#spearman correlation between the abundance of OTU of interest and the number of 
#Pro, Syn sequences in each biogeotraces sample 


library(psych)

cor_02ba11823bfe93a0250096d95e07e745e40736bc <- 
    corr.test(data %>% filter(OTU_ID == "02ba11823bfe93a0250096d95e07e745e40736bc") %>% select(propSeqOfUbiOTU), 
    data %>% filter(OTU_ID == "02ba11823bfe93a0250096d95e07e745e40736bc") %>% select(propProSynSeq), method = "spearman", adjust = "none")

cor_08d46233bb26721e2336eb640c2dca0dd94497d3 <-
    corr.test(data %>% filter(OTU_ID == "08d46233bb26721e2336eb640c2dca0dd94497d3") %>% select(propSeqOfUbiOTU), 
    data %>% filter(OTU_ID == "08d46233bb26721e2336eb640c2dca0dd94497d3") %>% select(propProSynSeq), method = "spearman", adjust = "none")

cor_1237c64eb05da167644a4bff9d681bde9b9daa54 <-
    corr.test(data %>% filter(OTU_ID == "1237c64eb05da167644a4bff9d681bde9b9daa54") %>% select(propSeqOfUbiOTU), 
    data %>% filter(OTU_ID == "1237c64eb05da167644a4bff9d681bde9b9daa54") %>% select(propProSynSeq), method = "spearman", adjust = "none")

cor_1618797b285d4650872ef897004644b2990e6d73 <- 
    corr.test(data %>% filter(OTU_ID == "1618797b285d4650872ef897004644b2990e6d73") %>% select(propSeqOfUbiOTU), 
    data %>% filter(OTU_ID == "1618797b285d4650872ef897004644b2990e6d73") %>% select(propProSynSeq), method = "spearman", adjust = "none")

cor_3dcd4e930322140ea6d4c075826c2dfbb7124553 <- 
    corr.test(data %>% filter(OTU_ID == "3dcd4e930322140ea6d4c075826c2dfbb7124553") %>% select(propSeqOfUbiOTU), 
    data %>% filter(OTU_ID == "3dcd4e930322140ea6d4c075826c2dfbb7124553") %>% select(propProSynSeq), method = "spearman", adjust = "none")


cor_6e447d49df572a0622adad2e68409bcf687108e8 <- 
    corr.test(data %>% filter(OTU_ID == "6e447d49df572a0622adad2e68409bcf687108e8") %>% select(propSeqOfUbiOTU), 
    data %>% filter(OTU_ID == "6e447d49df572a0622adad2e68409bcf687108e8") %>% select(propProSynSeq), method = "spearman", adjust = "none")

cor_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- 
    corr.test(data %>% filter(OTU_ID == "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d") %>% select(propSeqOfUbiOTU), 
    data %>% filter(OTU_ID == "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d") %>% select(propProSynSeq), method = "spearman", adjust = "none")

cor_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- 
    corr.test(data %>% filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf") %>% select(propSeqOfUbiOTU), 
    data %>% filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf") %>% select(propProSynSeq), method = "spearman", adjust = "none")

cor_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f <-
    corr.test(data %>% filter(OTU_ID == "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f") %>% select(propSeqOfUbiOTU), 
    data %>% filter(OTU_ID == "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f") %>% select(propProSynSeq), method = "spearman", adjust = "none")

cor_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d <-
    corr.test(data %>% filter(OTU_ID == "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d") %>% select(propSeqOfUbiOTU), 
    data %>% filter(OTU_ID == "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d") %>% select(propProSynSeq), method = "spearman", adjust = "none")





#these are the 10 ubiquitious OTUs
data %>% distinct(OTU_ID)

cor_02ba11823bfe93a0250096d95e07e745e40736bc_pvalue <- cor_02ba11823bfe93a0250096d95e07e745e40736bc$p %>% 
  as.data.frame() %>% mutate(OTU = "02ba11823bfe93a0250096d95e07e745e40736bc")

cor_08d46233bb26721e2336eb640c2dca0dd94497d3_pvalue <- cor_08d46233bb26721e2336eb640c2dca0dd94497d3$p %>% 
  as.data.frame() %>% mutate(OTU = "08d46233bb26721e2336eb640c2dca0dd94497d3")

cor_1237c64eb05da167644a4bff9d681bde9b9daa54_pvalue <- cor_1237c64eb05da167644a4bff9d681bde9b9daa54$p %>% 
  as.data.frame() %>% mutate(OTU = "1237c64eb05da167644a4bff9d681bde9b9daa54")

cor_1618797b285d4650872ef897004644b2990e6d73_pvalue <- cor_1618797b285d4650872ef897004644b2990e6d73$p %>% 
  as.data.frame() %>% mutate(OTU = "1618797b285d4650872ef897004644b2990e6d73")

cor_3dcd4e930322140ea6d4c075826c2dfbb7124553_pvalue <- cor_3dcd4e930322140ea6d4c075826c2dfbb7124553$p %>% 
  as.data.frame() %>% mutate(OTU = "3dcd4e930322140ea6d4c075826c2dfbb7124553")

cor_6e447d49df572a0622adad2e68409bcf687108e8_pvalue <- cor_6e447d49df572a0622adad2e68409bcf687108e8$p %>% 
  as.data.frame() %>% mutate(OTU = "6e447d49df572a0622adad2e68409bcf687108e8")

cor_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d_pvalue <- cor_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d$p %>% 
  as.data.frame() %>% mutate(OTU = "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d")

cor_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf_pvalue <- cor_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf$p %>% 
  as.data.frame() %>% mutate(OTU = "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf")

cor_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f_pvalue <- cor_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f$p %>% 
  as.data.frame() %>% mutate(OTU = "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f")

cor_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d_pvalue <- cor_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d$p %>% 
  as.data.frame() %>% mutate(OTU = "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d")



#these are the 10 ubiquitious OTUs
data %>% distinct(OTU_ID)


cor_02ba11823bfe93a0250096d95e07e745e40736bc <- cor_02ba11823bfe93a0250096d95e07e745e40736bc$r %>% 
  as.data.frame() %>% mutate(OTU = "02ba11823bfe93a0250096d95e07e745e40736bc")

cor_08d46233bb26721e2336eb640c2dca0dd94497d3 <- cor_08d46233bb26721e2336eb640c2dca0dd94497d3$r %>% 
  as.data.frame() %>% mutate(OTU = "08d46233bb26721e2336eb640c2dca0dd94497d3")

cor_1237c64eb05da167644a4bff9d681bde9b9daa54 <- cor_1237c64eb05da167644a4bff9d681bde9b9daa54$r %>% 
  as.data.frame() %>% mutate(OTU = "1237c64eb05da167644a4bff9d681bde9b9daa54")

cor_1618797b285d4650872ef897004644b2990e6d73 <- cor_1618797b285d4650872ef897004644b2990e6d73$r %>% 
  as.data.frame() %>% mutate(OTU = "1618797b285d4650872ef897004644b2990e6d73")

cor_3dcd4e930322140ea6d4c075826c2dfbb7124553 <- cor_3dcd4e930322140ea6d4c075826c2dfbb7124553$r %>% 
  as.data.frame() %>% mutate(OTU = "3dcd4e930322140ea6d4c075826c2dfbb7124553")

cor_6e447d49df572a0622adad2e68409bcf687108e8 <- cor_6e447d49df572a0622adad2e68409bcf687108e8$r %>% 
  as.data.frame() %>% mutate(OTU = "6e447d49df572a0622adad2e68409bcf687108e8")

cor_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- cor_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d$r %>% 
  as.data.frame() %>% mutate(OTU = "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d")

cor_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- cor_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf$r %>% 
  as.data.frame() %>% mutate(OTU = "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf")

cor_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f <- cor_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f$r %>% 
  as.data.frame() %>% mutate(OTU = "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f")

cor_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d <- cor_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d$r %>% 
  as.data.frame() %>% mutate(OTU = "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d")


#these are the 10 ubiquitious OTUs
data %>% distinct(OTU_ID)

cor <- bind_rows(cor_02ba11823bfe93a0250096d95e07e745e40736bc,
  cor_08d46233bb26721e2336eb640c2dca0dd94497d3,
  cor_1237c64eb05da167644a4bff9d681bde9b9daa54,
  cor_1618797b285d4650872ef897004644b2990e6d73,
  cor_3dcd4e930322140ea6d4c075826c2dfbb7124553,
  cor_6e447d49df572a0622adad2e68409bcf687108e8,
  cor_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d,
  cor_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf,
  cor_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f,
  cor_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d)

pvalue <- bind_rows(cor_02ba11823bfe93a0250096d95e07e745e40736bc_pvalue,
                 cor_08d46233bb26721e2336eb640c2dca0dd94497d3_pvalue,
                 cor_1237c64eb05da167644a4bff9d681bde9b9daa54_pvalue,
                 cor_1618797b285d4650872ef897004644b2990e6d73_pvalue,
                 cor_3dcd4e930322140ea6d4c075826c2dfbb7124553_pvalue,
                 cor_6e447d49df572a0622adad2e68409bcf687108e8_pvalue,
                 cor_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d_pvalue,
                 cor_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf_pvalue,
                 cor_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f_pvalue,
                 cor_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d_pvalue)

#renames variable names
head(cor)
cor
colnames(cor)
colnames(cor)[1] <- "Correlation between OTU and Pro., Syn. relative abundance"
colnames(cor)[2] <- "OTU_ID"

head(pvalue)
colnames(pvalue)
colnames(pvalue)[1] <- "P-value"
colnames(pvalue)[2] <- "OTU_ID"

head(cor)
head(tax)

#adds taxonomy to the cor dataframe
nrow(cor)
cor <- cor %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(cor)

cor %>% filter(is.na(Taxon)) %>% nrow()
head(cor)

cor %>% distinct(Taxon)

#makes a taxonomy variable without "k__, p__..."
cor <- cor %>% mutate(newTaxon = str_replace_all(Taxon, "([a-z]__)", ""))

#gets rid of ";" in fixed taxonomy variable
cor <- cor %>% mutate(newTaxon = str_replace_all(newTaxon, ";", ""))
#ubi %>% distinct(Taxon, newTaxon) %>% View()

#makes a shortened taxonomy variable
cor <- cor %>% mutate(shortTaxon = str_extract(newTaxon, "[a-zA-Z0-9\\-\\_]*\\s*$"))
#ubi %>% distinct(Taxon, shortTaxon) %>% View()

cor %>% distinct(shortTaxon)

#gets rid of extra spaces at the end of shortened taxonomies
str_extract(cor$shortTaxon, "\\s$")
cor <- cor %>% mutate(shortTaxon = str_replace_all(shortTaxon, "\\s$", ""))

cor %>% distinct(shortTaxon)
#cor %>% distinct(Taxon, shortTaxon) %>% View()

cor %>% head()
colnames(cor)

#gets rid of unnecessary variables
cor <- cor %>% select(1,2,5)

cor %>% head()

colnames(cor)
colnames(cor)[3] <- "Taxonomy"

cor %>% head()
cor %>% distinct(OTU_ID) %>% nrow()
nrow(cor)

head(cor)
head(pvalue)

#add p-value of correlations
nrow(cor)
cor <- cor %>% left_join(pvalue, by = c("OTU_ID"))
nrow(cor)

cor
colnames(cor)

#puts variables in the right order
cor <- cor %>% select(OTU_ID, Taxonomy, `Correlation between OTU and Pro., Syn. relative abundance`, `P-value`)

cor
head(cor)

#calculated adjusted p-value 
#from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6099145/
cor$`Adjusted p-value` <- p.adjust(cor$`P-value`, method="fdr")

cor
cor %>% select(`P-value`, `Adjusted p-value`)
cor %>% filter(is.na(`P-value`))

data %>% filter(OTU_ID == "1618797b285d4650872ef897004644b2990e6d73")
corr.test(data %>% filter(OTU_ID == "1618797b285d4650872ef897004644b2990e6d73") %>% select(propSeqOfUbiOTU), 
          data %>% filter(OTU_ID == "1618797b285d4650872ef897004644b2990e6d73") %>% select(propProSynSeq), method = "spearman", adjust = "none")

write_csv(cor, "ubiquitiousOTUsAcrossCultures_howTheirAbundanceCorrelatesWithProSynAbundanceAcrossBiogeotraces.csv")


