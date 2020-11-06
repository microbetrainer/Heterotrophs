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

biomG_ubi %>% head()
biomG_ubi %>% distinct(biogeotracesSample) %>% nrow()

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
missingBiogeo <- biomG %>% anti_join(biomG_proSyn, by = c("biogeotracesSample")) %>% distinct(biogeotracesSample) 

biomG_proSyn %>% head()
missingBiogeo %>% head()

#makes numProSynSeq value 0 for the biogeotraces samples that are missing from biomG_proSyn
missingBiogeo <- missingBiogeo %>% mutate(numProSynSeq = 0)
missingBiogeo %>% head()

#adds missing biogeotraces samples to biomG_proSyn
nrow(biomG_proSyn)
biomG_proSyn <- bind_rows(biomG_proSyn, missingBiogeo)
nrow(biomG_proSyn)

#now all of the biogeotraces samples are in biomG_proSyn
biomG_proSyn %>% head()
biomG_proSyn %>% nrow()
biomG_proSyn %>% distinct(biogeotracesSample) %>% nrow()

#these are the 10 ubiquitious OTUs
biomG_ubi %>% distinct(OTU_ID)

#makes a dataframe for each OTU in biomG_ubi
ubi_02ba11823bfe93a0250096d95e07e745e40736bc <- biomG_ubi %>% filter(OTU_ID == "02ba11823bfe93a0250096d95e07e745e40736bc")
ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 <- biomG_ubi %>% filter(OTU_ID == "08d46233bb26721e2336eb640c2dca0dd94497d3")
ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 <- biomG_ubi %>% filter(OTU_ID == "1237c64eb05da167644a4bff9d681bde9b9daa54")
ubi_1618797b285d4650872ef897004644b2990e6d73 <- biomG_ubi %>% filter(OTU_ID == "1618797b285d4650872ef897004644b2990e6d73")
ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 <- biomG_ubi %>% filter(OTU_ID == "3dcd4e930322140ea6d4c075826c2dfbb7124553")
ubi_6e447d49df572a0622adad2e68409bcf687108e8 <- biomG_ubi %>% filter(OTU_ID == "6e447d49df572a0622adad2e68409bcf687108e8")
ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- biomG_ubi %>% filter(OTU_ID == "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d")
ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- biomG_ubi %>% filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf")
ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f <- biomG_ubi %>% filter(OTU_ID == "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f")
ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d <- biomG_ubi %>% filter(OTU_ID == "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d")

#these are the 10 ubiquitious OTUs
biomG_ubi %>% distinct(OTU_ID)

#adds the number of Pro, Syn sequences in each biogeotraces sample to each dataframe
#so that there is a row for each biogeotraces sample for each ubiquitious OTU
ubi_02ba11823bfe93a0250096d95e07e745e40736bc <- ubi_02ba11823bfe93a0250096d95e07e745e40736bc %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 <- ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 <- ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
ubi_1618797b285d4650872ef897004644b2990e6d73 <- ubi_1618797b285d4650872ef897004644b2990e6d73 %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 <- ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
ubi_6e447d49df572a0622adad2e68409bcf687108e8 <- ubi_6e447d49df572a0622adad2e68409bcf687108e8 %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f <- ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))
ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d <- ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d %>% full_join(biomG_proSyn, by = c("biogeotracesSample"))

#these are the 10 ubiquitious OTUs
biomG_ubi %>% distinct(OTU_ID)

#now, there is a row for each biogeotraces sample for each ubiquitious OTU
ubi_02ba11823bfe93a0250096d95e07e745e40736bc %>% nrow()
ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 %>% nrow() 
ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 %>% nrow() 
ubi_1618797b285d4650872ef897004644b2990e6d73 %>% nrow() 
ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 %>% nrow() 
ubi_6e447d49df572a0622adad2e68409bcf687108e8 %>% nrow() 
ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d %>% nrow() 
ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf %>% nrow() 
ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f %>% nrow() 
ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d %>% nrow() 

#these are the 10 ubiquitious OTUs
biomG_ubi %>% distinct(OTU_ID)

#makes a variable with the OTU ID in each dataframe
ubi_02ba11823bfe93a0250096d95e07e745e40736bc <- ubi_02ba11823bfe93a0250096d95e07e745e40736bc %>% 
  mutate(OTU_ID = "02ba11823bfe93a0250096d95e07e745e40736bc")
ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 <- ubi_08d46233bb26721e2336eb640c2dca0dd94497d3 %>% 
  mutate(OTU_ID = "08d46233bb26721e2336eb640c2dca0dd94497d3")
ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 <- ubi_1237c64eb05da167644a4bff9d681bde9b9daa54 %>% 
  mutate(OTU_ID = "1237c64eb05da167644a4bff9d681bde9b9daa54")
ubi_1618797b285d4650872ef897004644b2990e6d73 <- ubi_1618797b285d4650872ef897004644b2990e6d73 %>% 
  mutate(OTU_ID = "1618797b285d4650872ef897004644b2990e6d73")
ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 <- ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553 %>% 
  mutate(OTU_ID = "3dcd4e930322140ea6d4c075826c2dfbb7124553")
ubi_6e447d49df572a0622adad2e68409bcf687108e8 <- ubi_6e447d49df572a0622adad2e68409bcf687108e8 %>% 
  mutate(OTU_ID = "6e447d49df572a0622adad2e68409bcf687108e8")
ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d %>% 
  mutate(OTU_ID = "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d")
ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf %>% 
  mutate(OTU_ID = "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf")
ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f <- ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f %>% 
  mutate(OTU_ID = "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f")
ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d <- ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d %>% 
  mutate(OTU_ID = "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d")

#these are the 10 ubiquitious OTUs
biomG_ubi %>% distinct(OTU_ID)

#binds dataframes together
data <- bind_rows(ubi_02ba11823bfe93a0250096d95e07e745e40736bc,
  ubi_08d46233bb26721e2336eb640c2dca0dd94497d3,
  ubi_1237c64eb05da167644a4bff9d681bde9b9daa54, 
  ubi_1618797b285d4650872ef897004644b2990e6d73, 
  ubi_3dcd4e930322140ea6d4c075826c2dfbb7124553, 
  ubi_6e447d49df572a0622adad2e68409bcf687108e8, 
  ubi_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d, 
  ubi_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf, 
  ubi_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f, 
  ubi_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d) 

data %>% distinct(OTU_ID) %>% nrow()

#makes numSeqOfUbiOTU value 0 if the heterotroph OTU is not present in the biogeotraces
#sample
data %>% filter(is.na(numSeqOfUbiOTU))
data <- data %>% mutate(numSeqOfUbiOTU = ifelse(is.na(numSeqOfUbiOTU), 0, numSeqOfUbiOTU))

data %>% head()
data %>% filter(is.na(OTU_ID))
data %>% filter(is.na(biogeotracesSample))
data %>% filter(is.na(numSeqOfUbiOTU))
data %>% filter(is.na(numProSynSeq))

str(data)
data %>% head()

data %>% nrow()
data %>% distinct(OTU_ID) %>% nrow()
data %>% group_by(OTU_ID) %>% summarize(n = n()) %>% filter(n != 610)

data %>% head()


data %>% str()
data %>% head()

#gathers data into long form
dataG <- data %>% gather(numSeqOfUbiOTU:numProSynSeq, key = "type", value = "value")

dataG %>% str()
dataG %>% head()

dataG %>% group_by(OTU_ID, biogeotracesSample) %>% summarize(n = n()) %>% filter(n != 2)

dataG %>% head()

dataG %>% summarize(min = min(value), max = max(value))

#dataG %>% ggplot(aes(x = biogeotracesSample, y = value, color = type)) + geom_point() + 
  #facet_wrap(~OTU_ID, scales = 'free') + scale_y_log10()

dataG %>% group_by(type) %>% summarize(mean = mean(value))
dataG %>% group_by(type) %>% summarize(min = min(value), max = max(value))

data %>% head()
data %>% filter(numSeqOfUbiOTU != 0) %>% summarize(min = min(numSeqOfUbiOTU))
data %>% nrow()
data %>% filter(numSeqOfUbiOTU == 1) %>% nrow()

data %>% head()
data %>% filter(numProSynSeq != 0) %>% summarize(min = min(numProSynSeq))
data %>% nrow()
data %>% filter(numProSynSeq == 1) %>% nrow()

#data %>% ggplot(aes(x = biogeotracesSample, y = numSeqOfUbiOTU/numProSynSeq)) + geom_point() + 
  #facet_wrap(~OTU_ID, scales = 'free')

#loads the correlation, p-value, and adjusted p-value for each heterotroph OTU and Pro, Syn 
#sequence abundances across biogeotraces sample
#correlations were only run on instances in which both the heterotroph OTU and the Pro, Syn 
#OTU were present in the biogeotraces sample
#from howDoBiogeotracesOTUsCorrelateWithPro.R (checked)
cor <- read_csv("ubiquitiousOTUsAcrossCultures_howTheirAbundanceCorrelatesWithProSynAbundanceAcrossBiogeotraces.csv")

head(cor)
nrow(cor)
cor

head(data)
data %>% filter(is.na(numSeqOfUbiOTU))
data %>% filter(is.na(numProSynSeq))

data %>% distinct(OTU_ID)
data %>% anti_join(cor, by = c("OTU_ID")) %>% distinct(OTU_ID)
cor %>% anti_join(data, by = c("OTU_ID")) %>% distinct(OTU_ID)

#adds correlation data to sequence abundance of heterotroph and Pro/Syn data
nrow(data)
data <- data %>% left_join(cor, by = c("OTU_ID"))
nrow(data)

data %>% filter(is.na(Taxonomy))

head(data)

data %>% distinct(OTU_ID)

#makes NA correlation values "" so I can do a string operation with them
data %>% distinct(OTU_ID, `Correlation between OTU and Pro., Syn. relative abundance`)
data %>% filter(is.na(`Correlation between OTU and Pro., Syn. relative abundance`)) %>% distinct(OTU_ID)
data <- data %>% mutate(`Correlation between OTU and Pro., Syn. relative abundance` = 
                          ifelse(is.na(`Correlation between OTU and Pro., Syn. relative abundance`), "", `Correlation between OTU and Pro., Syn. relative abundance`))
data %>% distinct(OTU_ID, `Correlation between OTU and Pro., Syn. relative abundance`)


##adds asterisks based on the p-value
data %>% distinct(OTU_ID, `P-value`) %>% arrange(`P-value`)

data <- data %>% mutate(significance = ifelse(`P-value` <= 0.05, "*", ""))
data <- data %>% mutate(significance = ifelse(`P-value` <= 0.01, "**", significance))
data <- data %>% mutate(significance = ifelse(`P-value` <= 0.001, "***", significance))
data <- data %>% mutate(significance = ifelse(`P-value` <= 0.0001, "****", significance))

data %>% distinct(OTU_ID, `P-value`, significance) %>% arrange(`P-value`)

##adds asterisks based on the adjusted p-value
data %>% distinct(OTU_ID, `Adjusted p-value`) %>% arrange(`Adjusted p-value`)

data <- data %>% mutate(significance_adjusted = ifelse(`Adjusted p-value` <= 0.05, "*", ""))
data <- data %>% mutate(significance_adjusted = ifelse(`Adjusted p-value` <= 0.01, "**", significance_adjusted))
data <- data %>% mutate(significance_adjusted = ifelse(`Adjusted p-value` <= 0.001, "***", significance_adjusted))
data <- data %>% mutate(significance_adjusted = ifelse(`Adjusted p-value` <= 0.0001, "****", significance_adjusted))

data %>% distinct(OTU_ID, `Adjusted p-value`, significance_adjusted) %>% arrange(`Adjusted p-value`)

data %>% distinct(OTU_ID, significance, significance_adjusted)
data %>% filter(significance != significance_adjusted) %>% 
  distinct(OTU_ID, `P-value`, `Adjusted p-value`, significance, significance_adjusted)

data %>% distinct(OTU_ID, `Correlation between OTU and Pro., Syn. relative abundance`)

#makes variable with taxonomy and correlation value
data <- data %>% mutate(taxCorr = str_c(Taxonomy, `Correlation between OTU and Pro., Syn. relative abundance`, sep = ":\n"))

data %>% distinct(OTU_ID, taxCorr, significance)

#replaces NA significance values with "NA" so that I can do a string operation with them
data <- data %>% mutate(significance = ifelse(is.na(significance), "", significance))
data <- data %>% mutate(significance_adjusted = ifelse(is.na(significance_adjusted), "", significance_adjusted))

data %>% distinct(OTU_ID, taxCorr, significance, significance_adjusted)

#combines taxonomy, correlation, and adjusted significance
data <- data %>% mutate(taxCorr = str_c(taxCorr, significance_adjusted))

data %>% distinct(OTU_ID, taxCorr)

data %>% distinct(OTU_ID, `P-value`, `Adjusted p-value`, significance_adjusted, taxCorr)

str(data)

data %>% filter(is.na(numSeqOfUbiOTU))
data %>% filter(is.na(numProSynSeq))

data %>% filter(numSeqOfUbiOTU == 0)
data %>% filter(numProSynSeq == 0)

data %>% filter(numSeqOfUbiOTU == 0 & numProSynSeq == 0) %>% nrow()

#excludes instances when both the heterotroph OTU and Pro/Syn were absent
#from the biogeotraces sample
nrow(data)
data <- data %>% filter(!(numSeqOfUbiOTU == 0 & numProSynSeq == 0))
nrow(data)

6100-5455

#gets rid of ":\n" if it is at the end of taxCorr
data %>% distinct(OTU_ID, taxCorr)
data <- data %>% mutate(taxCorr = ifelse(str_detect(taxCorr, ":\n$"), str_replace(taxCorr, ":\n", ""), taxCorr))
data %>% distinct(OTU_ID, taxCorr)

#adds number of sequences, including Pro and Syn sequences in each biogeotraces 
#sample to data
head(data)
nrow(data)
data <- data %>% left_join(numSeqsBioGeoSample, by = c("biogeotracesSample"))
nrow(data)

head(data)
data %>% filter(is.na(totalSeqsBioGeo))

str(data$numSeqOfUbiOTU)
str(data$numProSynSeq)

#add psuedocount so I can take log
data %>% filter(numSeqOfUbiOTU == 0 & numProSynSeq == 0) %>% nrow()
#data <- data %>% mutate(numSeqOfUbiOTU = ifelse(numSeqOfUbiOTU == 0, 1, numSeqOfUbiOTU))
#data <- data %>% mutate(numProSynSeq = ifelse(numProSynSeq == 0, 1, numProSynSeq))

str(data$totalSeqsBioGeo)

#normalizes number of ubiquitious heterotroph OTUs and number of pro, syn sequences 
#by total number of sequences, including pro and syn in each biogeotraces sample
data <- data %>% mutate(propSeqOfUbiOTU = numSeqOfUbiOTU/totalSeqsBioGeo)
data <- data %>% mutate(propProSynSeq = numProSynSeq/totalSeqsBioGeo)

data %>% filter(is.na(propSeqOfUbiOTU))
data %>% filter(is.na(propProSynSeq))

data %>% filter(propSeqOfUbiOTU == 0)
data %>% filter(propProSynSeq == 0)

head(data)
data %>% filter(propSeqOfUbiOTU != 0) %>% summarize(min = min(propSeqOfUbiOTU))
data %>% filter(propProSynSeq != 0) %>% summarize(min = min(propProSynSeq))
10^-5
10^-5 < 0.000176
10^-5 < 0.000309

#add psuedocount so I can take log
data <- data %>% mutate(propSeqOfUbiOTU = ifelse(propSeqOfUbiOTU == 0, 10^-5, propSeqOfUbiOTU))
data <- data %>% mutate(propProSynSeq = ifelse(propProSynSeq == 0, 10^-5, propProSynSeq))

data %>% filter(is.na(propSeqOfUbiOTU))
data %>% filter(is.na(propProSynSeq))

data %>% filter(propSeqOfUbiOTU == 0)
data %>% filter(propProSynSeq == 0)

#data %>% filter(propSeqOfUbiOTU != 0) %>% summarize(min = min(propSeqOfUbiOTU))
#data %>% filter(propProSynSeq != 0) %>% summarize(min = min(propProSynSeq))

#why do the plots look so similar??
data %>% ggplot(aes(x = log10(propProSynSeq), y = log10(propSeqOfUbiOTU))) + geom_point() + 
  facet_wrap(~taxCorr, scales = 'free') + 
  labs(x = "Log10 relative abundance of Pro. and Syn. sequences", y = "Log10 relative abundance of heterotroph") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=14))

ggsave("ubiquitiousOTUsAcrossCultures_howTheirAbundanceCorrelatesWithProSynAbundanceAcrossBiogeotraces.png", dpi = 300, height = 8, width = 10)
ggsave("ubiquitiousOTUsAcrossCultures_howTheirAbundanceCorrelatesWithProSynAbundanceAcrossBiogeotraces.svg", dpi = 300, height = 8, width = 10)



