library(tidyverse)
library(lubridate)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads the sean and biogeotraces sequences present in each cluster
#from clusterSeanAndBiogeotracesSeqs.txt (checked)
biomG <- read_csv("combinedBiogeotracesBiom.csv")

#28022 clusters (correct)
#1184489 sequences (correct)
head(biomG)
biomG %>% distinct(OTU_ID) %>% nrow()
biomG %>% distinct(sample) %>% nrow()
biomG %>% nrow()

str(biomG)

biomG %>% distinct(present)

#each sequence is only in one cluster
head(biomG)
biomG %>% group_by(sample) %>% summarize(n = n()) %>% filter(n != 1)

#loads the abundance of each sean sequence in each culture
#this includes pro and syn 
#sequences with prop < .002 have been excluded
#from heatmap_noContamination_noLowAbundance_noMock_noProSyn_without1223_propReads.R
seqG <- read_csv("forSeanComparisonToBiogeotraces.csv")
seqG %>% head()

#this includes Pro and Syn sequences and the JW3 and 1223 samples
seqG %>% distinct(sequence) %>% nrow()

seqG %>% filter(sample == "JW3")
seqG %>% filter(sample == "1223")

#excludes JW3 and 1223
seqG <- seqG %>% filter(sample != "JW3")
seqG <- seqG %>% filter(sample != "1223")

#includes Pro and Syn
seqG %>% distinct(sample) %>% nrow()
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

#correct number of samples
seqG %>% distinct(sample) %>% nrow()

#loads the sequences after excluding pro and syn 
#sequences with prop < .002 have been excluded
#from getSequencesForTree_noContamination_noProSyn_withoutSample1223.R
seqsNoProSyn <- read.table("sequencesForTree_noContamination_noMock_noProSyn_no1223.txt")

#includes Pro and Syn
seqG %>% distinct(sequence) %>% nrow()

#doesn't include Pro and Syn 
seqsNoProSyn %>% head()
seqsNoProSyn %>% distinct(V1) %>% nrow()

head(seqG)
head(seqsNoProSyn)

#excludes Pro and Syn sequences from seqG
seqG <- seqG %>% semi_join(seqsNoProSyn, by = c("sequence" = "V1"))

#correct number of sequences and samples
#no pro, syn sequences in seqG now
seqG %>% distinct(sequence) %>% nrow()
seqG %>% distinct(sample) %>% nrow()

#biomG includes Pro and Syn sequences in sean samples
biomG %>% head()
biomG %>% filter(str_detect(sample, "sean"))
biomG %>% filter(str_detect(sample, "sean")) %>% nrow()

#makes the sean sequence IDs in biomG match seqG
biomG %>% filter(str_detect(sample, "sean")) %>% head()
biomG <- biomG %>% mutate(sample = ifelse(str_detect(sample, "sean"), str_replace(sample, "X", ""), sample))

biomG %>% filter(str_detect(sample, "sean")) %>% head()
biomG <- biomG %>% mutate(sample = ifelse(str_detect(sample, "sean"), str_c("contig_", sample), sample))

biomG %>% filter(str_detect(sample, "sean")) %>% head()
biomG <- biomG %>% mutate(sample = ifelse(str_detect(sample, "sean"), str_replace(sample, "sean", ""), sample))

#these are the sean sequences
biomG %>% filter(str_detect(sample, "contig"))
biomG %>% filter(str_detect(sample, "contig")) %>% nrow()

#finds the sean sequences in biomG that do not have matches in seqG
#these are the sean pro and syn sequences
head(biomG)
head(seqG)
biomG %>% filter(str_detect(sample, "contig")) %>% left_join(seqG, by = c("sample" = "sequence")) %>% filter(is.na(abundance)) %>% nrow()

#these are the clusters arranged by how many sequences they have across 
#biogeotraces and sean samples
biomG %>% group_by(OTU_ID) %>% summarize(n = n()) %>% arrange(desc(n))

#there are clusters in sean samples that have multiple sequences across 
#sean samples
biomG %>% filter(str_detect(sample, "contig")) %>% group_by(OTU_ID) %>% summarize(n = n()) %>% arrange(desc(n))

#gets just the clusters that sean sequences are in 
seanBiom <- biomG %>% filter(str_detect(sample, "contig"))

#gets just the clusters that biogeotraces sequences are in
biogeoBiom <- biomG %>% filter(!str_detect(sample, "contig"))

nrow(biomG)
nrow(seanBiom)
nrow(biogeoBiom)
nrow(seanBiom)+nrow(biogeoBiom)==nrow(biomG)

head(seqG)
head(seanBiom)

#seqG does not include Pro and Syn
seqG %>% distinct(sequence) %>% nrow()
seanBiom %>% distinct(sample) %>% nrow()

#adds the cluster each sean sequence is assigned to to the abundance of
#sean sequences across sean samples
nrow(seqG)
merged <- seqG %>% left_join(seanBiom, by = c("sequence" = "sample"))
nrow(merged)

head(merged)
merged %>% filter(is.na(present))

#gets rid of unnecessary variable
merged <- merged %>% select(-present)

#doesn't include pro and syn
head(merged)
merged %>% distinct(sequence) %>% nrow()

merged %>% filter(is.na(OTU_ID))

head(merged)

head(biogeoBiom)

#gets rid of unnecessary variable
biogeoBiom <- biogeoBiom %>% select(-present)

#there are 28008 OTUs present in biogeotraces
biogeoBiom %>% distinct(OTU_ID) %>% nrow()

#merged has one fewer OTU because this OTU has the pro and syn seqs
merged %>% distinct(OTU_ID) %>% nrow()
seanBiom %>% distinct(OTU_ID) %>% nrow()

merged %>% distinct(sequence) %>% nrow()
seanBiom %>% distinct(sample) %>% nrow()

#there are 27925 OTUs present in biogeotraces that are not present in sean samples
biogeoBiom %>% anti_join(merged, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

#there are 83 OTUs present in biogeotraces that are also present in sean samples
biogeoBiom %>% semi_join(merged, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

27925+83

#there are 97 OTUs present in sean samples
merged %>% distinct(OTU_ID) %>% nrow() 

#there are 14 OTUs present in sean samples that are not present in biogeotraces
merged %>% anti_join(biogeoBiom, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

83+14

head(merged)
head(biogeoBiom)

#adds the biogeotraces sequences that are in OTUs that are present in sean samples 
#to the abundance of OTUs in sean samples
nrow(merged)
merged <- merged %>% left_join(biogeoBiom, by = c("OTU_ID"))
nrow(merged)

head(merged)

#these are the 14 OTUs that are in sean samples that are not present in biogeotraces
merged %>% filter(is.na(sample.y)) %>% distinct(OTU_ID) %>% nrow()

merged %>% distinct(OTU_ID) %>% nrow()

head(merged)

#all of the samples have at least one sequence that corresponds to one of the biogeotraces sequences
merged %>% distinct(sample.x) %>% nrow()
merged %>% filter(!is.na(sample.y)) %>% distinct(sample.x) %>% nrow()

#renames merged column names
head(merged)
colnames(merged)
colnames(merged)[1] <- "sample" 
colnames(merged)[5] <- "biogeotracesSequence" 
head(merged)

#extracts the biogeotraces sample from the biogeotraces sequence ID
merged <- merged %>% mutate(biogeotracesSample = str_extract(biogeotracesSequence, "[a-zA-Z0-9]*"))
head(merged)

#the only biogeotraces samples that are NA are for the sean OTUs that weren't present 
#in biogeotraces
merged %>% filter(is.na(biogeotracesSample)) %>% distinct(biogeotracesSequence)
merged %>% filter(is.na(biogeotracesSample)) %>% distinct(OTU_ID) %>% nrow()

#there are multiple biogeotraces sequences in the biogeotraces samples which makes sense
merged %>% group_by(biogeotracesSample) %>% distinct(biogeotracesSequence) %>% summarize(n = n())
merged %>% group_by(biogeotracesSample) %>% distinct(biogeotracesSequence) %>% summarize(n = n()) %>% arrange(n)

head(seqG)

#these are the sean sequences that are present in multiple sean samples
seqG %>% group_by(sequence) %>% distinct(sample) %>% summarize(n = n()) %>% filter(n > 1)

head(merged)

#these are the sean sequences that are present in multiple sean samples
merged %>% group_by(sequence) %>% distinct(sample) %>% summarize(n = n()) %>% filter(n > 1)

#there are NAs when a sequence in a sean sample did not correspond to any of the biogeotraces sequences
head(merged)
merged %>% filter(is.na(biogeotracesSample)) %>% nrow()
merged %>% filter(is.na(biogeotracesSample)) %>% distinct(sequence) %>% nrow()
merged %>% filter(is.na(biogeotracesSample)) %>% distinct(OTU_ID) %>% nrow()
merged %>% filter(!is.na(biogeotracesSample)) %>% distinct(sequence) %>% nrow()
merged %>% filter(!is.na(biogeotracesSample)) %>% distinct(OTU_ID) %>% nrow()

#excludes rows in merged for OTUs/sequences in sean samples that were not found in biogeotraces
head(merged)
merged %>% filter(!is.na(biogeotracesSample))
merged_withAllSeanOTUs <- merged
merged <- merged %>% filter(!is.na(biogeotracesSample))

#all of the rows are distinct
merged %>% nrow()
merged %>% distinct() %>% nrow()

#these are the sequences and OTUs that are present in both sean and biogeotraces samples
merged %>% head()
merged %>% distinct(sequence) %>% nrow()
merged %>% distinct(OTU_ID) %>% nrow()

#these are the matches between OTUs for a given sean sample in which there are multiple 
#sean sequences present in the OTU in the sean sample
head(merged)
merged %>% group_by(sample, OTU_ID, biogeotracesSequence, biogeotracesSample) %>% summarize(n = n()) %>% arrange(desc(n))
merged %>% group_by(sample, OTU_ID, biogeotracesSequence, biogeotracesSample) %>% summarize(n = n()) %>% arrange(desc(n)) %>% filter(sample == "JW7")
merged %>% group_by(sample, OTU_ID, biogeotracesSequence, biogeotracesSample) %>% distinct(sequence) %>% summarize(n = n()) %>% arrange(desc(n))
merged %>% group_by(sample, OTU_ID, biogeotracesSequence, biogeotracesSample) %>% distinct(sequence) %>% summarize(n = n()) %>% arrange(desc(n)) %>% filter(sample == "JW7")

merged %>% filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf" & sample == "JW4" & biogeotracesSequence == "SRR5720219.17338980.2biogeotraces")
merged %>% filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf" & sample == "JW4")
merged %>% filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf" & sample == "JW4") %>% distinct(sequence)

merged %>% group_by(sample, OTU_ID, biogeotracesSequence, biogeotracesSample, sequence) %>% 
  summarize(n = n()) %>% filter(n != 1)

#these are the OTUs that have multiple sean sequences in a sean sample
head(merged)
merged %>% group_by(OTU_ID, sample) %>% distinct(sequence) %>% summarize(n = n()) %>% arrange(desc(n))

#there are not multiple sequence abundances in sean samples
head(merged)
merged %>% group_by(sample, sequence) %>% distinct(abundance) %>% summarize(n = n()) %>% filter(n > 1)

merged %>% group_by(sample, sequence, biogeotracesSequence) %>% summarize(n = n()) %>% filter(n > 1)

#these are the sean sequences that are in OTUs that are also in biogeotraces
#after excluding Pro and Syn
head(merged)
merged %>% ungroup() %>% distinct(sequence) %>% nrow()

#gets distinct matches between sean samples and biogeotraces sequences based on OTUs
#I do not care about sean sequences and their abundances
merged <- merged %>% ungroup() 
head(merged)
merged <- merged %>% distinct(sample, OTU_ID, biogeotracesSequence, biogeotracesSample)
head(merged)

merged %>% filter(is.na(biogeotracesSample))

#none of the biogeotraces sequences are present in multiple OTUs
head(merged)
merged %>% group_by(sample, biogeotracesSequence, biogeotracesSample) %>% summarize(n = n()) %>% filter(n != 1)
merged %>% group_by(biogeotracesSequence) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% filter(n != 1)


##loads metadata for biogeotraces samples
#from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP109831&o=acc_s%3Aa&s=SRR6507277,SRR6507278,SRR6507279,SRR6507280,SRR5720219,SRR5720220,SRR5720221,SRR5720222,SRR5720223,SRR5720224,SRR5720225,SRR5720226,SRR5720227,SRR5720228,SRR5720229,SRR5720230,SRR5720231,SRR5720232,SRR5720233,SRR5720234,SRR5720235,SRR5720236,SRR5720237,SRR5720238,SRR5720239,SRR5720240,SRR5720241,SRR5720242,SRR5720243,SRR5720244,SRR5720245,SRR5720246,SRR5720247,SRR5720248,SRR5720249,SRR5720250,SRR5720251,SRR5720252,SRR5720253,SRR5720254,SRR5720255,SRR5720256,SRR5720257,SRR5720258,SRR5720259,SRR5720260,SRR5720261,SRR5720262,SRR5720263,SRR5720264,SRR5720265,SRR5720266,SRR5720267,SRR5720268,SRR5720269,SRR5720270,SRR5720271,SRR5720272,SRR5720273,SRR5720274,SRR5720275,SRR5720276,SRR5720277,SRR5720278,SRR5720279,SRR5720280,SRR5720281,SRR5720282,SRR5720283,SRR5720284,SRR5720285,SRR5720286,SRR5720287,SRR5720288,SRR5720289,SRR5720290,SRR5720291,SRR5720292,SRR5720293,SRR5720294,SRR5720295,SRR5720296,SRR5720297,SRR5720298,SRR5720299,SRR5720300,SRR5720301,SRR5720302,SRR5720303,SRR5720304,SRR5720305,SRR5720306,SRR5720307,SRR5720308,SRR5720309,SRR5720310,SRR5720311,SRR5720312,SRR5720313,SRR5720314,SRR5720315,SRR5720316,SRR5720317,SRR5720318,SRR5720319,SRR5720320,SRR5720321,SRR5720322,SRR5720323,SRR5720324,SRR5720325,SRR5720326,SRR5720327,SRR5720328,SRR5720329,SRR5720330,SRR5720331,SRR5720332,SRR5720333,SRR5720334,SRR5720335,SRR5720336,SRR5720337,SRR5720338,SRR5720339,SRR5720340,SRR5720341,SRR5720342,SRR5720343,SRR5720344
meta1 <- read_csv("SraRunTable.csv")

#from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP110813
meta2 <- read_csv("SraRunTable2.csv")

#cleans up metadata bottle ID and collection date
meta1 %>% head()
meta1$bottle_id <- as.character(meta1$bottle_id)
meta1$collection_date <- dmy(meta1$collection_date)
meta1 %>% head()

#combines meta1 and meta2
meta <- bind_rows(meta1, meta2)

#there is only one row for each sample in the metadata
head(meta)
meta %>% group_by(Run) %>% summarize(n = n()) %>% filter(n > 1)

#all of the biogeotraces samples have corresponding metadata in meta
head(merged)
head(meta)
merged %>% anti_join(meta, by = c("biogeotracesSample" = "Run"))

merged %>% filter(is.na(biogeotracesSample))

#there are biogeotraces samples in meta that are not in merged
nrow(meta)
merged %>% distinct(biogeotracesSample) %>% nrow()

#there must not have been any matches to sean samples for these biogeotraces samples
meta %>% anti_join(merged, by = c("Run" = "biogeotracesSample")) %>% select(Run)

#adds metadata of each biogeotraces data to merged
nrow(merged)
merged <- merged %>% left_join(meta, by = c("biogeotracesSample" = "Run"))
nrow(merged)

merged %>% filter(is.na(lat_lon)) %>% nrow()

#makes depth variable numeric
str(merged$Depth)
merged$Depth <- parse_number(merged$Depth)
str(merged$Depth)

#gets rid of unnecessary variables
head(merged)
colnames(merged)
merged %>% distinct(cruise_station)
merged %>% distinct(geotraces_section)
toPlot <- merged %>% select(sample, OTU_ID, biogeotracesSample, collection_date, cruise_id, Depth, geo_loc_name, lat_lon)

#I don't care about how many sequences in each biogeotraces sample match to each sean sample 
#I just care about how many OTUs match
colnames(toPlot)
toPlot %>% nrow()
toPlot %>% distinct() %>% nrow()
toPlot <- toPlot %>% distinct()
toPlot %>% nrow()

#loads metadata for each sean culture
#from cleanUpEnvVariablesOfCulturesForPCA.R
cultureMeta <- read_csv("cleanedUpCultureEnvVariablesForPCA.csv")
nrow(cultureMeta)

#adds "culture_" to the beginning of each cultureMeta variable so that I can 
#distinguish from biogeotraces metadata
head(cultureMeta)
colnames(cultureMeta)
colnames(cultureMeta) <- str_c("culture_", colnames(cultureMeta))
colnames(cultureMeta)


#all of the sean culture samples have corresponding metadata 
head(toPlot)
head(cultureMeta)
toPlot %>% anti_join(cultureMeta, by = c("sample" = "culture_CULTURE_originalName")) %>% distinct(sample)

#adds metadata to sean culture samples
nrow(toPlot)
toPlot <- toPlot %>% left_join(cultureMeta, by = c("sample" = "culture_CULTURE_originalName"))
nrow(toPlot)

#there are NAs for the culture metadata in toPlot when there is not metadata available for that culture 
#and metadata variable
cultureMeta %>% filter(is.na(culture_ECOTYPE)) %>% distinct(culture_CULTURE) %>% arrange(culture_CULTURE)
toPlot %>% filter(is.na(culture_ECOTYPE)) %>% distinct(sample) %>% arrange(sample)

#there are not multiple matches between sean samples, OTUs, and biogeotraces samples which is good
toPlot %>% head()
toPlot %>% group_by(sample, OTU_ID, biogeotracesSample) %>% summarize(n = n()) %>% filter(n != 1)

head(toPlot)

toPlot %>% group_by(sample, cruise_id) %>% summarize(n = n())
toPlot %>% group_by(sample, cruise_id) %>% distinct(OTU_ID) %>% summarize(n = n())

#makes a variable for general cruise
toPlot %>% distinct(cruise_id)
toPlot <- toPlot %>% mutate(cruiseGeneral = str_extract(cruise_id, "BATS|HOT"))
toPlot <- toPlot %>% mutate(cruiseGeneral = ifelse(is.na(cruiseGeneral), cruise_id, cruiseGeneral))

toPlot %>% distinct(cruiseGeneral) %>% arrange(cruiseGeneral)

toPlot %>% group_by(sample, cruiseGeneral) %>% summarize(n = n())
toPlot %>% group_by(sample, cruiseGeneral) %>% distinct(OTU_ID) %>% summarize(n = n())

str(toPlot)

toPlot %>% select(Depth)
str(toPlot$Depth)

#arrange depths
depthOrder <- toPlot %>% arrange(Depth) %>% distinct(Depth)
depthOrder
depthOrder <- depthOrder %>% as.list()
depthOrder
str(depthOrder)
depthOrder <- lapply(depthOrder, as.character)
depthOrder

str(toPlot$Depth)
toPlot$factorDepth <- as.character(toPlot$Depth)
str(toPlot$factorDepth)

str(depthOrder$Depth)
toPlot$factorDepth <- factor(toPlot$factorDepth, levels = depthOrder$Depth)
toPlot$factorDepth
depthOrder$Depth

str(toPlot$factorDepth)

#toPlot %>% distinct(factorDepth) %>% View()

#https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
toPlot %>% select(Depth, factorDepth) %>% head()
toPlot %>% mutate(factorDepth = as.numeric(as.character(factorDepth))) %>% select(Depth, factorDepth) %>% filter(Depth != factorDepth)


#biomG has all of the biogeotraces and sean samples while merged has the sean samples
#and only the biogeotraces samples that have OTUs that are also in Sean samples
merged %>% distinct(biogeotracesSample) %>% nrow()
merged %>% filter(!str_detect(biogeotracesSample, "SRR")) %>% nrow()
merged %>% distinct(sample) %>% nrow()
biomG %>% head()
biomG %>% filter(!str_detect(sample, "SRR"))

#makes a dataframe with all of the biogeotraces OTUs
numOTUs_biogeotraces <- biomG

#renames variables
colnames(numOTUs_biogeotraces)
colnames(numOTUs_biogeotraces)[2] <- "sequence"

#gets rid of unnecessary variable
numOTUs_biogeotraces %>% head()
numOTUs_biogeotraces %>% distinct(present)
numOTUs_biogeotraces <- numOTUs_biogeotraces %>% select(-present)

#gets just the OTUs that are in biogeotraces samples
numOTUs_biogeotraces %>% filter(str_detect(sequence, "biogeotraces"))
numOTUs_biogeotraces <- numOTUs_biogeotraces %>% filter(str_detect(sequence, "biogeotraces"))

#makes a variable for the biogeotraces sample
numOTUs_biogeotraces %>% head()
numOTUs_biogeotraces <- numOTUs_biogeotraces %>% mutate(sample = str_extract(sequence, "[a-zA-Z0-9]*"))
numOTUs_biogeotraces %>% head()
numOTUs_biogeotraces %>% filter(is.na(sample))

#there are 610 biogeotraces samples which is correct
numOTUs_biogeotraces %>% distinct(sample) %>% nrow()

#this includes Pro and Syn OTUs
head(numOTUs_biogeotraces)
numOTUs_biogeotraces %>% distinct(OTU_ID) %>% nrow()

#makes dataframe for later calculating the mean relative abundance of each OTU across 
#depth ranges
avgRelAbun <- numOTUs_biogeotraces
head(avgRelAbun)

#gets the distinct OTUs in each sample
numOTUs_biogeotraces <- numOTUs_biogeotraces %>% distinct(sample, OTU_ID)

head(numOTUs_biogeotraces)

#merged only has the OTUs that were also present in biogeotraces samples
#no pro or syn
merged %>% head()
merged %>% filter(is.na(biogeotracesSequence)) %>% nrow()

#merged_withAllSeanOTUs has all of the non-Pro, non-Syn OTUs in the sean culture samples 
#whether or not they were also present in biogeotraces samples
head(merged_withAllSeanOTUs)
merged_withAllSeanOTUs %>% filter(is.na(biogeotracesSequence)) %>% nrow()
merged_withAllSeanOTUs %>% filter(is.na(biogeotracesSequence)) %>% distinct(OTU_ID) %>% nrow()
merged_withAllSeanOTUs %>% distinct(OTU_ID) %>% nrow()

#merged_withAllSeanOTUs has all of the OTUs that are present in Sean samples after 
#excluding Pro and Syn sequences
merged_withAllSeanOTUs %>% head()
merged_withAllSeanOTUs %>% distinct(sample) %>% nrow()
merged_withAllSeanOTUs %>% distinct(sequence) %>% nrow()
merged_withAllSeanOTUs %>% distinct(OTU_ID) %>% nrow()

merged_withAllSeanOTUs %>% filter(is.na(OTU_ID))

merged_withAllSeanOTUs %>% head()

#gets the distinct OTUs in each sean sample
numOTUs_sean <- merged_withAllSeanOTUs %>% distinct(sample, OTU_ID)

numOTUs_sean %>% head()
numOTUs_biogeotraces %>% head()
numOTUs_biogeotraces %>% filter(!str_detect(sample, "SRR"))

#correct number of samples
numOTUs_sean %>% distinct(sample) %>% nrow()
numOTUs_biogeotraces %>% distinct(sample) %>% nrow()

##correct number of OTUs
#doesn't include pro and syn
numOTUs_sean %>% distinct(OTU_ID) %>% nrow()
#includes pro and syn
numOTUs_biogeotraces %>% distinct(OTU_ID) %>% nrow()

#combines the number of OTUs in biogeotraces and sean culture samples dataframes
numOTUs <- bind_rows(numOTUs_biogeotraces, numOTUs_sean)

#loads taxonomies of OTUs
#from clusterSeanAndBiogeotracesSeqs.txt (checked)
tax <- read_tsv("sequencesForTree_noContamination_noMock_withProSyn_no1223_biogeotracesTaxonomy.tsv")

#excludes first row because it just has comments
tax[1,]
tax <- tax[-1,]
tax %>% head()

#gets rid of unnecessary variable
tax <- tax %>% select(1,2)
tax %>% head()

#the pro, syn OTU that I excluded from the sean samples must also be in biogeotraces
numOTUs %>% distinct(OTU_ID) %>% nrow()
tax %>% distinct(`Feature ID`) %>% nrow()

numOTUs %>% head()
tax %>% head()

#adds taxonomy of each OTU to numOTUs
nrow(numOTUs)
numOTUs <- numOTUs %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(numOTUs)

#adds taxonomy of each OTU to avgRelAbun
head(avgRelAbun)
nrow(avgRelAbun)
avgRelAbun <- avgRelAbun %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(avgRelAbun)

numOTUs %>% filter(is.na(Taxon))
avgRelAbun %>% filter(is.na(Taxon))

numOTUs %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(OTU_ID)
numOTUs %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)
numOTUs %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(Taxon)

#this makes sense because I already excluded the pro and syn sequences in sean samples
numOTUs %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% filter(!str_detect(sample, "SRR"))

avgRelAbun %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)
avgRelAbun %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(Taxon)

#excludes pro and syn OTUs from avgRelAbun
avgRelAbun <- avgRelAbun %>% filter(!str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc"))
avgRelAbun %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(Taxon)

#avgRelAbun only has biogeotraces samples and has all of them
head(avgRelAbun)
avgRelAbun %>% filter(!str_detect(sample, "SRR"))
avgRelAbun %>% distinct(sample) %>% nrow()

#gets the number of non-pro, non-syn sequences in each sample
head(avgRelAbun)
numSeqBySample <- avgRelAbun %>% group_by(sample) %>% distinct(sequence) %>% summarize(numSeqSample = n())
nrow(numSeqBySample)
head(numSeqBySample)

#calculates the number of sequences of each otu in each sample
head(avgRelAbun)
avgRelAbun <- avgRelAbun %>% group_by(OTU_ID, sample) %>% distinct(sequence) %>% summarize(numSeq = n())
head(avgRelAbun)

head(avgRelAbun)
head(numSeqBySample)

#adds total number of sequences in each biogeotraces sample to avgRelAbun
nrow(avgRelAbun)
avgRelAbun <- avgRelAbun %>% ungroup() %>% left_join(numSeqBySample, by = c("sample"))
nrow(avgRelAbun)

head(avgRelAbun)
avgRelAbun %>% filter(is.na(numSeqSample))

#calculates the relative abundnace of each OTU in each sample
avgRelAbun <- avgRelAbun %>% mutate(relAbun = numSeq/numSeqSample)
head(avgRelAbun)

avgRelAbun %>% group_by(sample) %>% summarize(totalProp = sum(relAbun)) %>% filter(totalProp != 1)

avgRelAbun %>% group_by(sample) %>% summarize(meanProp = mean(relAbun))

#excludes pro and syn OTUs from numOTUs
numOTUs %>% distinct(OTU_ID) %>% nrow()
numOTUs <- numOTUs %>% filter(!str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc"))
numOTUs %>% distinct(OTU_ID) %>% nrow()
numOTUs %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)

#gets rid of unnecessary variable
numOTUs %>% head()
numOTUs <- numOTUs %>% select(-Taxon)
numOTUs %>% head()

#calculates the number of OTUs by sample
#head(numOTUs)
#numOTUs <- numOTUs %>% group_by(sample) %>% distinct(OTU_ID) %>% summarize(numOTUBiogeotraces = n())
#head(numOTUs)

nrow(numOTUs)

#same number of rows which is good because this includes all biogeotraces samples 
#not just ones with OTUs also found in sean samples
numOTUs %>% filter(str_detect(sample, "SRR")) %>% distinct(sample) %>% nrow()
nrow(meta)

numOTUs %>% filter(!str_detect(sample, "SRR")) %>% distinct(sample) %>% nrow()

head(meta)
head(numOTUs)

#adds OTUs present in samples to metadata of biogeotraces samples so that every biogeotraces 
#sample is included not just ones that had OTUs that were also present in sean culture samples
nrow(meta)
meta <- meta %>% left_join(numOTUs, by = c("Run" = "sample"))
nrow(meta)

head(meta)
meta %>% filter(is.na(OTU_ID)) %>% nrow()

#makes variable for shortened cruise ID in meta
head(meta)
colnames(meta)
meta <- meta %>% mutate(cruiseGeneral = str_extract(cruise_id, "BATS|HOT"))
meta <- meta %>% mutate(cruiseGeneral = ifelse(is.na(cruiseGeneral), cruise_id, cruiseGeneral))
meta %>% distinct(cruiseGeneral) %>% arrange(cruiseGeneral)

#each biogeotraces sample is only present in meta once but has multiple OTUs
meta %>% nrow()
meta %>% distinct(Run) %>% summarize(n = n())

#calculates total number of non-Pro, non-Syn OTUs in biogeotraces samples grouped by short cruise ID
head(meta)
meta <- meta %>% ungroup()
cruiseGeneralSummary <- meta %>% group_by(cruiseGeneral) %>% distinct(OTU_ID) %>% summarize(sumOTUBiogeotraces = n())
cruiseGeneralSummary

toPlot %>% head()
toPlot %>% filter(is.na(biogeotracesSample)) %>% nrow()
toPlot %>% distinct(sample) %>% nrow()

toPlot %>% filter(is.na(OTU_ID)) %>% nrow()

#calculates the number of non-Pro, non-Syn OTUs in each sean culture and biogeotraces sample
numOTUs %>% head()
numOTUSummary <- numOTUs %>% group_by(sample) %>% distinct(OTU_ID) %>% summarize(numOTUInSample = n())
head(numOTUSummary)
nrow(numOTUSummary)

head(toPlot)
toPlot %>% group_by(sample, cruiseGeneral) %>% distinct(OTU_ID) %>% summarize(n = n())
toPlot %>% group_by(sample, cruiseGeneral) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% nrow()
toPlot %>% group_by(sample, cruiseGeneral) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(cruiseGeneralSummary, by = c("cruiseGeneral")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  nrow()
toPlot %>% group_by(sample, cruiseGeneral) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(cruiseGeneralSummary, by = c("cruiseGeneral")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  filter(is.na(sumOTUBiogeotraces)|is.na(numOTUInSample))
toPlot %>% group_by(sample, cruiseGeneral) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(cruiseGeneralSummary, by = c("cruiseGeneral")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  ungroup() %>% distinct(cruiseGeneral, sumOTUBiogeotraces)


#not all of the biogeotraces OTUs match up to sean culture samples which is why the sumProps 
#of each shortened cruise ID don't add up to 1 or more than 1
toPlot %>% group_by(sample, cruiseGeneral) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(cruiseGeneralSummary, by = c("cruiseGeneral")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>% ungroup() %>%
  mutate(prop = n/(sumOTUBiogeotraces+numOTUInSample-n)) %>% ungroup() %>% group_by(cruiseGeneral) %>% summarize(sumProp = sum(prop))

toPlot %>% distinct(cruiseGeneral)
toPlot %>% distinct(culture_CRUISE)

#makes a variable for shortened cruise ID for culture samples
toPlot <- toPlot %>% mutate(culture_cruiseGeneral = str_extract(culture_CRUISE, "BATS|HOT|EqPac"))
toPlot <- toPlot %>% mutate(culture_cruiseGeneral = ifelse(is.na(culture_cruiseGeneral), culture_CRUISE, culture_cruiseGeneral))

toPlot %>% distinct(culture_CRUISE, culture_cruiseGeneral)

#arranges sean culture samples by culture's shortened cruise ID
toPlot %>% arrange(culture_cruiseGeneral) %>% distinct(culture_cruiseGeneral)
order <- toPlot %>% arrange(culture_cruiseGeneral) %>% distinct(sample)
order
order <- order %>% select(sample)
order
order <- order %>% as.list()
order
order <- order$sample
order

toPlot$sample <- factor(toPlot$sample, levels = order)
levels(toPlot$sample)

#toPlot %>% group_by(sample, cruiseGeneral) %>% summarize(n = n()) %>% 
  #left_join(cruiseGeneralSummary, by = c("cruiseGeneral")) %>% 
  #left_join(numOTUSummary, by = c("sample")) %>%
  #ggplot() + geom_bar(aes(x = sample, y = n/(sumOTUBiogeotraces+numOTUInSample-n), fill = cruiseGeneral), stat = 'identity', position = 'dodge') + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  #labs(fill = "BioGeoTraces cruise", x = "")

#ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_cruise_OTU.png", dpi = 600, height = 6, width = 14)
#ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_cruise_OTU.svg", dpi = 600, height = 6, width = 14)


#toPlot %>% distinct(sample, culture_cruiseGeneral) %>% arrange(sample) %>% View()


#gets the total number of distinct OTUs across biogeotraces samples, grouped by location
meta %>% distinct(geo_loc_name)
geo_loc_nameSummary <- meta %>% group_by(geo_loc_name) %>% distinct(OTU_ID) %>% summarize(sumOTUBiogeotraces = n())
geo_loc_nameSummary

toPlot %>% group_by(sample, geo_loc_name) %>% distinct(OTU_ID) %>% summarize(n = n())
toPlot %>% group_by(sample, geo_loc_name) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% nrow() 
toPlot %>% group_by(sample, geo_loc_name) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(geo_loc_nameSummary, by = c("geo_loc_name")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  ungroup() %>% nrow()
toPlot %>% group_by(sample, geo_loc_name) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(geo_loc_nameSummary, by = c("geo_loc_name")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  ungroup() %>% filter(is.na(sumOTUBiogeotraces)|is.na(numOTUInSample))
toPlot %>% group_by(sample, geo_loc_name) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(geo_loc_nameSummary, by = c("geo_loc_name")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  ungroup() %>% distinct(geo_loc_name, sumOTUBiogeotraces)

toPlot %>% distinct(geo_loc_name)
toPlot %>% distinct(`culture_PLACE OF ORIGIN`)

#arranges biogeotraces locations alphabetically
toPlot$geo_loc_name <- factor(toPlot$geo_loc_name, levels = c("Atlantic Ocean", "Atlantic Ocean: Sargasso Sea\\, BATS", "Pacific Ocean", "Pacific Ocean: North Pacific Subtropical Gyre\\, Station ALOHA"))

toPlot %>% distinct(geo_loc_name)
levels(toPlot$geo_loc_name)

toPlot %>% distinct(`culture_PLACE OF ORIGIN`)

#puts sean culture samples in order by their place of origin
toPlot %>% arrange(`culture_PLACE OF ORIGIN`) %>% distinct(`culture_PLACE OF ORIGIN`)
order <- toPlot %>% arrange(`culture_PLACE OF ORIGIN`) %>% distinct(sample)
order
order <- order %>% select(sample)
order
order <- order %>% as.list()
order
order <- order$sample
order

toPlot$sample <- factor(toPlot$sample, levels = order)

toPlot %>% group_by(sample, geo_loc_name) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(geo_loc_nameSummary, by = c("geo_loc_name")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>% ungroup() %>%
  mutate(sample = factor(sample, levels = order)) %>% 
  mutate(geo_loc_name = str_replace(geo_loc_name, "\\\\", "")) %>%
  ggplot() + geom_bar(aes(x = sample, y = n/(sumOTUBiogeotraces+numOTUInSample-n), fill = geo_loc_name), stat = 'identity', position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", fill = "BioGeoTraces region", y = "Number of OTUs shared / total number of OTUs in culture and BioGeotraces samples")

ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_location_OTU.png", dpi = 600, height = 10, width = 30)
ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_location_OTU.svg", dpi = 600, height = 10, width = 30)


levels(toPlot$sample)
#toPlot %>% arrange(sample) %>% distinct(sample, `culture_PLACE OF ORIGIN`) %>% View()


#toPlot %>% distinct(sample, `culture_PLACE OF ORIGIN`) %>% arrange(sample) %>% View()

#there are a couple of ourlier depths
toPlot %>% ggplot(aes(x = Depth)) + geom_bar(color = 'red')

#toPlot %>% distinct(Depth) %>% View()
#toPlot %>% group_by(Depth) %>% summarize(n = n()) %>% View()

str(toPlot$Depth)

#makes a variable for depth range for the biogeotraces samples
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth <= 50, "0-50m", NA))
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth > 50 & Depth <= 100, "50-100m", depthRange))
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth > 100 & Depth <= 150, "100-150m", depthRange))
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth > 150 & Depth <= 200, "150-200m", depthRange))
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth > 200 & Depth <= 250, "200-250m", depthRange))
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth > 250 & Depth <= 300, "250-300m", depthRange))
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth > 300 & Depth <= 350, "300-350m", depthRange))
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth > 350 & Depth <= 400, "350-400m", depthRange))
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth > 400 & Depth <= 450, "400-450m", depthRange))
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth > 450 & Depth <= 500, "450-500m", depthRange))

toPlot %>% filter(is.na(depthRange)) %>% distinct(Depth) %>% arrange(Depth)

toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth >= 1000 & Depth <= 1100, "1000-1100m", depthRange))
toPlot <- toPlot %>% mutate(depthRange = ifelse(Depth >= 4550 & Depth <= 5650, "4550-5650m", depthRange))

toPlot %>% filter(is.na(depthRange)) %>% nrow()
toPlot %>% distinct(depthRange)
toPlot %>% group_by(depthRange) %>% distinct(biogeotracesSample) %>% summarize(n = n())

#makes depth of biogeotraces samples numeric in meta
str(meta$Depth)
meta$Depth <- parse_number(meta$Depth)
str(meta$Depth)
meta %>% filter(is.na(Depth)) %>% nrow()

#makes a variable for depth range for the biogeotraces samples in meta also
meta <- meta %>% mutate(depthRange = ifelse(Depth <= 50, "0-50m", NA))
meta <- meta %>% mutate(depthRange = ifelse(Depth > 50 & Depth <= 100, "50-100m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(Depth > 100 & Depth <= 150, "100-150m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(Depth > 150 & Depth <= 200, "150-200m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(Depth > 200 & Depth <= 250, "200-250m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(Depth > 250 & Depth <= 300, "250-300m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(Depth > 300 & Depth <= 350, "300-350m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(Depth > 350 & Depth <= 400, "350-400m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(Depth > 400 & Depth <= 450, "400-450m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(Depth > 450 & Depth <= 500, "450-500m", depthRange))

meta %>% filter(is.na(depthRange)) %>% distinct(Depth) %>% arrange(Depth)

meta <- meta %>% mutate(depthRange = ifelse(Depth >= 1000 & Depth <= 1100, "1000-1100m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(Depth >= 4550 & Depth <= 5650, "4550-5650m", depthRange))

meta %>% filter(is.na(depthRange)) %>% nrow()
meta %>% distinct(depthRange)

#calculates the total number of OTUs across biogeotraces samples, based on depth range
depthRangeSummary <- meta %>% group_by(depthRange) %>% distinct(OTU_ID) %>% summarize(sumOTUBiogeotraces = n())
depthRangeSummary

toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n())
toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% nrow()
toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  nrow()
toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  ungroup() %>% filter(is.na(sumOTUBiogeotraces)|is.na(numOTUInSample))
toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  ungroup() %>% distinct(depthRange, sumOTUBiogeotraces)

toPlot %>% distinct(depthRange)
#toPlot %>% distinct(culture_DEPTH) %>% arrange(culture_DEPTH) %>% View()

str(toPlot$culture_DEPTH)

#makes a variable for depth range for the sean culture samples
toPlot <- toPlot %>% mutate(culture_depthRange = ifelse(culture_DEPTH <= 50, "0-50m", NA))
toPlot <- toPlot %>% mutate(culture_depthRange = ifelse(culture_DEPTH > 50 & culture_DEPTH <= 100, "50-100m", culture_depthRange))
toPlot <- toPlot %>% mutate(culture_depthRange = ifelse(culture_DEPTH > 100 & culture_DEPTH <= 150, "100-150m", culture_depthRange))
toPlot <- toPlot %>% mutate(culture_depthRange = ifelse(culture_DEPTH > 150 & culture_DEPTH <= 200, "150-200m", culture_depthRange))

toPlot %>% filter(is.na(culture_depthRange)) %>% distinct(culture_DEPTH)

#toPlot %>% distinct(culture_DEPTH, culture_depthRange) %>% arrange(culture_DEPTH) %>% View()

toPlot %>% distinct(culture_depthRange)

#puts culture depth ranges in order
#is this okay to do???
toPlot$culture_depthRange <- factor(toPlot$culture_depthRange, levels = c("0-50m", "50-100m", "100-150m", "150-200m", NA))

levels(toPlot$culture_depthRange)

toPlot %>% distinct(culture_depthRange) %>% arrange(culture_depthRange)

##puts sean culture samples in order by culture depth range
toPlot %>% arrange(culture_depthRange) %>% distinct(culture_depthRange)
order <- toPlot %>% arrange(culture_depthRange) %>% distinct(sample)
order
order <- order %>% select(sample)
order
order <- order %>% as.list()
order
order <- order$sample
order

length(order)
toPlot %>% distinct(sample) %>% nrow()

toPlot$sample <- factor(toPlot$sample, levels = order)

toPlot %>% distinct(sample) %>% nrow()
toPlot %>% ungroup() %>% distinct(depthRange) %>% arrange(depthRange)
toPlot %>% filter(is.na(depthRange)) %>% nrow()
toPlot %>% ungroup() %>% distinct(depthRange) %>% nrow()
depthRangeOrder <- c("0-50m", "50-100m", "100-150m", "150-200m", "200-250m", "250-300m", "300-350m", "350-400m", "400-450m", "450-500m", "1000-1100m", "4550-5650m")
depthRangeOrder %>% length()
toPlot %>% ungroup() %>% distinct(depthRange) %>% filter(depthRange %in% depthRangeOrder) %>% nrow()

toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>% ungroup() %>%
  mutate(sample = factor(sample, levels = order)) %>% 
  mutate(depthRange = factor(depthRange, levels = depthRangeOrder)) %>%
  ggplot() + geom_bar(aes(x = sample, y = n/(sumOTUBiogeotraces+numOTUInSample-n), fill = depthRange), stat = 'identity', position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", fill = "BioGeoTraces depth", y = "Number of OTUs shared / total number of OTUs in culture and BioGeotraces samples")

ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU.png", dpi = 600, height = 10, width = 30)
ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU.svg", dpi = 600, height = 10, width = 30)

toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% nrow()
toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>% nrow()
toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>% filter(is.na(sumOTUBiogeotraces)) %>% nrow()
toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>% filter(is.na(numOTUInSample)) %>% nrow()


grouped <- toPlot %>% group_by(sample, depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>% ungroup() %>%
  mutate(prop = n/(sumOTUBiogeotraces+numOTUInSample-n))

#not all of the cultures are present in each depth range
head(grouped)
grouped %>% group_by(depthRange) %>% summarize(n = n())

#gets rid of unnecessary variables so I can spread data
head(grouped)
grouped <- grouped %>% select(sample, depthRange, prop)
head(grouped)

#spread data so that if a culture is not present in a depth range, it
#has a prop of 0
grouped_s <- grouped %>% spread(key = sample, value = prop, fill = 0)
dim(grouped_s)

#gathers data
colnames(grouped_s)
grouped <- grouped_s %>% gather(`0601`:`WH8109`, key = "sample", value = "prop")

#now there is an entry for ecah culture for each depth range
head(grouped)
grouped %>% group_by(depthRange) %>% summarize(n = n()) %>% filter(n != 74)

str(grouped)

grouped %>% 
  group_by(depthRange) %>% summarize(meanProp = mean(prop), stdev = sd(prop)) %>% 
  mutate(depthRange = factor(depthRange, levels = depthRangeOrder)) %>% 
  ggplot() + geom_bar(aes(x = depthRange, y=meanProp), stat = 'identity') + 
  geom_errorbar(aes(x = depthRange, ymin=meanProp-stdev, ymax=meanProp+stdev), stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "BioGeoTraces depth", y = "Mean proportion of OTUs shared between culture and BioGeotraces samples") + 
  theme(text = element_text(size=14))

ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU_grouped.png", dpi = 600, height = 8, width = 10)
ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU_grouped.svg", dpi = 600, height = 8, width = 10)

head(grouped)
totalPropBySample <- grouped %>% group_by(sample) %>% summarize(totalSampleProp = sum(prop))

head(grouped)
head(totalPropBySample)

nrow(grouped)
grouped <- grouped %>% left_join(totalPropBySample, by = c("sample"))
nrow(grouped)

grouped %>% filter(is.na(totalSampleProp))

head(grouped)

grouped <- grouped %>% mutate(normProp = prop/totalSampleProp)

grouped %>% group_by(sample) %>% summarize(totalProp = sum(normProp)) %>% filter(totalProp != 1)

grouped %>% 
  group_by(depthRange) %>% summarize(meanProp = mean(normProp), stdev = sd(normProp)) %>% 
  mutate(depthRange = factor(depthRange, levels = depthRangeOrder)) %>% 
  ggplot() + geom_bar(aes(x = depthRange, y=meanProp), stat = 'identity') + 
  geom_errorbar(aes(x = depthRange, ymin=meanProp-stdev, ymax=meanProp+stdev), stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "BioGeoTraces depth", y = "Mean proportion of OTUs shared between culture and BioGeotraces samples") + 
  theme(text = element_text(size=14))

head(toPlot)
toPlot %>% filter(is.na(sample)) %>% nrow()
toPlot %>% filter(is.na(OTU_ID)) %>% nrow()
toPlot %>% filter(is.na(biogeotracesSample)) %>% nrow()

toPlot %>% group_by(depthRange) %>% distinct(OTU_ID) %>% summarize(n = n())
toPlot %>% group_by(depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange"))

depthRangeSummary

toPlot %>% group_by(depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% ungroup() %>%
  mutate(depthRange = factor(depthRange, levels = depthRangeOrder)) %>%
  ggplot() + geom_bar(aes(x = depthRange, y = n/(sumOTUBiogeotraces)), stat = 'identity', position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", fill = "BioGeoTraces depth", y = "Proportion of OTUs in BioGEOTRACES samples also present in cultures") + 
  theme(text = element_text(size=14))

ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU_allCulturesGrouped.png", dpi = 600, height = 8, width = 10)
ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU_allCulturesGrouped.svg", dpi = 600, height = 8, width = 10)


head(avgRelAbun)
head(meta)

#gets the depthRange of each run
meta_depthRange <- meta %>% distinct(Run, Depth, depthRange)
head(meta_depthRange)
nrow(meta_depthRange)

head(avgRelAbun)
head(meta_depthRange)

#adds depthRange to the relative abundance of each non-pro, non-syn otu in 
#each biogeotraces sample
nrow(avgRelAbun)
avgRelAbun <- avgRelAbun %>% left_join(meta_depthRange, by = c("sample" = "Run"))
nrow(avgRelAbun)
head(avgRelAbun)

avgRelAbun %>% distinct(OTU_ID) %>% nrow()
avgRelAbun %>% distinct(sample) %>% nrow()

#calculates mean relative abundance of non-pro, non-syn otus by depth range
head(avgRelAbun)
avgRelAbun_depth <- avgRelAbun %>% group_by(depthRange) %>% summarize(avgRelAbun = mean(relAbun))
head(avgRelAbun)

toPlot %>% group_by(depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(avgRelAbun_depth, by = c("depthRange")) %>% ungroup()

toPlot %>% group_by(depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(avgRelAbun_depth, by = c("depthRange")) %>% ungroup() %>% filter(is.na(avgRelAbun))

toPlot %>% group_by(depthRange) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary, by = c("depthRange")) %>% 
  left_join(avgRelAbun_depth, by = c("depthRange")) %>% ungroup() %>%
  mutate(depthRange = factor(depthRange, levels = depthRangeOrder)) %>%
  ggplot() + geom_bar(aes(x = depthRange, y = (n/(sumOTUBiogeotraces))*avgRelAbun), stat = 'identity', position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", fill = "BioGeoTraces depth", y = "Proportion of OTUs in BioGEOTRACES samples also present in cultures\nnormalized by mean relative abundance of BioGEOTRACES OTUs") + 
  theme(text = element_text(size=14))

ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU_allCulturesGrouped_normalizedRelAbun.png", dpi = 600, height = 8, width = 10)
ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU_allCulturesGrouped_normalizedRelAbun.svg", dpi = 600, height = 8, width = 10)



#depthRangeSummary
#meta %>% group_by(depthRange) %>% distinct(OTU_ID) %>% summarize(sumOTUBiogeotraces = n())

#calculates the total number of OTUs across biogeotraces samples, based on depth range
str(meta$Depth)
meta %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% distinct(depthRange_broad)

#calculates the total number of OTUs in 0-200m and deeper than 200m
depthRangeSummary_broad <- meta %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% 
  group_by(depthRange_broad) %>% distinct(OTU_ID) %>% summarize(sumOTUBiogeotraces = n())
depthRangeSummary_broad

str(toPlot$Depth)
toPlot %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% distinct(depthRange_broad)
toPlot %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% 
  group_by(depthRange_broad) %>% distinct(OTU_ID) %>% summarize(n = n())

toPlot %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% 
  group_by(depthRange_broad) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary_broad, by = c("depthRange_broad")) %>% ungroup()

#makes a variable for broad depth range in avgRelAbun
head(avgRelAbun)
avgRelAbun <- avgRelAbun %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) 
avgRelAbun %>% distinct(depthRange_broad)
head(avgRelAbun)

#calculates mean relative abundance of non-pro, non-syn otus by broad depth range
head(avgRelAbun)
avgRelAbun_depth <- avgRelAbun %>% group_by(depthRange_broad) %>% summarize(avgRelAbun = mean(relAbun))
head(avgRelAbun_depth)

toPlot %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% 
  group_by(depthRange_broad) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary_broad, by = c("depthRange_broad")) %>% 
  left_join(avgRelAbun_depth, by = c("depthRange_broad")) %>% ungroup()

toPlot %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% 
  group_by(depthRange_broad) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary_broad, by = c("depthRange_broad")) %>% 
  left_join(avgRelAbun_depth, by = c("depthRange_broad")) %>% ungroup() %>% filter(is.na(avgRelAbun))

toPlot %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% 
  group_by(depthRange_broad) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary_broad, by = c("depthRange_broad")) %>% ungroup() %>%
  ggplot() + geom_bar(aes(x = depthRange_broad, y = n/(sumOTUBiogeotraces)), stat = 'identity', position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", fill = "BioGeoTraces depth", y = "Proportion of OTUs in BioGEOTRACES samples also present in cultures") + 
  theme(text = element_text(size=14)) 

ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU_allCulturesGrouped_200mOrDeeper.png", dpi = 600, height = 8, width = 10)
ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU_allCulturesGrouped_200mOrDeeper.svg", dpi = 600, height = 8, width = 10)



toPlot %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% 
  group_by(depthRange_broad) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary_broad, by = c("depthRange_broad")) %>% ungroup() %>%
  mutate(`Proportion of OTUs in BioGEOTRACES samples also present in cultures` = n/(sumOTUBiogeotraces)) %>% 
  select(depthRange_broad, `Proportion of OTUs in BioGEOTRACES samples also present in cultures`)


forFisher <- toPlot %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% 
  group_by(depthRange_broad) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary_broad, by = c("depthRange_broad")) %>% ungroup() 

forFisher

matrix(c(forFisher$n[1], forFisher$sumOTUBiogeotraces[1]-forFisher$n[1], 
         forFisher$n[2], forFisher$sumOTUBiogeotraces[2]-forFisher$n[2]),2,2)

fisher.test(matrix(c(forFisher$n[1], forFisher$sumOTUBiogeotraces[1]-forFisher$n[1], 
                     forFisher$n[2], forFisher$sumOTUBiogeotraces[2]-forFisher$n[2]),2,2), alternative = "two.sided")

toPlot %>% mutate(depthRange_broad = ifelse(Depth <=200, "0-200m", "Deeper than 200m")) %>% 
  group_by(depthRange_broad) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(depthRangeSummary_broad, by = c("depthRange_broad")) %>% 
  left_join(avgRelAbun_depth, by = c("depthRange_broad")) %>% ungroup() %>%
  ggplot() + geom_bar(aes(x = depthRange_broad, y = (n/(sumOTUBiogeotraces))*avgRelAbun), stat = 'identity', position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", fill = "BioGeoTraces depth", y = "Proportion of OTUs in BioGEOTRACES samples also present in cultures\nnormalized by mean relative abundance of BioGEOTRACES OTUs") + 
  theme(text = element_text(size=14)) 

ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU_allCulturesGrouped_200mOrDeeper_normalizedRelAbun.png", dpi = 600, height = 8, width = 10)
ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_depth_OTU_allCulturesGrouped_200mOrDeeper_normalizedRelAbun.svg", dpi = 600, height = 8, width = 10)



#toPlot %>% arrange(sample) %>% distinct(sample, culture_depthRange) %>% View()


levels(toPlot$sample)
#toPlot %>% distinct(sample, culture_depthRange) %>% arrange(sample) %>% View()

#extracts year from collection date of biogeotraces samples
toPlot %>% distinct(collection_date)
str(toPlot$collection_date)
toPlot <- toPlot %>% mutate(year = str_extract(collection_date, "[0-9]*"))
toPlot %>% distinct(year)
toPlot %>% filter(is.na(year)) %>% nrow()

#also extracts year from collection date of biogeotraces samples for meta
meta %>% distinct(collection_date)
str(meta$collection_date)
meta <- meta %>% mutate(year = str_extract(collection_date, "[0-9]*"))
meta %>% distinct(year)
toPlot %>% filter(is.na(year)) %>% nrow()

toPlot %>% distinct(year) %>% arrange(year)
meta %>% distinct(year) %>% arrange(year)

#calculates the total number of sequences across biogeotraces samples, grouped by year
meta
yearSummary <- meta %>% group_by(year) %>% distinct(OTU_ID) %>% summarize(sumOTUBiogeotraces = n())
yearSummary

toPlot %>% group_by(sample, year) %>% distinct(OTU_ID) %>% summarize(n = n())
toPlot %>% group_by(sample, year) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% nrow()
toPlot %>% group_by(sample, year) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(yearSummary, by = c("year")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  nrow()
toPlot %>% group_by(sample, year) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(yearSummary, by = c("year")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  ungroup() %>% filter(is.na(sumOTUBiogeotraces)|is.na(numOTUInSample))
toPlot %>% group_by(sample, year) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(yearSummary, by = c("year")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  ungroup() %>% distinct(year, sumOTUBiogeotraces)

colnames(toPlot)

#makes variable for year of isolation for sean culture samples
toPlot %>% distinct(`culture_DATE ISOLATED_Edited`)
toPlot <- toPlot %>% mutate(culture_year = str_extract(`culture_DATE ISOLATED_Edited`, "[0-9]*"))

toPlot %>% distinct(`culture_DATE ISOLATED_Edited`, culture_year)
toPlot %>% distinct(culture_year)

#the only NA years also have NA date isolated
toPlot %>% filter(is.na(culture_year)) %>% distinct(`culture_DATE ISOLATED_Edited`)

toPlot %>% distinct(culture_year) %>% arrange(culture_year)

#puts the sean culture samples in order by their year of isolation
toPlot %>% arrange(culture_year) %>% distinct(sample)
order <- toPlot %>% arrange(culture_year) %>% distinct(sample)
order
order <- order %>% select(sample)
order
order <- order %>% as.list()
order
order <- order$sample
order
length(order)

toPlot$sample <- factor(toPlot$sample, levels = order)

toPlot %>% group_by(sample, year) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(yearSummary, by = c("year")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>% ungroup() %>%
  mutate(sample = factor(sample, levels = order)) %>%
  ggplot() + geom_bar(aes(x = sample, y = n/(sumOTUBiogeotraces+numOTUInSample-n), fill = year), stat = 'identity', position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", fill = "BioGeoTraces year", y = "Number of OTUs shared / total number of OTUs in culture and BioGeotraces samples")

ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_year_OTU.png", dpi = 600, height = 10, width = 30)
ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_year_OTU.svg", dpi = 600, height = 10, width = 30)


toPlot %>% arrange(sample) %>% distinct(sample, culture_year) %>% View()


#toPlot %>% distinct(sample, culture_year) %>% arrange(sample) %>% View()


#makes a variable for whether biogeotraces sample is in Pacific or Atlantic
meta %>% distinct(geo_loc_name)
meta <- meta %>% mutate(altOrPac = str_extract(geo_loc_name, "Pacific Ocean|Atlantic Ocean"))
meta %>% distinct(altOrPac)
meta %>% distinct(geo_loc_name, altOrPac)

#calculates the total number of sequences across biogeotraces samples, based on whether biogeotraces 
#is in Pacific or Atlantic
meta <- meta %>% ungroup()
meta %>% head()
altPacSummary <- meta %>% group_by(altOrPac) %>% distinct(OTU_ID) %>% summarize(sumOTUBiogeotraces = n())
altPacSummary

#makes a variable for whether biogeotraces sample is in Pacific or Atlantic for toPlot also
toPlot %>% distinct(geo_loc_name)
toPlot <- toPlot %>% mutate(altOrPac = str_extract(geo_loc_name, "Pacific Ocean|Atlantic Ocean"))
toPlot %>% distinct(altOrPac)
toPlot %>% distinct(geo_loc_name, altOrPac)


toPlot %>% group_by(sample, altOrPac) %>% distinct(OTU_ID) %>% summarize(n = n())
toPlot %>% group_by(sample, altOrPac) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% nrow()
toPlot %>% group_by(sample, altOrPac) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(altPacSummary, by = c("altOrPac")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  nrow()
toPlot %>% group_by(sample, altOrPac) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(altPacSummary, by = c("altOrPac")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  ungroup() %>% filter(is.na(sumOTUBiogeotraces)|is.na(numOTUInSample))
toPlot %>% group_by(sample, altOrPac) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(altPacSummary, by = c("altOrPac")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>%
  ungroup() %>% distinct(altOrPac, sumOTUBiogeotraces)


###need to check from here down
toPlot %>% distinct(`culture_PLACE OF ORIGIN`)

#makes a variable for whether sean culture samples is in Pacific or Atlantic
toPlot <- toPlot %>% mutate(culture_altOrPac = str_extract(`culture_PLACE OF ORIGIN`, "Atlantic|Pacific"))
toPlot %>% distinct(`culture_PLACE OF ORIGIN`, culture_altOrPac)

toPlot <- toPlot %>% mutate(culture_altOrPac = ifelse(`culture_PLACE OF ORIGIN` == "Gulf Stream", "Atlantic", culture_altOrPac))
toPlot %>% distinct(`culture_PLACE OF ORIGIN`, culture_altOrPac)

toPlot <- toPlot %>% mutate(culture_altOrPac = ifelse(`culture_PLACE OF ORIGIN` == "Sargasso Sea", "Atlantic", culture_altOrPac))
toPlot %>% distinct(`culture_PLACE OF ORIGIN`, culture_altOrPac)

toPlot <- toPlot %>% mutate(culture_altOrPac = ifelse(`culture_PLACE OF ORIGIN` == "Arabian Sea", "Arabian Sea", culture_altOrPac))
toPlot %>% distinct(`culture_PLACE OF ORIGIN`, culture_altOrPac)

toPlot <- toPlot %>% mutate(culture_altOrPac = ifelse(`culture_PLACE OF ORIGIN` == "Mediterranean Sea", "Mediterranean Sea", culture_altOrPac))
toPlot %>% distinct(`culture_PLACE OF ORIGIN`, culture_altOrPac)

toPlot <- toPlot %>% mutate(culture_altOrPac = ifelse(`culture_PLACE OF ORIGIN` == "Red Sea", "Red Sea", culture_altOrPac))
toPlot %>% distinct(`culture_PLACE OF ORIGIN`, culture_altOrPac)

toPlot %>% distinct(altOrPac)
toPlot %>% distinct(culture_altOrPac)

#arranges toPlot by altOrPac
toPlot <- toPlot %>% arrange(altOrPac)

toPlot %>% distinct(altOrPac)
toPlot %>% distinct(culture_altOrPac)

#arranges sean culture samples by where they come from
toPlot %>% arrange(culture_altOrPac) %>% distinct(culture_altOrPac)
toPlot %>% arrange(culture_altOrPac) %>% distinct(sample)
order <- toPlot %>% arrange(culture_altOrPac) %>% distinct(sample)
order
order <- order %>% select(sample)
order
order <- order %>% as.list()
order
order <- order$sample
order

length(order)
toPlot %>% distinct(sample) %>% nrow()

toPlot$sample <- factor(toPlot$sample, levels = order)

toPlot

toPlot %>% group_by(sample, altOrPac) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% 
  left_join(altPacSummary, by = c("altOrPac")) %>% 
  left_join(numOTUSummary, by = c("sample")) %>% ungroup() %>% 
  mutate(sample = factor(sample, levels = order)) %>% 
  ggplot() + geom_bar(aes(x = sample, y = n/(sumOTUBiogeotraces+numOTUInSample-n), fill = altOrPac), stat = 'identity', position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", fill = "BioGeoTraces location", y = "Number of OTUs shared / total number of OTUs in culture and BioGeotraces samples")

ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_atlanticOrPacific_OTU.png", dpi = 600, height = 10, width = 30)
ggsave("seanSeqsNoProSyn_ComparedToBiogeotraces_atlanticOrPacific_OTU.svg", dpi = 600, height = 10, width = 30)

levels(toPlot$sample)
toPlot %>% arrange(sample) %>% distinct(sample, culture_altOrPac) %>% View()


###is this max depth correct???? 
###why are so many OTUs in cultures matching up to TAN1109 
###compared to other cruises???
toPlot %>% summarize(min = min(Depth), max = max(Depth))
#toPlot %>% filter(Depth == 5601) %>% View()
toPlot %>% filter(cruise_id == "TAN1109")
toPlot %>% filter(cruise_id == "TAN1109") %>% distinct(OTU_ID)
toPlot %>% filter(cruise_id == "TAN1109") %>% distinct(sample)




