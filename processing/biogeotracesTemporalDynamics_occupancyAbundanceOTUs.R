library(tidyverse)
library(lubridate)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

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

#loads the abundance of each sean sequence in each culture
#this includes pro and syn 
#sequences with prop < .002 have been excluded
#from heatmap_noContamination_noLowAbundance_noMock_noProSyn_without1223_propReads.R
seqG <- read_csv("forSeanComparisonToBiogeotraces.csv")

#this includes Pro and Syn sequences and the JW3 and 1223 samples
head(seqG)
seqG %>% distinct(sequence) %>% nrow()

seqG %>% filter(sample == "JW3")
seqG %>% filter(sample == "1223")

#excludes JW3 and 1223
seqG <- seqG %>% filter(sample != "JW3")
seqG <- seqG %>% filter(sample != "1223")

#includes Pro and Syn
seqG %>% distinct(sequence) %>% nrow()
seqG %>% distinct(sample) %>% nrow()
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

head(seqsNoProSyn)

#includes Pro and Syn
seqG %>% distinct(sequence) %>% nrow()

#doesn't include Pro and Syn 
seqsNoProSyn %>% distinct(V1) %>% nrow()

head(seqG)
head(seqsNoProSyn)

#excludes Pro and Syn sequences from seqG
seqG <- seqG %>% semi_join(seqsNoProSyn, by = c("sequence" = "V1"))

#correct number of sequences and samples
head(seqG)
seqG %>% distinct(sequence) %>% nrow()
seqG %>% distinct(sample) %>% nrow()

#biomG includes Pro and Syn sequences in sean samples
head(biomG)
biomG %>% filter(str_detect(sample, "sean"))
biomG %>% filter(str_detect(sample, "sean")) %>% nrow()

##makes the sean sequence IDs in biomG match seqG
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
biomG %>% filter(str_detect(sample, "contig")) %>% head()
seqG %>% head()
biomG %>% filter(str_detect(sample, "contig")) %>% left_join(seqG, by = c("sample" = "sequence")) %>% filter(is.na(abundance)) %>% nrow()


#these are the clusters arranged by how many sequences they have across 
#biogeotraces and sean samples
head(biomG)
biomG %>% group_by(OTU_ID) %>% summarize(n = n()) %>% arrange(desc(n))

#these are the clusters in sean samples that have multiple sequences across 
#sean samples
biomG %>% filter(str_detect(sample, "contig")) %>% head()
biomG %>% filter(str_detect(sample, "contig")) %>% group_by(OTU_ID) %>% summarize(n = n()) %>% arrange(desc(n))

#gets just the clusters that sean sequences are in 
head(biomG)
seanBiom <- biomG %>% filter(str_detect(sample, "contig"))

#gets just the clusters that biogeotraces sequences are in
biogeoBiom <- biomG %>% filter(!str_detect(sample, "contig"))

nrow(biomG)
nrow(seanBiom)
nrow(biogeoBiom)
nrow(seanBiom)+nrow(biogeoBiom)==nrow(biomG)

#doesn't include pro and syn
seqG %>% distinct(sequence) %>% nrow()
#includes pro and syn
seanBiom %>% distinct(sample) %>% nrow()

head(seqG)
head(seanBiom)

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
merged %>% distinct(sequence) %>% nrow()

head(merged)

head(biogeoBiom)
biogeoBiom %>% distinct(present)

#gets rid of unnecessary variable
biogeoBiom <- biogeoBiom %>% select(-present)

#there are 28008 OTUs present in biogeotraces
biogeoBiom %>% distinct(OTU_ID) %>% nrow()

head(biogeoBiom)
head(merged)

#there are 27925 OTUs present in biogeotraces that are not present in sean samples
biogeoBiom %>% anti_join(merged, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

#there are 83 OTUs present in biogeotraces that are also present in sean samples
biogeoBiom %>% semi_join(merged, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

27925+83

#there are 97 OTUs present in sean samples
#after excluding pro and syn
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
merged %>% filter(is.na(sample.y)) %>% nrow()
merged %>% filter(is.na(sample.y)) %>% distinct(OTU_ID) %>% nrow()

merged %>% distinct(OTU_ID) %>% nrow()

head(merged)

#all of the samples have at least one OTU that is in one of the biogeotraces samples
merged %>% distinct(sample.x) %>% nrow()
merged %>% filter(!is.na(sample.y)) %>% distinct(sample.x) %>% nrow()

#renames merged column names
head(merged)
colnames(merged)
colnames(merged)[1]
colnames(merged)[1] <- "sample" 
colnames(merged)[5]
colnames(merged)[5] <- "biogeotracesSequence" 
head(merged)

#extracts the biogeotraces sample from the biogeotraces sequence ID
merged <- merged %>% mutate(biogeotracesSample = str_extract(biogeotracesSequence, "[a-zA-Z0-9]*"))
head(merged)

#the only biogeotraces samples that are NA are for the sean OTUs that weren't present 
#in biogeotraces
merged %>% filter(is.na(biogeotracesSample)) %>% distinct(biogeotracesSequence)
merged %>% filter(is.na(biogeotracesSample)) %>% distinct(OTU_ID) %>% nrow()

#there are multiple biogeotraces sequences in many of the biogeotraces samples which makes sense
head(merged)
merged %>% group_by(biogeotracesSample) %>% distinct(biogeotracesSequence) %>% summarize(n = n())
merged %>% group_by(biogeotracesSample) %>% distinct(biogeotracesSequence) %>% summarize(n = n()) %>% arrange(n)

head(seqG)

#these are the sean sequences that are present in multiple sean samples
seqG %>% group_by(sequence) %>% distinct(sample) %>% summarize(n = n()) %>% filter(n > 1)

head(merged)

#these are the sean sequences that are present in multiple sean samples
merged %>% group_by(sequence) %>% distinct(sample) %>% summarize(n = n()) %>% filter(n > 1)

#there are NAs when an OTU in a sean sample was not found in any of the biogeotraces 
#samples
head(merged)
merged %>% filter(is.na(biogeotracesSample)) %>% nrow()
merged %>% filter(is.na(biogeotracesSample)) %>% distinct(sequence) %>% nrow()
merged %>% filter(is.na(biogeotracesSample)) %>% distinct(OTU_ID) %>% nrow()
merged %>% filter(!is.na(biogeotracesSample)) %>% distinct(sequence) %>% nrow()
merged %>% filter(!is.na(biogeotracesSample)) %>% distinct(OTU_ID) %>% nrow()

#I don't want to exclude rows in merged for OTUs/sequences in sean samples 
#that were not found in biogeotraces
head(merged)
merged %>% filter(!is.na(biogeotracesSample))
#merged_withAllSeanOTUs <- merged
#merged <- merged %>% filter(!is.na(biogeotracesSample))

merged %>% nrow()
merged %>% distinct() %>% nrow()

#these are the sequences and OTUs that are present in sean samples and some 
#are also present in biogeotraces samples
#doesn't include pro and syn
merged %>% head()
merged %>% distinct(sequence) %>% nrow()
merged %>% distinct(OTU_ID) %>% nrow()

#these are the matches between OTUs for a given sean sample in which there are multiple 
#sean sequences present in the OTU in the sean sample
head(merged)
merged %>% group_by(sample, OTU_ID, biogeotracesSequence, biogeotracesSample) %>% summarize(n = n()) %>% arrange(desc(n))
merged %>% group_by(sample, OTU_ID, biogeotracesSequence, biogeotracesSample) %>% summarize(n = n()) %>% arrange(desc(n)) %>% filter(sample == "JW7")
merged %>% group_by(sample, OTU_ID, biogeotracesSequence, biogeotracesSample) %>% distinct(sequence) %>% summarize(n = n()) %>% arrange(desc(n))

merged %>% filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf" & sample == "JW4" & biogeotracesSequence == "SRR5720219.17338980.2biogeotraces")
merged %>% filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf" & sample == "JW4") %>% distinct(sequence)

#there are the OTUs that have multiple sean sequences in a sean sample
head(merged)
merged %>% group_by(OTU_ID, sample) %>% distinct(sequence) %>% summarize(n = n()) %>% arrange(desc(n))

#there are not multiple sequence abundances in sean samples
head(merged)
merged %>% group_by(sample, sequence) %>% distinct(abundance) %>% summarize(n = n()) %>% arrange(desc(n)) %>% filter(n > 1)

#these are the sean sequences
#after excluding Pro and Syn
head(merged)
merged %>% ungroup() %>% distinct(sequence) %>% nrow()

#gets distinct matches between sean samples and biogeotraces sequences based on OTUs
#I do not care about sean sequences and their abundances
head(merged)
merged <- merged %>% ungroup() 
merged %>% filter(is.na(sequence))
merged <- merged %>% distinct(sample, OTU_ID, biogeotracesSequence, biogeotracesSample)
merged

#these are the OTUs that are present in sean samples that are not also 
#in biogeotraces samples
merged %>% filter(is.na(biogeotracesSample))
merged %>% filter(is.na(biogeotracesSample)) %>% distinct(OTU_ID) %>% nrow()

#none of the biogeotraces sequences are present in multiple OTUs unless it is for 
#an OTU that was not found in any biogeotraces samples
head(merged)
merged %>% group_by(sample, biogeotracesSequence, biogeotracesSample) %>% summarize(n = n()) %>% 
  filter(n != 1) %>% ungroup() %>% distinct(biogeotracesSequence)
head(merged)
merged %>% group_by(biogeotracesSequence) %>% distinct(OTU_ID) %>% summarize(n = n()) %>% filter(n != 1)


##loads metadata for biogeotraces samples
#from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP109831&o=acc_s%3Aa&s=SRR6507277,SRR6507278,SRR6507279,SRR6507280,SRR5720219,SRR5720220,SRR5720221,SRR5720222,SRR5720223,SRR5720224,SRR5720225,SRR5720226,SRR5720227,SRR5720228,SRR5720229,SRR5720230,SRR5720231,SRR5720232,SRR5720233,SRR5720234,SRR5720235,SRR5720236,SRR5720237,SRR5720238,SRR5720239,SRR5720240,SRR5720241,SRR5720242,SRR5720243,SRR5720244,SRR5720245,SRR5720246,SRR5720247,SRR5720248,SRR5720249,SRR5720250,SRR5720251,SRR5720252,SRR5720253,SRR5720254,SRR5720255,SRR5720256,SRR5720257,SRR5720258,SRR5720259,SRR5720260,SRR5720261,SRR5720262,SRR5720263,SRR5720264,SRR5720265,SRR5720266,SRR5720267,SRR5720268,SRR5720269,SRR5720270,SRR5720271,SRR5720272,SRR5720273,SRR5720274,SRR5720275,SRR5720276,SRR5720277,SRR5720278,SRR5720279,SRR5720280,SRR5720281,SRR5720282,SRR5720283,SRR5720284,SRR5720285,SRR5720286,SRR5720287,SRR5720288,SRR5720289,SRR5720290,SRR5720291,SRR5720292,SRR5720293,SRR5720294,SRR5720295,SRR5720296,SRR5720297,SRR5720298,SRR5720299,SRR5720300,SRR5720301,SRR5720302,SRR5720303,SRR5720304,SRR5720305,SRR5720306,SRR5720307,SRR5720308,SRR5720309,SRR5720310,SRR5720311,SRR5720312,SRR5720313,SRR5720314,SRR5720315,SRR5720316,SRR5720317,SRR5720318,SRR5720319,SRR5720320,SRR5720321,SRR5720322,SRR5720323,SRR5720324,SRR5720325,SRR5720326,SRR5720327,SRR5720328,SRR5720329,SRR5720330,SRR5720331,SRR5720332,SRR5720333,SRR5720334,SRR5720335,SRR5720336,SRR5720337,SRR5720338,SRR5720339,SRR5720340,SRR5720341,SRR5720342,SRR5720343,SRR5720344
meta1 <- read_csv("SraRunTable.csv")

#from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP110813
meta2 <- read_csv("SraRunTable2.csv")

str(meta1$bottle_id)
str(meta1$collection_date)

#cleans up metadata bottle ID and collection date
meta1$bottle_id <- as.character(meta1$bottle_id)
meta1$collection_date <- dmy(meta1$collection_date)

str(meta1$bottle_id)
str(meta1$collection_date)

#combines meta1 and meta2
meta <- bind_rows(meta1, meta2)

#there is only one row for each sample in the metadata
head(meta)
nrow(meta)
meta %>% group_by(Run) %>% summarize(n = n()) %>% filter(n > 1)

#all of the biogeotraces samples have corresponding metadata in meta
#the rows in merged that are for OTUs that weren't found in biogeotraces samples 
#do not match to meta
head(merged)
head(meta)
merged %>% anti_join(meta, by = c("biogeotracesSample" = "Run")) %>% distinct(biogeotracesSample)

#these are the OTUs that were not found in any biogeotraces samples
merged %>% filter(is.na(biogeotracesSample)) %>% distinct(OTU_ID) %>% nrow()

#there are 10 biogeotraces samples in meta that are not in merged
nrow(meta)
merged %>% filter(!is.na(biogeotracesSample)) %>% distinct(biogeotracesSample) %>% nrow()

#there must not have been any matches to sean samples for these biogeotraces samples
meta %>% anti_join(merged, by = c("Run" = "biogeotracesSample")) %>% select(Run)

#adds metadata of each biogeotraces data to merged
nrow(merged)
merged <- merged %>% left_join(meta, by = c("biogeotracesSample" = "Run"))
nrow(merged)

#these are the OTUs in sean cultures that were not in biogeotraces samples
merged %>% filter(is.na(lat_lon)) %>% distinct(biogeotracesSample)

#makes depth variable numeric
str(merged$Depth)
merged$Depth <- parse_number(merged$Depth)
str(merged$Depth)

#gets rid of unnecessary variables
head(merged)
colnames(merged)
toPlot <- merged %>% select(sample, OTU_ID, biogeotracesSequence, biogeotracesSample, collection_date, 
                            cruise_id, Depth, geo_loc_name, lat_lon)


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

#there are not multiple matches between sean samples, sequences, OTUs, and biogeotraces samples which is good
head(toPlot)
toPlot %>% filter(is.na(biogeotracesSample))
toPlot %>% group_by(sample, OTU_ID, biogeotracesSequence, biogeotracesSample) %>% summarize(n = n()) %>% filter(n != 1)

#makes a variable for general cruise
head(toPlot)
toPlot %>% distinct(cruise_id)
toPlot <- toPlot %>% mutate(cruiseGeneral = str_extract(cruise_id, "BATS|HOT"))
toPlot <- toPlot %>% mutate(cruiseGeneral = ifelse(is.na(cruiseGeneral), cruise_id, cruiseGeneral))

toPlot %>% distinct(cruiseGeneral) %>% arrange(cruiseGeneral)
toPlot %>% filter(is.na(cruiseGeneral)) %>% distinct(cruise_id)

str(toPlot)

#arrange depths
str(toPlot$Depth)
depthOrder <- toPlot %>% arrange(Depth) %>% distinct(Depth)
depthOrder
depthOrder <- depthOrder %>% as.list()
depthOrder
depthOrder <- lapply(depthOrder, as.character)
depthOrder

str(toPlot$Depth)
toPlot$factorDepth <- as.character(toPlot$Depth)

str(toPlot$factorDepth)
toPlot %>% filter(is.na(factorDepth)) %>% distinct(OTU_ID) %>% nrow()

toPlot$factorDepth <- factor(toPlot$factorDepth, levels = depthOrder$Depth)
toPlot %>% filter(is.na(factorDepth)) %>% distinct(OTU_ID) %>% nrow()
toPlot %>% filter(is.na(Depth)) %>% distinct(OTU_ID) %>% nrow()

toPlot$factorDepth
depthOrder$Depth

toPlot %>% select(Depth, factorDepth)

#toPlot %>% distinct(factorDepth) %>% View()

#https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
toPlot %>% select(Depth, factorDepth) %>% head()
toPlot %>% mutate(factorDepth = as.numeric(as.character(factorDepth))) %>% select(Depth, factorDepth) %>% filter(Depth != factorDepth)

#biomG has all of the biogeotraces samples while toPlot has the sean samples
#and only the biogeotraces samples that have OTUs that are also in Sean samples

#one of these is NA
toPlot %>% distinct(biogeotracesSample) %>% nrow()
toPlot %>% distinct(biogeotracesSample) %>% filter(is.na(biogeotracesSample))

biomG %>% distinct(sample) %>% nrow()

#this has all of the OTUs that are present in Sean samples whether or not they 
#were present in biogeotraces samples
#after excluding Pro and Syn sequences
#toPlot has some OTUs that were not present in biogeotraces samples
toPlot %>% filter(!is.na(biogeotracesSequence)) %>% distinct(OTU_ID) %>% nrow()
toPlot %>% filter(is.na(biogeotracesSequence)) %>% distinct(OTU_ID) %>% nrow()
toPlot %>% distinct(OTU_ID) %>% nrow()

head(toPlot)
toPlot %>% distinct(sample) %>% nrow()
toPlot %>% distinct(OTU_ID) %>% nrow()

toPlot %>% filter(is.na(OTU_ID)) %>% nrow()

#loads taxonomies of OTUs
#from clusterSeanAndBiogeotracesSeqs.txt (checked)
tax <- read_tsv("sequencesForTree_noContamination_noMock_withProSyn_no1223_biogeotracesTaxonomy.tsv")

#gets rid of first row because it just has comments
tax[1,]
tax <- tax[-1,]

head(tax)

#gets rid of unnecessary variable
tax <- tax %>% select(1,2)

nrow(tax)

head(toPlot)
head(tax)

#add taxonomy of OTUs to toPlot
nrow(toPlot)
toPlot <- toPlot %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(toPlot)

toPlot %>% filter(is.na(Taxon)) %>% nrow()

head(toPlot)

#no Pro or Syn OTUs which makes sense because I already excluded Pro and Syn sequences 
#from Sean samples
toPlot$Taxon %>% head()
toPlot %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)
toPlot %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(Taxon)

#gets rid of unnecessary variable
head(toPlot)
toPlot <- toPlot %>% select(-Taxon)

#toPlot includes just biogeotraces samples that have OTUs also found in sean samples

#one of these is NA
toPlot %>% distinct(biogeotracesSample) %>% nrow()
toPlot %>% distinct(biogeotracesSample) %>% filter(is.na(biogeotracesSample))

nrow(meta)

#correct number of sean samples
toPlot %>% distinct(sample) %>% nrow()

toPlot <- toPlot %>% ungroup()


###high occupancy, abundance OTUs

#loads the OTUs that have mean relative abundance and are present in many cultures 
#in occupancy abundance plot
#these OTU IDs are from the sean and biogeotraces OTUs combined sequences
#I got these OTU IDs in ubiquityAnalysisSeanBiogeotraces_clustered.R
#from ubiquityAnalysisSeanBiogeotraces_clustered.R
occAbun <- read_csv("occupancyAbundanceNonProNonSynOTUsInSeanCultures_seanAndBiogeotracesClusteredTogether.csv")

head(occAbun)
nrow(occAbun)

#tax is from the sean and biogeotraces OTUs combined sequences 
head(tax)
nrow(tax)

#adds taxonomy to occAbun
nrow(occAbun)
occAbun <- occAbun %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(occAbun)

occAbun %>% filter(is.na(Taxon))
head(occAbun)

occAbun %>% distinct(Taxon)

occAbun

#makes a taxonomy variable without "k__, p__..."
occAbun <- occAbun %>% mutate(newTaxon = str_replace_all(Taxon, "([a-z]__)", ""))

#gets rid of ";" in fixed taxonomy variable
occAbun <- occAbun %>% mutate(newTaxon = str_replace_all(newTaxon, ";", ""))
#occAbun %>% distinct(Taxon, newTaxon) %>% View()

#makes a shortened taxonomy variable
occAbun <- occAbun %>% mutate(shortTaxon = str_extract(newTaxon, "[a-zA-Z0-9\\-\\_]*\\s*$"))
#occAbun %>% distinct(Taxon, shortTaxon) %>% View()

occAbun %>% distinct(shortTaxon)

#gets rid of extra spaces at the of shortened taxonomies
str_extract(occAbun$shortTaxon, "\\s$")
occAbun <- occAbun %>% mutate(shortTaxon = str_replace_all(shortTaxon, "\\s$", ""))

occAbun %>% distinct(shortTaxon)
#occAbun %>% distinct(Taxon, shortTaxon) %>% View()

#makes a variable that has both the OTU and shortened taxonomy
head(occAbun)
occAbun <- occAbun %>% mutate(otuSequence = str_c(OTU_ID, shortTaxon, sep = ": "))

occAbun

##gets biogeotraces samples that have occupancy, abundance OTUs
toPlot %>% distinct(OTU_ID) %>% nrow()
toPlot %>% head()

toPlot <- toPlot %>% semi_join(occAbun, by = c("OTU_ID"))
toPlot %>% distinct(OTU_ID) %>% nrow()

nrow(occAbun)
toPlot %>% distinct(OTU_ID) %>% nrow()

#all of the occupancy, abundance sean OTUs were present in at least one biogeotraces
#sample
toPlot %>% head()
toPlot %>% filter(is.na(biogeotracesSample)) %>% nrow()

#gets distinct samples, OTU, biogeotraces sequences, and biogeotraces samples matches
colnames(toPlot)
toPlot <- toPlot %>% distinct(sample, OTU_ID, biogeotracesSequence, biogeotracesSample)


#gets the number of matches (number of biogeotraces sequences in the biogeotraces sample) 
#for each biogeotraces sample and occupancy, abundance OTU
head(toPlot)
toPlot
toPlot %>% group_by(OTU_ID, biogeotracesSample) %>% distinct(biogeotracesSequence) %>% summarize(n = n())

#how many biogeotraces samples each occupancy, abundance OTU is found in
head(toPlot)
toPlot %>% group_by(OTU_ID) %>% distinct(biogeotracesSample) %>% summarize(n = n())

#gets the number of biogeotraces sequences for each biogeotraces sample and OTU
head(toPlot)
toPlot <- toPlot %>% ungroup()
ubiSummary <- toPlot %>% group_by(OTU_ID, biogeotracesSample) %>% distinct(biogeotracesSequence) %>% 
  summarize(n = n()) %>% arrange(desc(n))

ubiSummary <- ubiSummary %>% ungroup()

head(ubiSummary)

ubiSummary %>% filter(!str_detect(biogeotracesSample, "SRR"))

#all of the occupancy, abundance OTUs were found in at least one biogeotraces sample
ubiSummary %>% filter(is.na(biogeotracesSample))

head(occAbun)
head(ubiSummary)
occAbun %>% distinct(OTU_ID)
ubiSummary %>% distinct(OTU_ID)

#adds otu/taxonomy variable to ubiSummary
nrow(ubiSummary)
ubiSummary <- ubiSummary %>% left_join(occAbun %>% select(OTU_ID, otuSequence), by = c("OTU_ID"))
nrow(ubiSummary)

head(ubiSummary)

ubiSummary %>% filter(is.na(otuSequence))

ubiSummary %>% distinct(OTU_ID, otuSequence)

#not all of the biogeotraces samples are in ubiSummary
ubiSummary %>% distinct(biogeotracesSample) %>% nrow()
meta %>% nrow()

#adds all of the biogeotraces samples to the dataframe
ubiSummary <- ubiSummary %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))

head(ubiSummary)

#now all of the biogeotraces samples are in ubiSummary
ubiSummary %>% distinct(biogeotracesSample) %>% nrow()

#the number of biogeotraces samples that each OTU is found in
ubiSummary %>% group_by(OTU_ID) %>% summarize(n = n())
ubiSummary %>% group_by(OTU_ID) %>% distinct(biogeotracesSample) %>% summarize(n = n())

nrow(meta)
head(ubiSummary)
ubiSummary %>% filter(is.na(biogeotracesSample)) %>% nrow()
ubiSummary %>% distinct(biogeotracesSample) %>% nrow()
ubiSummary %>% group_by(OTU_ID) %>% distinct(biogeotracesSample) %>% summarize(n = n())

ubiSummary %>% filter(is.na(lat_lon)) %>% nrow()


###gets number of non-Pro, non-Syn sequences for each biogeotraces sample

#biomG includes Pro and Syn sequences
biomG %>% head()
biomG %>% filter(!str_detect(sample, "SRR"))
biomG %>% filter(str_detect(sample, "contig")) %>% distinct(sample) %>% nrow()

numSeq <- biomG

#gets rid of unnecessary variable
numSeq %>% distinct(present)
numSeq <- numSeq %>% select(-present)
numSeq %>% head()

#makes a variable for the biogeotraces sample
numSeq <- numSeq %>% mutate(biogeotracesSample = str_extract(sample, "[a-zA-Z0-9]*"))
numSeq %>% head()

#gets just biogeotraces samples
numSeq %>% filter(!str_detect(biogeotracesSample, "SRR")) %>% distinct(biogeotracesSample)
numSeq <- numSeq %>% filter(str_detect(biogeotracesSample, "SRR"))

nrow(numSeq)
numSeq %>% distinct(biogeotracesSample) %>% nrow()
numSeq %>% distinct(OTU_ID) %>% nrow()

head(numSeq)
head(tax)

#adds taxonomy to numSeq
nrow(numSeq)
numSeq <- numSeq %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(numSeq)

head(numSeq)
numSeq %>% filter(is.na(Taxon))

numSeq %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)
numSeq %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(Taxon)

#excludes Pro and Syn OTUs from numSeq
numSeq <- numSeq %>% filter(!str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc"))

numSeq %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)

#all of the biogeotraces samples have at least one non-Pro, non-Syn OTU
numSeq %>% distinct(biogeotracesSample) %>% nrow()

numSeq <- numSeq %>% ungroup() 

#all of the biogeotraces sequences are only present once in biomG which 
#makes sense
head(numSeq)
numSeq %>% group_by(sample) %>% summarize(n = n()) %>% filter(n > 1)

#gets the number of non-Pro, non-Syn sequences in each biogeotraces sample
numSeq <- numSeq %>% group_by(biogeotracesSample) %>% distinct(sample) %>% summarize(numSeq = n())

nrow(numSeq)
head(numSeq)
numSeq %>% filter(is.na(numSeq))
str(numSeq)

#write_csv(numSeq, "numSequencesInEachBiogeotraceSampleNoProSyn.csv")

head(ubiSummary)
head(numSeq)

#adds number of sequences in each biogeotraces sample to ubiSummary
nrow(ubiSummary)
ubiSummary <- ubiSummary %>% left_join(numSeq, by = c("biogeotracesSample"))
nrow(ubiSummary)

head(ubiSummary)
ubiSummary %>% filter(is.na(numSeq)) %>% nrow()

#n is the abundance of the OTU in the biogeotraces sample
#the NA rows are the rows for biogeotraces samples that didn't have
#the occupancy, abundance OTUs
head(ubiSummary)
ubiSummary %>% filter(is.na(n)) %>% nrow()

NA/2

#makes a variable for the number of sequences of an OTU present in the 
#biogeotraces sample divided by the total number of non-Pro, non-Syn sequences 
#in the sample
head(ubiSummary)
str(ubiSummary$n)
str(ubiSummary$numSeq)
ubiSummary <- ubiSummary %>% mutate(propOfBiogeotraces = n/numSeq)

ubiSummary %>% filter(is.na(propOfBiogeotraces)) %>% nrow()

ubiSummary %>% arrange(desc(propOfBiogeotraces)) %>% select(propOfBiogeotraces)

#gets rid of unnecessary variables
colnames(ubiSummary)
ubiSummary <- ubiSummary %>% select(OTU_ID, otuSequence, propOfBiogeotraces, biogeotracesSample, Depth, lat_lon)
colnames(ubiSummary)

head(ubiSummary)

meta %>% distinct(BioSample) %>% nrow()
ubiSummary %>% distinct(biogeotracesSample) %>% nrow()

#these are the biogeotraces samples that did not have the occupancy, abundnace culture OTUs
ubiSummary %>% filter(is.na(propOfBiogeotraces)) %>% nrow()

#replaces NA propOfBiogeotraces values with 0 
ubiSummary <- ubiSummary %>% mutate(propOfBiogeotraces = ifelse(is.na(propOfBiogeotraces), 0, propOfBiogeotraces))

ubiSummary %>% filter(propOfBiogeotraces == 0) %>% distinct(otuSequence)

ubiSummary %>% filter(propOfBiogeotraces == 0) %>% nrow()

ubiSummary %>% distinct(OTU_ID, otuSequence)

#these are how many biogeotraces samples each OTU was in
ubiSummary %>% group_by(OTU_ID, otuSequence) %>% summarize(n = n())
ubiSummary %>% group_by(OTU_ID, otuSequence) %>% distinct(biogeotracesSample) %>% summarize(n = n())

ubiSummary %>% distinct(otuSequence)
ubiSummary %>% filter(is.na(otuSequence)) %>% nrow()

#replaces NA in rows for biogeotraces samples that didn't have one of the OTUs 
#with "All BioGeoTraces sites"
ubiSummary <- ubiSummary %>% mutate(otuSequence = ifelse(is.na(otuSequence), "All BioGeoTraces sites", otuSequence))
ubiSummary %>% distinct(otuSequence)
ubiSummary %>% filter(otuSequence == "All BioGeoTraces sites") %>% nrow()

ubiSummary %>% filter(is.na(propOfBiogeotraces))

#replaces ":" with ":\n" in otuSequence IDs
ubiSummary %>% distinct(otuSequence)
ubiSummary <- ubiSummary %>% mutate(otuSequence = str_replace(otuSequence, ": ", ":\n"))
ubiSummary %>% distinct(otuSequence)

ubiSummary %>% filter(propOfBiogeotraces == 0) %>% distinct(otuSequence)

colnames(meta)

##re-adds geolocation name and collection date to ubiSummary
meta <- meta %>% select(Run, geo_loc_name, collection_date)
head(meta)
nrow(meta)

ubiSummary %>% head()

nrow(ubiSummary)
ubiSummary <- ubiSummary %>% left_join(meta, by = c("biogeotracesSample" = "Run"))
nrow(ubiSummary)

ubiSummary %>% filter(is.na(geo_loc_name))

ubiSummary %>% distinct(OTU_ID, otuSequence)

#130/610 of the biogeotraces samples are in HOT or BATS
meta %>% distinct(Run) %>% nrow()
meta %>% filter(geo_loc_name %in% 
  c("Atlantic Ocean: Sargasso Sea\\, BATS", "Pacific Ocean: North Pacific Subtropical Gyre\\, Station ALOHA")) %>% 
  distinct(Run) %>% nrow()

ubiSummaryFull <- ubiSummary

ubiSummary %>% filter(geo_loc_name %in% 
  c("Atlantic Ocean: Sargasso Sea\\, BATS", "Pacific Ocean: North Pacific Subtropical Gyre\\, Station ALOHA")) %>% 
  distinct(biogeotracesSample) %>% nrow()

#gets just the rows in ubiSummary for the HOT or BATS biogeotraces samples
ubiSummary <- ubiSummary %>% filter(geo_loc_name %in% 
     c("Atlantic Ocean: Sargasso Sea\\, BATS", "Pacific Ocean: North Pacific Subtropical Gyre\\, Station ALOHA"))

ubiSummary %>% distinct(geo_loc_name)

#now some of the OTUs are missing
ubiSummary %>% distinct(OTU_ID)
ubiSummaryFull %>% distinct(OTU_ID)
missing <- ubiSummaryFull %>% anti_join(ubiSummary, by = c("OTU_ID")) %>% distinct(OTU_ID, otuSequence)

ubiSummary %>% distinct(OTU_ID)
missing

head(ubiSummary)
missing

#adds missing OTUs to ubiSummary
nrow(ubiSummary)
ubiSummary <- bind_rows(ubiSummary, missing)
nrow(ubiSummary)

ubiSummary %>% distinct(OTU_ID)

#one of these is NA
ubiSummary %>% distinct(biogeotracesSample) %>% nrow()
ubiSummary %>% filter(is.na(biogeotracesSample))

head(ubiSummary)
colnames(ubiSummary)

#gets rid of unnecessary variables so that I can spread data
ubiSummary <- ubiSummary %>% select(OTU_ID:biogeotracesSample)
head(ubiSummary)

ubiSummary %>% group_by(OTU_ID, biogeotracesSample) %>% summarize(n = n()) %>% filter(n != 1)

#spreads data so that there is a value of 0 for propOfBiogeotraces if the OTU is not 
#present in a biogeotraces sample
ubiSummary_spread <- ubiSummary %>% spread(key = biogeotracesSample, value = propOfBiogeotraces, fill = 0)

colnames(ubiSummary_spread)

#regathers data
ubiSummary <- ubiSummary_spread %>% gather("SRR5720219":"<NA>", key = "biogeotracesSample", value = "propOfBiogeotraces")

head(ubiSummary)

ubiSummary %>% filter(biogeotracesSample == "<NA>")
ubiSummary <- ubiSummary %>% filter(biogeotracesSample != "<NA>")

ubiSummary %>% filter(is.na(propOfBiogeotraces))

ubiSummary %>% distinct(OTU_ID)

#there is a propOfBiogeotraces value for each of the 130 HOT and BATS biogeotraces 
#samples for each OTU
ubiSummary %>% group_by(OTU_ID) %>% summarize(n = n())

#gets rid of NA OTU because it is no longer needed 
#I just needed it to spread data
ubiSummary <- ubiSummary %>% filter(!is.na(OTU_ID))

ubiSummary %>% group_by(OTU_ID) %>% summarize(n = n())

#re-adds biogeotraces geolocation and collection date info 
#to ubiSummary

meta %>% filter(geo_loc_name %in% 
  c("Atlantic Ocean: Sargasso Sea\\, BATS", "Pacific Ocean: North Pacific Subtropical Gyre\\, Station ALOHA")) %>% 
  distinct(Run) %>% nrow()

nrow(ubiSummary)
ubiSummary <- ubiSummary %>% left_join(meta, by = c("biogeotracesSample" = "Run"))
nrow(ubiSummary)

ubiSummary %>% distinct(geo_loc_name, collection_date)

#makes a variable for the biogeotraces collection year
ubiSummary <- ubiSummary %>% mutate(year = str_extract(collection_date, "[0-9]*"))

ubiSummary %>% distinct(geo_loc_name, collection_date, year)
ubiSummary %>% distinct(geo_loc_name, year)
ubiSummary %>% distinct(geo_loc_name)

ubiSummary %>% head()
ubiSummary %>% distinct(OTU_ID) %>% nrow()
ubiSummary %>% group_by(OTU_ID) %>% summarize(n = n())
ubiSummary %>% distinct(biogeotracesSample) %>% nrow()

ubiSummary %>% distinct(OTU_ID)
ubiSummary %>% distinct(geo_loc_name)

ubiSummary %>% group_by(geo_loc_name, year) %>% distinct(biogeotracesSample) %>% summarize(n = n())

ubiSummary

ubiSummary %>% 
  filter(geo_loc_name == "Atlantic Ocean: Sargasso Sea\\, BATS") %>%
  ggplot() + geom_boxplot(aes(x = year, y = propOfBiogeotraces)) +
  labs(y = expression(paste("Proportion of non-Pro., non-Syn. sequences in BioGeoTraces sample"))) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence, scales = 'free') +
  labs(title = "BATs: OTUs with high mean relative abundance across cultures that\nwere also present in a large proportion of cultures")

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_byYearBATS.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_byYearBATS.svg", dpi = 600, height = 10, width = 20)


ubiSummary %>% 
  filter(geo_loc_name == "Pacific Ocean: North Pacific Subtropical Gyre\\, Station ALOHA") %>%
  ggplot() + geom_boxplot(aes(x = year, y = propOfBiogeotraces)) +
  labs(y = expression(paste("Proportion of non-Pro., non-Syn. sequences in BioGeoTraces sample"))) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence, scales = 'free') +
  labs(title = "HOT: OTUs with high mean relative abundance across cultures that\nwere also present in a large proportion of cultures")

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_byYearHOT.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_byYearHOT.svg", dpi = 600, height = 10, width = 20)

