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
biomG %>% filter(present == 1) %>% group_by(sample) %>% summarize(n = n()) %>% filter(n != 1)

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

seqG %>% distinct(sequence) %>% nrow()
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

#all of the samples have at least one sequence that corresponds to one of the biogeotrace sequences
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

#there are NAs when a sequence in a sean sample did not correspond to any of the biogeotrace sequences
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

#these are the sequences and OTUs that are present in samples and some 
#are also present in biogeotraces samples
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

#I care about how many sequences in each biogeotraces sample match to each sean sample 
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
toPlot <- toPlot %>% mutate(cruiseGeneral = str_extract(cruise_id, "BATS|HOT|PE|KN"))
toPlot <- toPlot %>% mutate(cruiseGeneral = ifelse(is.na(cruiseGeneral), cruise_id, cruiseGeneral))

toPlot %>% distinct(cruiseGeneral) %>% arrange(cruiseGeneral)

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
toPlot %>% filter(is.na(factorDepth)) %>% nrow()

toPlot$factorDepth <- factor(toPlot$factorDepth, levels = depthOrder$Depth)
toPlot %>% filter(is.na(factorDepth)) %>% nrow()
toPlot %>% filter(is.na(Depth)) %>% nrow()

toPlot$factorDepth
depthOrder$Depth

toPlot %>% select(Depth, factorDepth)

#toPlot %>% distinct(factorDepth) %>% View()

#https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
toPlot %>% select(Depth, factorDepth) %>% head()
toPlot %>% mutate(factorDepth = as.numeric(as.character(factorDepth))) %>% select(Depth, factorDepth) %>% filter(Depth != factorDepth)

#biomG has all of the biogeotraces samples while toPlot has the sean samples
#and only the biogeotraces samples that have OTUs that are also in Sean samples
toPlot %>% distinct(biogeotracesSample) %>% nrow()
biomG %>% distinct(sample) %>% nrow()

#this has all of the OTUs that are present in Sean samples whether or not they 
#were present in biogeotraces samples
#after excluding Pro and Syn sequences
#toPlot has some OTUs that were not present in biogeotraces samples
toPlot %>% filter(!is.na(biogeotracesSequence)) %>% distinct(OTU_ID) %>% nrow()
toPlot %>% filter(is.na(biogeotracesSequence)) %>% distinct(OTU_ID) %>% nrow()

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

#add taxonomy of OTUs to toPlot
nrow(toPlot)
toPlot <- toPlot %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(toPlot)

toPlot %>% filter(is.na(Taxon)) %>% nrow()

head(toPlot)

#no Pro or Syn OTUs which makes sense because I already excluded Pro and Syn sequences 
#from Sean samples
toPlot %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)
toPlot %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(Taxon)

#gets rid of unnecessary variable
head(toPlot)
toPlot <- toPlot %>% select(-Taxon)

#toPlot includes just biogeotraces samples that have OTUs also found in sean samples
toPlot %>% distinct(biogeotracesSample) %>% nrow()
nrow(meta)

#correct number of sean samples
toPlot %>% distinct(sample) %>% nrow()


###ubiquity: sequences that are common in cultures

head(toPlot)
toPlot %>% filter(is.na(biogeotracesSequence))
toPlot %>% distinct(OTU_ID) %>% nrow()
toPlot %>% distinct(sample) %>% nrow()

toPlot %>% filter(is.na(OTU_ID)) %>% nrow()

#gets the OTUs that are found in the most cultures
head(toPlot)
toPlot <- toPlot %>% ungroup()
ubi <- toPlot %>% group_by(OTU_ID) %>% distinct(sample) %>% summarize(numSamples = n())
head(ubi)
nrow(ubi)
write_csv(ubi %>% arrange(desc(numSamples)) %>% slice(1:10), 
    "mostUbiquitiousNonProNonSynOTUsInSeanCultures_seanAndBiogeotracesClusteredTogether.csv")
ubi <- ubi %>% arrange(desc(numSamples)) %>% slice(1:5)
ubi



head(tax)
nrow(tax)

#adds taxonomy to the most ubiquitious OTUS
nrow(ubi)
ubi <- ubi %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(ubi)

ubi %>% filter(is.na(Taxon))
head(ubi)

ubi %>% distinct(Taxon)

#makes a taxonomy variable without "k__, p__..."
ubi <- ubi %>% mutate(newTaxon = str_replace_all(Taxon, "([a-z]__)", ""))

#gets rid of ";" in fixed taxonomy variable
ubi <- ubi %>% mutate(newTaxon = str_replace_all(newTaxon, ";", ""))
#ubi %>% distinct(Taxon, newTaxon) %>% View()

#makes a shortened taxonomy variable
ubi <- ubi %>% mutate(shortTaxon = str_extract(newTaxon, "[a-zA-Z0-9\\-\\_]*\\s*$"))
#ubi %>% distinct(Taxon, shortTaxon) %>% View()

ubi %>% distinct(shortTaxon)

#gets rid of extra spaces at the of shortened taxonomies
str_extract(ubi$shortTaxon, "\\s$")
ubi <- ubi %>% mutate(shortTaxon = str_replace_all(shortTaxon, "\\s$", ""))

ubi %>% distinct(shortTaxon)
#ubi %>% distinct(Taxon, shortTaxon) %>% View()

#makes a variable that has both the OTU and shortened taxonomy
ubi <- ubi %>% mutate(otuSequence = str_c(OTU_ID, shortTaxon, sep = ": "))

ubi

##gets biogeotraces samples that each ubiquitious OTU across sean samples 
##matches up to
toPlot %>% distinct(OTU_ID)
toPlot %>% head()
toPlot_forClass <- toPlot

toPlot_occAbun <- toPlot

toPlot <- toPlot %>% semi_join(ubi, by = c("OTU_ID"))
toPlot %>% distinct(OTU_ID)

nrow(ubi)
toPlot %>% distinct(OTU_ID) %>% nrow()

#all of the most ubiquitious sean OTUs were present in at least one biogeotraces
#sample
toPlot %>% head()
toPlot %>% filter(is.na(biogeotracesSample)) %>% nrow()

#gets distinct samples, OTU, biogeotraces sequences, and biogeotraces samples matches
colnames(toPlot)
toPlot <- toPlot %>% distinct(sample, OTU_ID, biogeotracesSequence, biogeotracesSample)


#gets the number of matches (number of biogeotraces sequences in the biogeotraces sample) 
#for each biogeotraces sample and ubiquitious OTU
head(toPlot)
toPlot %>% group_by(OTU_ID, biogeotracesSample) %>% distinct(biogeotracesSequence) %>% summarize(n = n())

#how many biogeotraces samples each ubiquitious OTU is found in
head(toPlot)
toPlot %>% group_by(OTU_ID) %>% distinct(biogeotracesSample) %>% summarize(n = n())

#gets the number of biogeotraces sequences for each biogeotraces sample and OTU
head(toPlot)
toPlot <- toPlot %>% ungroup()
ubiSummary <- toPlot %>% group_by(OTU_ID, biogeotracesSample) %>% distinct(biogeotracesSequence) %>% summarize(n = n()) %>% arrange(desc(n))

ubiSummary <- ubiSummary %>% ungroup()

head(ubiSummary)

ubiSummary %>% filter(!str_detect(biogeotracesSample, "SRR"))

#all of the ubiquitious OTUs were found in at least one biogeotraces sample
ubiSummary %>% filter(is.na(biogeotracesSample))

head(ubi)
head(ubiSummary)
ubi %>% distinct(OTU_ID)
ubiSummary %>% distinct(OTU_ID)

#adds otu/taxonomy variable to ubiSummary
nrow(ubiSummary)
ubiSummary <- ubiSummary %>% left_join(ubi %>% select(OTU_ID, otuSequence), by = c("OTU_ID"))
nrow(ubiSummary)

ubiSummary %>% filter(is.na(otuSequence))

ubiSummary %>% distinct(OTU_ID, otuSequence)

#makes a separate dataframe for each OTU
ubiSummary_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- ubiSummary %>% filter(OTU_ID == "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d")
ubiSummary_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- ubiSummary %>% filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf")
ubiSummary_6e447d49df572a0622adad2e68409bcf687108e8 <- ubiSummary %>% filter(OTU_ID == "6e447d49df572a0622adad2e68409bcf687108e8")
ubiSummary_02ba11823bfe93a0250096d95e07e745e40736bc <- ubiSummary %>% filter(OTU_ID == "02ba11823bfe93a0250096d95e07e745e40736bc")
ubiSummary_1618797b285d4650872ef897004644b2990e6d73 <- ubiSummary %>% filter(OTU_ID == "1618797b285d4650872ef897004644b2990e6d73")


#adds all of the biogeotraces samples to each OTU dataframe
ubiSummary_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- ubiSummary_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- ubiSummary_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_6e447d49df572a0622adad2e68409bcf687108e8 <- ubiSummary_6e447d49df572a0622adad2e68409bcf687108e8 %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_02ba11823bfe93a0250096d95e07e745e40736bc <- ubiSummary_02ba11823bfe93a0250096d95e07e745e40736bc %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_1618797b285d4650872ef897004644b2990e6d73 <- ubiSummary_1618797b285d4650872ef897004644b2990e6d73 %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))

#combines all of the dataframes
ubiSummary <- bind_rows(ubiSummary_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d, 
                        ubiSummary_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf, 
                        ubiSummary_6e447d49df572a0622adad2e68409bcf687108e8, 
                        ubiSummary_02ba11823bfe93a0250096d95e07e745e40736bc, 
                        ubiSummary_1618797b285d4650872ef897004644b2990e6d73)

head(ubiSummary)

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

write_csv(numSeq, "numSequencesInEachBiogeotraceSampleNoProSyn.csv")

head(ubiSummary)
head(numSeq)

#adds number of sequences in each biogeotraces sample to ubiSummary
nrow(ubiSummary)
ubiSummary <- ubiSummary %>% left_join(numSeq, by = c("biogeotracesSample"))
nrow(ubiSummary)

head(ubiSummary)
ubiSummary %>% filter(is.na(numSeq)) %>% nrow()

#n is the abundance of the OTU across all biogeotraces samples 
#these are for the rows for biogeotraces samples that didn't have
#the ubiquitious OTUs
ubiSummary %>% filter(is.na(n)) %>% nrow()

NA/2

#makes a variable for the number of sequences of an OTU present in the 
#biogeotraces sample divided by the total number of non-Pro, non-Syn sequences 
#in the sample
head(ubiSummary)
ubiSummary <- ubiSummary %>% mutate(propOfBiogeotraces = n/numSeq)

ubiSummary %>% filter(is.na(propOfBiogeotraces)) %>% nrow()

ubiSummary %>% arrange(desc(propOfBiogeotraces)) %>% select(propOfBiogeotraces)

#gets rid of unnecessary variables
colnames(ubiSummary)
ubiSummary <- ubiSummary %>% select(OTU_ID, otuSequence, propOfBiogeotraces, biogeotracesSample, Depth, lat_lon)
colnames(ubiSummary)

head(ubiSummary)

##makes lat, long variables

ubiSummary %>% select(lat_lon) %>% head()
ubiSummary <- ubiSummary %>% mutate(lat = str_extract(lat_lon, "[0-9]*\\.*[0-9]*\\s(N|S)"))
ubiSummary %>% filter(is.na(lat))
ubiSummary %>% select(lat_lon, lat) %>% head()
ubiSummary %>% select(lat_lon, lat) %>% filter(str_detect(lat_lon, "S")) %>% head()

ubiSummary <- ubiSummary %>% mutate(long = str_extract(lat_lon, "[0-9]*\\.*[0-9]*\\s(E|W)"))
ubiSummary %>% filter(is.na(long))
ubiSummary %>% select(lat_lon, lat, long) %>% head()
ubiSummary %>% select(lat_lon, lat, long) %>% filter(str_detect(lat_lon, "E")) %>% head()

#converts lat and log variables numeric
ubiSummary <- ubiSummary %>% mutate(lat_2 = parse_number(lat))
ubiSummary <- ubiSummary %>% mutate(long_2 = parse_number(long))

ubiSummary %>% select(lat_lon, lat, long, lat_2, long_2) %>% head()

#makes lat, long values negative if their direction is S or W
ubiSummary <- ubiSummary %>% mutate(lat_2 = ifelse(str_detect(lat, "S"), -1*lat_2, lat_2))
ubiSummary <- ubiSummary %>% mutate(long_2 = ifelse(str_detect(long, "W"), -1*long_2, long_2))

ubiSummary %>% select(lat_lon, lat, long, lat_2, long_2) %>% head()
ubiSummary %>% select(lat_lon, lat, long, lat_2, long_2) %>% filter(str_detect(lat_lon, "S")) %>% head()
ubiSummary %>% select(lat_lon, lat, long, lat_2, long_2) %>% filter(str_detect(lat_lon, "E")) %>% head()

#gets rid of unnecessary variables
ubiSummary <- ubiSummary %>% select(-c(lat_lon, lat, long))

#renames lat, long variables
colnames(ubiSummary)
colnames(ubiSummary)[6:7]
colnames(ubiSummary)[6:7] <- c("lat", "long")

meta %>% distinct(BioSample) %>% nrow()
ubiSummary %>% distinct(biogeotracesSample) %>% nrow()

#these are the biogeotraces samples that did not have the ubiquitious culture OTUs
ubiSummary %>% filter(is.na(propOfBiogeotraces))

#replaces NA propOfBiogeotraces values with 0 
ubiSummary <- ubiSummary %>% mutate(propOfBiogeotraces = ifelse(is.na(propOfBiogeotraces), 0, propOfBiogeotraces))

ubiSummary %>% filter(propOfBiogeotraces == 0) %>% distinct(otuSequence)

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


##puts the OTUs in order of how many samples they were found in
ubi
order <- ubi %>% arrange(desc(numSamples))
order
order <- order %>% select(otuSequence)
order
order <- order$otuSequence
order <- order %>% as.list()
order
length(order)
order[5]
order[6]
order[6] <- "All BioGeoTraces sites"
length(order)
order
order <- str_replace(order, ": ", ":\n")
order
str(order)
order <- order %>% as.list()
order

head(ubiSummary)
ubiSummary %>% distinct(otuSequence)
ubiSummary$otuSequence <- factor(ubiSummary$otuSequence, levels = order)
levels(ubiSummary$otuSequence)
ubi

ubiSummary %>% filter(propOfBiogeotraces == 0) %>% distinct(otuSequence)

library(ggworldmap)

# project the data
proj="robin"
#long_0 = -90
ubiSummary_2 <- project(ubiSummary, proj)
head(ubiSummary_2)

###is it okay to get rid of the long_0 in the projection??
# plot the world
ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence) +
  labs(title = "OTUs present in the most cultures,\narranged by number of cultures they were found in")
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousOTUs.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousOTUs.svg", dpi = 600, height = 10, width = 20)



ubiSummary$otuSequence %>% levels()

ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(otuSequence == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf:\nRhodobacteraceae") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(otuSequence == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf:\nRhodobacteraceae") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousSeqs_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf_Rhodobacteraceae.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousSeqs_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf_Rhodobacteraceae.svg", dpi = 600, height = 5, width = 10)


ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(otuSequence == "1618797b285d4650872ef897004644b2990e6d73:\nThalassospira") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(otuSequence == "1618797b285d4650872ef897004644b2990e6d73:\nThalassospira") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousSeqs_1618797b285d4650872ef897004644b2990e6d73_Thalassospira.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousSeqs_1618797b285d4650872ef897004644b2990e6d73_Thalassospira.svg", dpi = 600, height = 5, width = 10)


ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(otuSequence == "6e447d49df572a0622adad2e68409bcf687108e8:\nGammaproteobacteria") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(otuSequence == "6e447d49df572a0622adad2e68409bcf687108e8:\nGammaproteobacteria") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousSeqs_6e447d49df572a0622adad2e68409bcf687108e8_Gammaproteobacteria.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousSeqs_6e447d49df572a0622adad2e68409bcf687108e8_Gammaproteobacteria.svg", dpi = 600, height = 5, width = 10)


ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(otuSequence == "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d:\nAlteromonadaceae") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(otuSequence == "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d:\nAlteromonadaceae") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousSeqs_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d_Alteromonadaceae.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousSeqs_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d_Alteromonadaceae.svg", dpi = 600, height = 5, width = 10)


ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(otuSequence == "02ba11823bfe93a0250096d95e07e745e40736bc:\nMethylophaga") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(otuSequence == "02ba11823bfe93a0250096d95e07e745e40736bc:\nMethylophaga") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousSeqs_02ba11823bfe93a0250096d95e07e745e40736bc_Methylophaga.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousSeqs_02ba11823bfe93a0250096d95e07e745e40736bc_Methylophaga.svg", dpi = 600, height = 5, width = 10)



###occupancy, abundance plots 

#the OTU IDs in occAbunOTU are different than in toPlot_occAbun because 
#occAbunOTU has OTUs based on clustering of just sean culture sequences 
#while toPlot_occAbun has OTUs based on clustering of sean culture and 
#biogeotraces sequences
#from occupancyAbundancePlots.R
occAbunOTU <- read_csv("OTUs_fromOccupancyAbundancePlots.csv")
occAbunOTU %>% head()
occAbunOTU %>% nrow()
toPlot_occAbun %>% distinct(OTU_ID) %>% nrow()

#biomGSeanOTU has OTUs just based on clustering of sean sequences
#from occupancyAbundancePlots.R
biomGSeanOTU <- read_csv("biomG_fromOccupancyAbundancePlots.csv")

biomGSeanOTU %>% head()
occAbunOTU %>% head()

biomGSeanOTU %>% semi_join(occAbunOTU, by = "OTU_ID") %>% distinct(OTU_ID)

#gets just the rows in biomGSeanOTU that are for OTUs in occAbunOTU
biomGSeanOTU <- biomGSeanOTU %>% semi_join(occAbunOTU, by = "OTU_ID") 
biomGSeanOTU %>% distinct(OTU_ID) %>% nrow()

#seanBiom has OTUs based on clustering of sean and biogeotraces sequences
seanBiom %>% distinct(present)
head(seanBiom)
head(biomGSeanOTU)

##makes sean sequence IDs in biomGSeanOTU match the format of seanBiom
biomGSeanOTU <- biomGSeanOTU %>% mutate(sample = str_replace(sample, "seq", ""))
head(biomGSeanOTU)
biomGSeanOTU <- biomGSeanOTU %>% mutate(sample = str_c("contig_", sample))
head(biomGSeanOTU$sample)

head(seanBiom)
head(biomGSeanOTU)

#the sequences are grouped into 10 clusters just based on sean data 
#and 9 clusters based on sean and biogeotraces data
biomGSeanOTU %>% distinct(sample) %>% nrow()
biomGSeanOTU %>% distinct(OTU_ID) %>% nrow()
seanBiom %>% semi_join(biomGSeanOTU, by = c("sample")) %>% distinct(sample) %>% nrow()
seanBiom %>% semi_join(biomGSeanOTU, by = c("sample")) %>% distinct(OTU_ID) %>% nrow()

#gets the sequences of interest because the OTU IDs don't correspond
seanBiom <- seanBiom %>% semi_join(biomGSeanOTU, by = c("sample"))

head(seanBiom)
seanBiom %>% distinct(present)
seanBiom %>% distinct(OTU_ID) %>% nrow()
seanBiom %>% distinct(sample) %>% nrow()



#gets the rows in toPlot that correspond to the OTUs in seanBiom
#9/9 of the OTUs of interest are also present in biogeotraces samples
head(toPlot_occAbun)
head(seanBiom)
toPlot_occAbun %>% distinct(OTU_ID) %>% nrow()
toPlot %>% distinct(OTU_ID) %>% nrow()

toPlot_occAbun <- toPlot_occAbun %>% semi_join(seanBiom, by = c("OTU_ID"))
toPlot_occAbun %>% distinct(OTU_ID) %>% nrow()
head(toPlot_occAbun)
toPlot_occAbun %>% filter(is.na(biogeotracesSequence)) %>% nrow()

toPlot_occAbun %>% distinct(OTU_ID) %>% 
  write_csv("occupancyAbundanceNonProNonSynOTUsInSeanCultures_seanAndBiogeotracesClusteredTogether.csv")

#gets distinct samples, OTU, biogeotraces sequences, and biogeotraces samples matches
colnames(toPlot_occAbun)
toPlot_occAbun %>% nrow()
toPlot_occAbun %>% distinct(sample, OTU_ID, biogeotracesSequence, biogeotracesSample) %>% nrow()
toPlot_occAbun <- toPlot_occAbun %>% distinct(sample, OTU_ID, biogeotracesSequence, biogeotracesSample)

#gets the number of matches (number of biogeotraces sequences in the biogeotraces sample) 
#for each biogeotraces sample and OTU
toPlot_occAbun %>% head()
toPlot_occAbun %>% group_by(OTU_ID, biogeotracesSample, biogeotracesSequence) %>% summarize(n = n())
toPlot_occAbun %>% group_by(OTU_ID, biogeotracesSample, biogeotracesSequence, sample) %>% summarize(n = n()) %>% filter(n != 1)
toPlot_occAbun %>% group_by(OTU_ID, biogeotracesSample) %>% distinct(biogeotracesSequence) %>% summarize(n = n())

#how many biogeotraces samples each OTU is found in
toPlot_occAbun %>% group_by(OTU_ID) %>% distinct(biogeotracesSample) %>% summarize(n = n())

#gets the number of biogeotraces sequences for each biogeotraces sample and OTU
toPlot_occAbun <- toPlot_occAbun %>% ungroup()
toPlot_occAbun %>% head()
ubiSummary_occAbun <- toPlot_occAbun %>% group_by(OTU_ID, biogeotracesSample) %>% distinct(biogeotracesSequence) %>% summarize(n = n()) %>% arrange(desc(n))

ubiSummary_occAbun <- ubiSummary_occAbun %>% ungroup()

head(ubiSummary_occAbun)

ubiSummary_occAbun %>% filter(!str_detect(biogeotracesSample, "SRR"))

ubiSummary_occAbun %>% filter(is.na(biogeotracesSample))
ubiSummary_occAbun %>% distinct(OTU_ID) %>% nrow()

occAbunOTU
head(ubiSummary_occAbun)
occAbunOTU %>% distinct(OTU_ID) %>% nrow()
ubiSummary_occAbun %>% distinct(OTU_ID) %>% nrow()

#ubiSummary_occAbun is based on clustering of sean and biogeotraces sequences

#tax is based on clustering of sean and biogeotraces sequences 
tax %>% distinct(`Feature ID`) %>% nrow()

head(ubiSummary_occAbun)
head(tax)

#adds taxonomies to ubiSummary_occAbun
nrow(ubiSummary_occAbun)
ubiSummary_occAbun <- ubiSummary_occAbun %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(ubiSummary_occAbun)

head(ubiSummary_occAbun)
ubiSummary_occAbun %>% filter(is.na(Taxon))



##makes a variable for shortened taxonomy
ubiSummary_occAbun %>% distinct(Taxon)

ubiSummary_occAbun <- ubiSummary_occAbun %>% mutate(newTaxon = str_replace_all(Taxon, "[a-z]__", ""))
ubiSummary_occAbun %>% distinct(newTaxon)

ubiSummary_occAbun <- ubiSummary_occAbun %>% mutate(newTaxon = str_replace_all(newTaxon, ";", ""))
ubiSummary_occAbun %>% distinct(newTaxon)

ubiSummary_occAbun <- ubiSummary_occAbun %>% mutate(newTaxon = str_replace_all(newTaxon, " *$", ""))
ubiSummary_occAbun %>% distinct(newTaxon)

ubiSummary_occAbun <- ubiSummary_occAbun %>% mutate(shortTaxa = str_extract(newTaxon, "[A-Za-z]*$"))
#ubiSummary_occAbun %>% distinct(Taxon, shortTaxa) %>% View()

ubiSummary_occAbun %>% distinct(Taxon, shortTaxa)

#makes a variable for combined OTU ID and shortened taxonomy
ubiSummary_occAbun <- ubiSummary_occAbun %>% mutate(otuSequence = str_c(OTU_ID, shortTaxa, sep = ": "))

ubiSummary_occAbun %>% distinct(otuSequence)
ubiSummary_occAbun %>% distinct(OTU_ID)


#makes a separate dataframe for each OTU
ubiSummary_occAbun_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- ubiSummary_occAbun %>% 
  filter(OTU_ID == "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d")
ubiSummary_occAbun_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f <- ubiSummary_occAbun %>% 
  filter(OTU_ID == "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f")
ubiSummary_occAbun_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- ubiSummary_occAbun %>% 
  filter(OTU_ID == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf")
ubiSummary_occAbun_6e447d49df572a0622adad2e68409bcf687108e8 <- ubiSummary_occAbun %>% 
  filter(OTU_ID == "6e447d49df572a0622adad2e68409bcf687108e8")
ubiSummary_occAbun_3dcd4e930322140ea6d4c075826c2dfbb7124553 <- ubiSummary_occAbun %>% 
  filter(OTU_ID == "3dcd4e930322140ea6d4c075826c2dfbb7124553")
ubiSummary_occAbun_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d <- ubiSummary_occAbun %>% 
  filter(OTU_ID == "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d")
ubiSummary_occAbun_1237c64eb05da167644a4bff9d681bde9b9daa54 <- ubiSummary_occAbun %>% 
  filter(OTU_ID == "1237c64eb05da167644a4bff9d681bde9b9daa54")
ubiSummary_occAbun_08d46233bb26721e2336eb640c2dca0dd94497d3 <- ubiSummary_occAbun %>% 
  filter(OTU_ID == "08d46233bb26721e2336eb640c2dca0dd94497d3")
ubiSummary_occAbun_1618797b285d4650872ef897004644b2990e6d73 <- ubiSummary_occAbun %>% 
  filter(OTU_ID == "1618797b285d4650872ef897004644b2990e6d73")



ubiSummary_occAbun_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d
ubiSummary_occAbun_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f
ubiSummary_occAbun_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf
ubiSummary_occAbun_6e447d49df572a0622adad2e68409bcf687108e8
ubiSummary_occAbun_3dcd4e930322140ea6d4c075826c2dfbb7124553
ubiSummary_occAbun_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d
ubiSummary_occAbun_1237c64eb05da167644a4bff9d681bde9b9daa54
ubiSummary_occAbun_08d46233bb26721e2336eb640c2dca0dd94497d3
ubiSummary_occAbun_1618797b285d4650872ef897004644b2990e6d73


#adds all of the biogeotraces samples to each dataframe
ubiSummary_occAbun_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d <- ubiSummary_occAbun_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_occAbun_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f <- ubiSummary_occAbun_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_occAbun_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf <- ubiSummary_occAbun_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_occAbun_6e447d49df572a0622adad2e68409bcf687108e8 <- ubiSummary_occAbun_6e447d49df572a0622adad2e68409bcf687108e8 %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_occAbun_3dcd4e930322140ea6d4c075826c2dfbb7124553 <- ubiSummary_occAbun_3dcd4e930322140ea6d4c075826c2dfbb7124553 %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_occAbun_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d <- ubiSummary_occAbun_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_occAbun_1237c64eb05da167644a4bff9d681bde9b9daa54 <- ubiSummary_occAbun_1237c64eb05da167644a4bff9d681bde9b9daa54 %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_occAbun_08d46233bb26721e2336eb640c2dca0dd94497d3 <- ubiSummary_occAbun_08d46233bb26721e2336eb640c2dca0dd94497d3 %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_occAbun_1618797b285d4650872ef897004644b2990e6d73 <- ubiSummary_occAbun_1618797b285d4650872ef897004644b2990e6d73 %>% 
  full_join(meta, by = c("biogeotracesSample" = "Run"))


ubiSummary_occAbun %>% distinct(OTU_ID)

#combines all of the dataframes
ubiSummary <- bind_rows(ubiSummary_occAbun_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d,
                        ubiSummary_occAbun_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f,
                        ubiSummary_occAbun_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf,
                        ubiSummary_occAbun_6e447d49df572a0622adad2e68409bcf687108e8,
                        ubiSummary_occAbun_3dcd4e930322140ea6d4c075826c2dfbb7124553,
                        ubiSummary_occAbun_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d,
                        ubiSummary_occAbun_1237c64eb05da167644a4bff9d681bde9b9daa54,
                        ubiSummary_occAbun_08d46233bb26721e2336eb640c2dca0dd94497d3,
                        ubiSummary_occAbun_1618797b285d4650872ef897004644b2990e6d73)

head(ubiSummary)

#the number of biogeotraces samples each OTU is found in
ubiSummary %>% group_by(OTU_ID) %>% summarize(n = n())
ubiSummary %>% group_by(OTU_ID) %>% distinct(biogeotracesSample) %>% summarize(n = n())

nrow(meta)
head(ubiSummary)
ubiSummary %>% filter(is.na(biogeotracesSample)) %>% nrow()
ubiSummary %>% distinct(biogeotracesSample) %>% nrow()
ubiSummary %>% group_by(OTU_ID) %>% distinct(biogeotracesSample) %>% summarize(n = n())

ubiSummary %>% filter(is.na(lat_lon)) %>% nrow()

nrow(numSeq)
numSeq %>% filter(is.na(numSeq))
str(numSeq)
head(numSeq)

#adds number of non-Pro, non-Syn sequences in each biogeotraces sample 
#to ubiSummary
nrow(ubiSummary)
ubiSummary <- ubiSummary %>% left_join(numSeq, by = c("biogeotracesSample"))
nrow(ubiSummary)

head(ubiSummary)
ubiSummary %>% filter(is.na(numSeq)) %>% nrow()

ubiSummary %>% filter(is.na(n)) %>% nrow()

str(ubiSummary$n)
str(ubiSummary$numSeq)

#makes a variable for the number of sequences of each OTU in each biogeotraces 
#sample divided by the total number of non-Pro, non-Syn sequences in the
#biogeotraces sample
ubiSummary <- ubiSummary %>% mutate(propOfBiogeotraces = n/numSeq)

#these are the rows in which the biogeotraces sample did not have the OTU
ubiSummary %>% filter(is.na(propOfBiogeotraces)) %>% nrow()
ubiSummary %>% filter(is.na(propOfBiogeotraces)) %>% distinct(OTU_ID)

ubiSummary %>% arrange(desc(propOfBiogeotraces)) %>% select(propOfBiogeotraces)

#gets rid of unnecessary variables
colnames(ubiSummary)
ubiSummary <- ubiSummary %>% select(OTU_ID, otuSequence, propOfBiogeotraces, biogeotracesSample, Depth, lat_lon)
colnames(ubiSummary)

head(ubiSummary)

##makes lat, long variables

ubiSummary <- ubiSummary %>% mutate(lat = str_extract(lat_lon, "[0-9]*\\.*[0-9]*\\s(N|S)"))
ubiSummary %>% select(lat_lon, lat)
ubiSummary %>% filter(str_detect(lat_lon, "S")) %>% select(lat_lon, lat)
ubiSummary %>% filter(is.na(lat))

ubiSummary <- ubiSummary %>% mutate(long = str_extract(lat_lon, "[0-9]*\\.*[0-9]*\\s(E|W)"))
ubiSummary %>% select(lat_lon, lat, long)
ubiSummary %>% filter(str_detect(lat_lon, "E")) %>% select(lat_lon, lat, long)
ubiSummary %>% filter(is.na(long))

ubiSummary <- ubiSummary %>% mutate(lat_2 = parse_number(lat))
ubiSummary <- ubiSummary %>% mutate(long_2 = parse_number(long))

ubiSummary <- ubiSummary %>% mutate(lat_2 = ifelse(str_detect(lat, "S"), -1*lat_2, lat_2))
ubiSummary <- ubiSummary %>% mutate(long_2 = ifelse(str_detect(long, "W"), -1*long_2, long_2))

ubiSummary %>% select(lat_lon, lat, long, lat_2, long_2)
ubiSummary %>% filter(str_detect(lat_lon, "S")) %>% select(lat_lon, lat, long, lat_2, long_2)
ubiSummary %>% filter(str_detect(lat_lon, "E")) %>% select(lat_lon, lat, long, lat_2, long_2)

#gets rid of unnecessary variables
ubiSummary <- ubiSummary %>% select(-c(lat_lon, lat, long))

#renames variables
colnames(ubiSummary)
colnames(ubiSummary)[6:7]
colnames(ubiSummary)[6:7] <- c("lat", "long")
colnames(ubiSummary)[6:7]

meta %>% distinct(BioSample) %>% nrow()
ubiSummary %>% distinct(biogeotracesSample) %>% nrow()

#these are for the biogeotraces samples that did not have the OTUs
#of interest
ubiSummary %>% filter(is.na(propOfBiogeotraces))
ubiSummary %>% filter(is.na(propOfBiogeotraces)) %>% nrow()

#replaces NA propOfBiogeotraces values with 0
ubiSummary <- ubiSummary %>% mutate(propOfBiogeotraces = ifelse(is.na(propOfBiogeotraces), 0, propOfBiogeotraces))

ubiSummary %>% filter(propOfBiogeotraces == 0) %>% nrow()
ubiSummary %>% filter(is.na(propOfBiogeotraces))

ubiSummary %>% distinct(OTU_ID, otuSequence)

#the number of biogeotraces samples each OTU is found in
ubiSummary %>% group_by(OTU_ID, otuSequence) %>% summarize(n = n())

seanBiom %>% distinct(OTU_ID)
ubiSummary %>% distinct(OTU_ID)
ubiSummary %>% distinct(otuSequence) 

#replaces ": " with ":\n" in otuSequence values
ubiSummary <- ubiSummary %>% mutate(otuSequence = str_replace(otuSequence, ": ", ":\n"))
ubiSummary %>% distinct(otuSequence) 

#replaces NA with "All BioGeoTraces sites" in otuSequence values because I will eventually 
#want to plot these rows as black dots on the map
ubiSummary %>% filter(is.na(otuSequence)) %>% nrow()
ubiSummary <- ubiSummary %>% mutate(otuSequence = ifelse(is.na(otuSequence), "All BioGeoTraces sites", otuSequence))
ubiSummary %>% distinct(otuSequence) 
ubiSummary %>% filter(otuSequence == "All BioGeoTraces sites") %>% nrow()

#the number of biogeotraces samples each OTU is found in
ubiSummary %>% group_by(otuSequence) %>% summarize(n = n())

ubiSummary %>% filter(is.na(propOfBiogeotraces))

ubiSummary %>% distinct(otuSequence)
ubiSummary %>% distinct(otuSequence) %>% arrange(otuSequence)
ubiSummary %>% distinct(otuSequence) %>% nrow()


order <- c("08d46233bb26721e2336eb640c2dca0dd94497d3:\nMuricauda",          
  "1237c64eb05da167644a4bff9d681bde9b9daa54:\nMaricaulis",         
  "1618797b285d4650872ef897004644b2990e6d73:\nThalassospira",      
  "3dcd4e930322140ea6d4c075826c2dfbb7124553:\nHyphomonas",     
  "6e447d49df572a0622adad2e68409bcf687108e8:\nGammaproteobacteria",
  "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d:\nAlteromonadaceae",
  "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf:\nRhodobacteraceae",  
  "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f:\nAlcanivorax",        
  "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d:\nMarinobacter", 
  "All BioGeoTraces sites")

length(order)
ubiSummary %>% distinct(otuSequence) %>% nrow()

#puts otuSequence in order
ubiSummary$otuSequence <- factor(ubiSummary$otuSequence, levels = order)
levels(ubiSummary$otuSequence)

ubiSummary %>% filter(propOfBiogeotraces == 0) %>% distinct(otuSequence)
ubiSummary %>% filter(is.na(propOfBiogeotraces)) %>% distinct(otuSequence)

# project the data
proj="robin"
#long_0 = -90
ubiSummary_2 <- project(ubiSummary, proj)

ubiSummary_2 %>% filter(propOfBiogeotraces != 0) %>% summarize(min = min(log10(propOfBiogeotraces)), max = max(log10(propOfBiogeotraces)))

###is it okay to get rid of the long_0 in the projection??
# plot the world
ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_viridis(direction = -1, limits = c(-5, 0), name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence) +
  labs(title = "OTUs with high mean relative abundance across cultures that\nwere also present in a large proportion of cultures")
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs.svg", dpi = 600, height = 10, width = 20)



ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_gradient(low = "#FEF0D9", high = "#D7301F", limits = c(-5, 0), name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence) +
  labs(title = "OTUs with high mean relative abundance across cultures that\nwere also present in a large proportion of cultures")
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_2.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_2.svg", dpi = 600, height = 10, width = 20)




levels(ubiSummary_2$otuSequence)


ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), 
              ubiSummary_2 %>% filter(otuSequence == "08d46233bb26721e2336eb640c2dca0dd94497d3:\nMuricauda") %>% 
                filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_08d46233bb26721e2336eb640c2dca0dd94497d3_Muricauda.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_08d46233bb26721e2336eb640c2dca0dd94497d3_Muricauda.svg", dpi = 600, height = 10, width = 20)


ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), 
              ubiSummary_2 %>% filter(otuSequence == "1237c64eb05da167644a4bff9d681bde9b9daa54:\nMaricaulis") %>% 
                filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_1237c64eb05da167644a4bff9d681bde9b9daa54_Maricaulis.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_1237c64eb05da167644a4bff9d681bde9b9daa54_Maricaulis.svg", dpi = 600, height = 10, width = 20)


ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), 
    ubiSummary_2 %>% filter(otuSequence == "1618797b285d4650872ef897004644b2990e6d73:\nThalassospira") %>% 
    filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_1618797b285d4650872ef897004644b2990e6d73_Thalassospira.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_1618797b285d4650872ef897004644b2990e6d73_Thalassospira.svg", dpi = 600, height = 10, width = 20)


ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), 
              ubiSummary_2 %>% filter(otuSequence == "3dcd4e930322140ea6d4c075826c2dfbb7124553:\nHyphomonas") %>% 
                filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_3dcd4e930322140ea6d4c075826c2dfbb7124553_Hyphomonas.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_3dcd4e930322140ea6d4c075826c2dfbb7124553_Hyphomonas.svg", dpi = 600, height = 10, width = 20)


ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), 
              ubiSummary_2 %>% filter(otuSequence == "6e447d49df572a0622adad2e68409bcf687108e8:\nGammaproteobacteria") %>% 
                filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_6e447d49df572a0622adad2e68409bcf687108e8_Gammaproteobacteria.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_6e447d49df572a0622adad2e68409bcf687108e8_Gammaproteobacteria.svg", dpi = 600, height = 10, width = 20)



ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), 
              ubiSummary_2 %>% filter(otuSequence == "9ad2160c72d9dbf9ddd62021bce7d5f25f55482d:\nAlteromonadaceae") %>% 
                filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d_Alteromonadaceae.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_9ad2160c72d9dbf9ddd62021bce7d5f25f55482d_Alteromonadaceae.svg", dpi = 600, height = 10, width = 20)


ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), 
              ubiSummary_2 %>% filter(otuSequence == "b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf:\nRhodobacteraceae") %>% 
                filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf_Rhodobacteraceae.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_b7a5ec575bd44a1d1fbf9621ebf819f9166dbcdf_Rhodobacteraceae.svg", dpi = 600, height = 10, width = 20)


ubiSummary %>% filter(is.na(propOfBiogeotraces)) %>% distinct(otuSequence)


ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), 
              ubiSummary_2 %>% filter(otuSequence == "cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f:\nAlcanivorax") %>% 
                filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f_Alcanivorax.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_cd0f80cf49d094a2ef45ad8f7b9949d8b50da23f_Alcanivorax.svg", dpi = 600, height = 10, width = 20)

ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), 
              ubiSummary_2 %>% filter(otuSequence == "f7a384d68d1737a3caf4a4e5b104e2810afb9f9d:\nMarinobacter") %>% 
                filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuSequence)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d_Marinobacter.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithOccupancyAbundanceOTUs_f7a384d68d1737a3caf4a4e5b104e2810afb9f9d_Marinobacter.svg", dpi = 600, height = 10, width = 20)



###ubiquity by class

head(tax)
nrow(tax)

#makes a variable for taxonomic class
tax <- tax %>% mutate(class = str_extract(Taxon, "c__[A-Za-z0-9 \\-\\_\\[\\]]*(?=;)"))
#tax %>% distinct(Taxon, class) %>% View()

#extracts class IDs that weren't extracted correctly before
tax <- tax %>% mutate(class = ifelse(is.na(class), str_extract(Taxon, "c__[A-Za-z0-9 \\-\\_\\[\\]]*"), class))
#tax %>% distinct(Taxon, class) %>% View()

#tax %>% filter(is.na(class)) %>% distinct(Taxon) %>% View()
tax %>% filter(is.na(class)) %>% filter(str_detect(Taxon, "c__"))

#gets rid of "c__" in class
tax %>% distinct(class) %>% head()
tax %>% filter(str_detect(class, ";"))
tax <- tax %>% mutate(class = str_replace(class, "c__", ""))
tax %>% distinct(class) %>% head()
tax %>% filter(str_detect(class, "c__"))

#replaces NA class values with "Unclassified at class level"
tax %>% distinct(class)
tax %>% filter(is.na(class))
tax %>% filter(is.na(class)) %>% nrow()
tax$class <- ifelse(is.na(tax$class), "Unclassified at class level", tax$class)
tax %>% distinct(class)
tax %>% filter(class == "Unclassified at class level") %>% nrow()

#gets rid of unnecessary variable
head(tax)
#tax <- tax %>% select(`Feature ID`,class)
#head(tax)

#there are distinct samples, OTU, biogeotraces sequences, and biogeotraces samples matches
toPlot_forClass %>% head()
toPlot_forClass %>% nrow()
toPlot_forClass %>% distinct(OTU_ID) %>% nrow()
colnames(toPlot_forClass) 
toPlot_forClass %>% distinct(sample, OTU_ID, biogeotracesSequence, biogeotracesSample) %>% nrow()

#adds taxonomy of each OTU to toPlot_forClass
nrow(toPlot_forClass)
toPlot_forClass <- toPlot_forClass %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(toPlot_forClass)

toPlot_forClass %>% filter(is.na(class)) %>% nrow()

toPlot_forClass %>% distinct(sample) %>% nrow()
toPlot_forClass %>% distinct(OTU_ID) %>% nrow()
toPlot_forClass %>% distinct(biogeotracesSample) %>% nrow()

#these are the OTUs that were not in any biogeotraces samples
toPlot_forClass %>% filter(is.na(biogeotracesSample))

#toPlot_forClass %>% distinct(class) %>% arrange(class) %>% View()

toPlot_forClass <- toPlot_forClass %>% ungroup() 


head(toPlot_forClass)
toPlot_forClass %>% distinct(class)

#toPlot_forClass has all of the OTUs
toPlot_forClass %>% distinct(OTU_ID) %>% nrow()

#gets the number of biogeotraces samples each OTU was found in
ubi_class <- toPlot_forClass %>% group_by(class) %>% distinct(sample) %>% summarize(numSamples = n())

ubi_class

#arrange classes by the number of samples they were found in
ubi_class <- ubi_class %>% arrange(desc(numSamples))
ubi_class

head(toPlot_forClass)

#these are the sean culture OTUs that are not present in any biogeotraces samples
toPlot_forClass %>% filter(is.na(biogeotracesSequence))
toPlot_forClass %>% filter(is.na(biogeotracesSequence)) %>% distinct(OTU_ID) %>% nrow()

#10 of the 11 sean culture classes are found in biogeotraces samples
ubi_class %>% nrow()
toPlot_forClass %>% distinct(class) %>% nrow()
toPlot_forClass %>% filter(!is.na(biogeotracesSequence)) %>% distinct(class) %>% nrow()

#[Leptospirae] is the sean culture class that is not present in any of the 
#biogeotraces samples
ubi_class %>% anti_join(toPlot_forClass %>% filter(!is.na(biogeotracesSequence)), by = c("class"))

head(toPlot_forClass)

toPlot_forClass %>% group_by(sample, OTU_ID, biogeotracesSequence, class, biogeotracesSample) %>% summarize(n = n()) %>% filter(n != 1)

#the number of sequences of each class in each biogeotraces sample
toPlot_forClass %>% group_by(biogeotracesSample, class) %>% distinct(biogeotracesSequence) %>% summarize(numSeqs = n())

#the number of OTUs of each class in each biogeotraces sample
toPlot_forClass %>% group_by(biogeotracesSample, class) %>% distinct(OTU_ID) %>% summarize(numOTUs = n())

#calculates the number of biogeotraces sequences of each class in each biogeotraces sample
head(toPlot_forClass) 
toPlot_forClass <- toPlot_forClass %>% ungroup()
ubiSummary_forClass <- toPlot_forClass %>% group_by(biogeotracesSample, class) %>% 
  distinct(biogeotracesSequence) %>% summarize(n = n()) %>% arrange(desc(n))

ubiSummary_forClass <- ubiSummary_forClass %>% ungroup()

##there are only NA biogeotraces sequences when the biogeotraces sample is NA
toPlot_forClass %>% group_by(biogeotracesSample, class) %>% 
  filter(is.na(biogeotracesSequence)) %>% distinct(biogeotracesSequence) %>% 
  summarize(n = n()) %>% arrange(desc(n))

toPlot_forClass %>% filter(is.na(biogeotracesSequence)) %>% distinct(biogeotracesSample)
toPlot_forClass %>% filter(is.na(biogeotracesSample)) %>% distinct(biogeotracesSequence)


#[Leptospirae] is not present in any of the biogeotraces samples
head(ubiSummary_forClass)
ubiSummary_forClass %>% filter(is.na(biogeotracesSample))
ubiSummary_forClass %>% filter(!is.na(biogeotracesSample)) %>% distinct(class)

#excludes rows for which the biogeotraces sample is NA
#these rows are from when the OTU wasn't present in a biogeotraces sample
ubiSummary_forClass <- ubiSummary_forClass %>% ungroup()
ubiSummary_forClass <- ubiSummary_forClass %>% filter(!is.na(biogeotracesSample))
ubiSummary_forClass %>% distinct(class) %>% nrow()

head(ubiSummary_forClass)

#makes a row for Leptospirae because it had no OTUs that were present in
#biogeotraces samples
Leptospirae <- data.frame(biogeotracesSample = NA, class = "[Leptospirae]", n = NA)
Leptospirae
head(ubiSummary_forClass)

#adds the row for Leptospirae to ubiSummary_forClass
nrow(ubiSummary_forClass)
nrow(Leptospirae)
ubiSummary_forClass <- bind_rows(ubiSummary_forClass, Leptospirae)
nrow(ubiSummary_forClass)

ubiSummary_forClass %>% distinct(class)

#makes a dataframe for each class
ubiSummary_Gammaproteobacteria <- ubiSummary_forClass %>% filter(class == "Gammaproteobacteria")
ubiSummary_Alphaproteobacteria <- ubiSummary_forClass %>% filter(class == "Alphaproteobacteria")
ubiSummary_Flavobacteriia <- ubiSummary_forClass %>% filter(class == "Flavobacteriia")
ubiSummary_Betaproteobacteria <- ubiSummary_forClass %>% filter(class == "Betaproteobacteria")
ubiSummary_Unclassified <- ubiSummary_forClass %>% filter(class == "Unclassified at class level")
ubiSummary_Phycisphaerae <- ubiSummary_forClass %>% filter(class == "Phycisphaerae")
ubiSummary_A712011 <- ubiSummary_forClass %>% filter(class == "A712011")
ubiSummary_Cytophagia <- ubiSummary_forClass %>% filter(class == "Cytophagia")
ubiSummary_Rhodothermi <- ubiSummary_forClass %>% filter(class == "[Rhodothermi]")
ubiSummary_Deltaproteobacteria <- ubiSummary_forClass %>% filter(class == "Deltaproteobacteria")
ubiSummary_Leptospirae <- ubiSummary_forClass %>% filter(class == "[Leptospirae]")


#adds all of the biogeotraces samples to each dataframe
ubiSummary_Gammaproteobacteria <- ubiSummary_Gammaproteobacteria %>% full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_Alphaproteobacteria <- ubiSummary_Alphaproteobacteria %>% full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_Flavobacteriia <- ubiSummary_Flavobacteriia %>% full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_Betaproteobacteria <- ubiSummary_Betaproteobacteria %>% full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_Unclassified <- ubiSummary_Unclassified %>% full_join(meta, by = c("biogeotracesSample" = "Run")) 
ubiSummary_Phycisphaerae <- ubiSummary_Phycisphaerae %>% full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_A712011 <- ubiSummary_A712011 %>% full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_Cytophagia <- ubiSummary_Cytophagia %>% full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_Rhodothermi <- ubiSummary_Rhodothermi %>% full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_Deltaproteobacteria <- ubiSummary_Deltaproteobacteria %>% full_join(meta, by = c("biogeotracesSample" = "Run"))
ubiSummary_Leptospirae <- ubiSummary_Leptospirae %>% full_join(meta, by = c("biogeotracesSample" = "Run"))

#makes a variable for class in each dataframe
ubiSummary_Gammaproteobacteria <- ubiSummary_Gammaproteobacteria %>% mutate(class = "Gammaproteobacteria")
ubiSummary_Alphaproteobacteria <- ubiSummary_Alphaproteobacteria %>% mutate(class = "Alphaproteobacteria")
ubiSummary_Flavobacteriia <- ubiSummary_Flavobacteriia %>% mutate(class = "Flavobacteriia")
ubiSummary_Betaproteobacteria <- ubiSummary_Betaproteobacteria %>% mutate(class = "Betaproteobacteria")
ubiSummary_Unclassified <- ubiSummary_Unclassified %>% mutate(class = "Unclassified at class level")
ubiSummary_Phycisphaerae <- ubiSummary_Phycisphaerae %>% mutate(class = "Phycisphaerae") 
ubiSummary_A712011 <- ubiSummary_A712011 %>% mutate(class = "A712011")
ubiSummary_Cytophagia <- ubiSummary_Cytophagia %>% mutate(class = "Cytophagia")
ubiSummary_Rhodothermi <- ubiSummary_Rhodothermi %>% mutate(class ="[Rhodothermi]")
ubiSummary_Deltaproteobacteria <- ubiSummary_Deltaproteobacteria %>% mutate(class = "Deltaproteobacteria")
ubiSummary_Leptospirae <- ubiSummary_Leptospirae %>% mutate(class = "[Leptospirae]") 

ubiSummary_Gammaproteobacteria %>% nrow()
ubiSummary_Alphaproteobacteria %>% nrow()
ubiSummary_Flavobacteriia %>% nrow()
ubiSummary_Betaproteobacteria %>% nrow()
ubiSummary_Unclassified %>% nrow()
ubiSummary_Phycisphaerae %>% nrow()
ubiSummary_A712011 %>% nrow()
ubiSummary_Cytophagia %>% nrow()
ubiSummary_Rhodothermi %>% nrow()
ubiSummary_Deltaproteobacteria %>% nrow()
ubiSummary_Leptospirae %>% nrow()

#gets rid of the extra, unnecessary row in ubiSummary_Leptospirae
head(ubiSummary_Leptospirae)
ubiSummary_Leptospirae %>% filter(is.na(biogeotracesSample))
ubiSummary_Leptospirae <- ubiSummary_Leptospirae %>% filter(!is.na(biogeotracesSample))

ubiSummary_Leptospirae %>% nrow()

#binds the dataframes of each class together
ubiSummary <- bind_rows(ubiSummary_Gammaproteobacteria,
                        ubiSummary_Alphaproteobacteria,
                        ubiSummary_Flavobacteriia,
                        ubiSummary_Betaproteobacteria,
                        ubiSummary_Unclassified,
                        ubiSummary_Phycisphaerae,
                        ubiSummary_A712011,
                        ubiSummary_Cytophagia,
                        ubiSummary_Rhodothermi,
                        ubiSummary_Deltaproteobacteria,
                        ubiSummary_Leptospirae)

ubiSummary %>% filter(is.na(biogeotracesSample)) %>% nrow()

#each class has all of the biogeotraces samples
ubiSummary %>% head()
ubiSummary %>% group_by(class) %>% summarize(n = n()) %>% filter(n != 610)
ubiSummary %>% group_by(class) %>% distinct(biogeotracesSample) %>% summarize(n = n()) %>% filter(n != 610)

nrow(meta)
ubiSummary %>% filter(is.na(biogeotracesSample)) %>% nrow()
ubiSummary %>% distinct(biogeotracesSample) %>% nrow()

ubiSummary %>% filter(is.na(lat_lon)) %>% nrow()

#there is only one row for each biogeotraces sample, class pair
head(ubiSummary)
ubiSummary %>% group_by(biogeotracesSample, class) %>% summarize(n = n()) %>% filter(n != 1)

#numSeq has the number of non-Pro, non-Syn sequences in each biogeotraces sample
head(numSeq)
nrow(numSeq)

#adds number of sequence in each biogeotraces sample to ubiSummary
nrow(ubiSummary)
ubiSummary <- ubiSummary %>% left_join(numSeq, by = c("biogeotracesSample"))
nrow(ubiSummary)

ubiSummary %>% filter(is.na(numSeq)) %>% nrow()

#makes a variable for the proportion of non-Pro, non-Syn sequences that each 
#class makes up in each biogeotraces sample
#n is the number of biogeotraces sequences of each class
#numSeq is the total number of non-Pro, non-Syn sequences in each biogeotraces sample
head(ubiSummary)
str(ubiSummary$n)
str(ubiSummary$numSeq)
ubiSummary <- ubiSummary %>% mutate(propOfBiogeotraces = n/numSeq)

#there are NA props when the class is not present in the sample
ubiSummary %>% filter(is.na(propOfBiogeotraces))

ubiSummary %>% arrange(desc(propOfBiogeotraces)) %>% select(propOfBiogeotraces)

#gets rid of unnecessary variables
head(ubiSummary)
colnames(ubiSummary)
ubiSummary <- ubiSummary %>% select(class, propOfBiogeotraces, biogeotracesSample, Depth, lat_lon)
head(ubiSummary)

#all of the rows are distinct
nrow(ubiSummary)
ubiSummary %>% distinct() %>% nrow()

##makes variables for lat, long

#lat
ubiSummary %>% select(lat_lon) %>% head()
ubiSummary <- ubiSummary %>% mutate(lat = str_extract(lat_lon, "[0-9]*\\.*[0-9]*\\s(N|S)"))
ubiSummary %>% select(lat_lon, lat) %>% head()
ubiSummary %>% filter(str_detect(lat_lon, "S")) %>% select(lat_lon, lat) %>% head()
ubiSummary %>% filter(is.na(lat))

#long
ubiSummary <- ubiSummary %>% mutate(long = str_extract(lat_lon, "[0-9]*\\.*[0-9]*\\s(E|W)"))
ubiSummary %>% select(lat_lon, lat, long) %>% head()
ubiSummary %>% filter(str_detect(lat_lon, "E")) %>% select(lat_lon, lat, long) %>% head()
ubiSummary %>% filter(is.na(long))

#makes lat, long values numeric
ubiSummary <- ubiSummary %>% mutate(lat_2 = parse_number(lat))
ubiSummary <- ubiSummary %>% mutate(long_2 = parse_number(long))

ubiSummary <- ubiSummary %>% mutate(lat_2 = ifelse(str_detect(lat, "S"), -1*lat_2, lat_2))
ubiSummary <- ubiSummary %>% mutate(long_2 = ifelse(str_detect(long, "W"), -1*long_2, long_2))

ubiSummary %>% select(lat_lon, lat, long, lat_2, long_2) %>% head()
ubiSummary %>% distinct(lat_lon, lat, long, lat_2, long_2) %>% head()

#gets rid of unnecessary variables
ubiSummary <- ubiSummary %>% select(-c(lat_lon, lat, long))

colnames(ubiSummary)
colnames(ubiSummary)[5:6]
colnames(ubiSummary)[5:6] <- c("lat", "long")

ubiSummary %>% filter(is.na(propOfBiogeotraces)) %>% nrow()

#replaces NA propOfBiogeotraces values with 0
ubiSummary <- ubiSummary %>% mutate(propOfBiogeotraces = ifelse(is.na(propOfBiogeotraces), 0, propOfBiogeotraces))

ubiSummary %>% filter(propOfBiogeotraces == 0) %>% nrow()

ubiSummary %>% distinct(class)
ubiSummary %>% head()

ubiSummary %>% group_by(class) %>% summarize(n = n()) %>% filter(n != 610)

ubiSummary %>% filter(is.na(class))

head(ubiSummary)

#makes a dataframe for all of the biogeotraces samples
allSamples <- ubiSummary %>% 
  filter(propOfBiogeotraces == 0) %>% distinct(biogeotracesSample, long, lat, propOfBiogeotraces)

allSamples %>% nrow()
allSamples %>% head()

#makes the class "All BioGeoTraces sites" because I want this to be the
#plot title
allSamples <- allSamples %>% mutate(class = "All BioGeoTraces sites")

allSamples %>% head()

#adds allSamples to ubiSummary
nrow(ubiSummary)
nrow(allSamples)
ubiSummary <- bind_rows(ubiSummary, allSamples)
nrow(ubiSummary)

head(ubiSummary)

#puts the classes in order by how many sean culture samples they were 
#found in
ubi_class
ubi_class %>% arrange(desc(numSamples))
order <- ubi_class %>% arrange(desc(numSamples)) %>% select(class)
order
order <- order$class
order
length(order)
order[12] <- "All BioGeoTraces sites"
order
length(order)
ubiSummary %>% distinct(class) %>% nrow()
head(ubiSummary)
ubiSummary$class <- factor(ubiSummary$class, levels = order)
levels(ubiSummary$class)
levels(ubiSummary$class) %>% length()
ubi_class %>% arrange(desc(numSamples))

ubiSummary %>% distinct(class)
ubiSummary %>% group_by(class) %>% summarize(n = n()) %>% filter(n != 610)

ubiSummary %>% filter(is.na(propOfBiogeotraces)) %>% distinct(class)
ubiSummary %>% filter(propOfBiogeotraces == 0) %>% distinct(class)

ubiSummary %>% distinct(class)

#[Leptospirae] was not found in any biogeotraces samples
ubiSummary %>% filter(class == "[Leptospirae]") %>% distinct(propOfBiogeotraces)

ubiSummary %>% head()
ubiSummary %>% filter(propOfBiogeotraces == 0) %>% head()
ubiSummary %>% filter(propOfBiogeotraces == 0) %>% distinct(biogeotracesSample, long, lat)
ubiSummary %>% filter(propOfBiogeotraces == 0) %>% distinct(biogeotracesSample, long, lat) %>% nrow()

levels(ubiSummary$class)


library(ggworldmap)

# project the data
proj="robin"
#long_0 = -90
ubiSummary_2 <- project(ubiSummary, proj)

###is it okay to get rid of the long_0 in the projection??
# plot the world
ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "[Leptospirae]"), size=2, alpha = 0) +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "All BioGeoTraces sites"), size=1.5, color = "black") +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class) +
  labs(title = "Taxonomic classes present in cultures,\narranged by number of cultures they were found in")
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses.svg", dpi = 600, height = 10, width = 20)


levels(ubiSummary_2$class)

ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "Alphaproteobacteria") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "Alphaproteobacteria") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Alphaproteobacteria.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Alphaproteobacteria.svg", dpi = 600, height = 5, width = 10)


ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "Gammaproteobacteria") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "Gammaproteobacteria") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Gammaproteobacteria.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Gammaproteobacteria.svg", dpi = 600, height = 5, width = 10)



ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "Flavobacteriia") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "Flavobacteriia") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Flavobacteriia.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Flavobacteriia.svg", dpi = 600, height = 5, width = 10)




ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "Cytophagia") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "Cytophagia") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Cytophagia.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Cytophagia.svg", dpi = 600, height = 5, width = 10)



ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "[Rhodothermi]") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "[Rhodothermi]") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Rhodothermi.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Rhodothermi.svg", dpi = 600, height = 5, width = 10)



ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "Phycisphaerae") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "Phycisphaerae") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Phycisphaerae.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Phycisphaerae.svg", dpi = 600, height = 5, width = 10)



ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "Betaproteobacteria") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "Betaproteobacteria") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Betaproteobacteria.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Betaproteobacteria.svg", dpi = 600, height = 5, width = 10)




ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "Unclassified at class level") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "Unclassified at class level") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Unclassified.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Unclassified.svg", dpi = 600, height = 5, width = 10)





ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "[Leptospirae]") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "[Leptospirae]") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Leptospirae.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Leptospirae.svg", dpi = 600, height = 5, width = 10)



ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "A712011") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "A712011") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_A712011.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_A712011.svg", dpi = 600, height = 5, width = 10)


ggplot() +
  #geom_point(aes(long, lat, color=propOfBiogeotraces), ubiSummary_2 %>% filter(class == "Deltaproteobacteria") %>% filter(propOfBiogeotraces == 0), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), ubiSummary_2 %>% filter(class == "Deltaproteobacteria") %>% filter(propOfBiogeotraces != 0), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_distiller(palette="Spectral", name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~class)
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Deltaproteobacteria.png", dpi = 600, height = 5, width = 10)
ggsave("seanSeqsNoProSyn_biogeotracesSamplesWithUbiquitousClasses_Deltaproteobacteria.svg", dpi = 600, height = 5, width = 10)





