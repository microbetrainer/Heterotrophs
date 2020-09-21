library(tidyverse)
library(lubridate)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads the sean and biogeotraces sequences present in each cluster
biomG <- read_csv("combinedBiogeotracesBiom.csv")

#28022 clusters 
#1184489 sequences
#which is the correct number of clusters and sequences
head(biomG)
biomG %>% distinct(OTU_ID) %>% nrow()
biomG %>% distinct(sample) %>% nrow()

str(biomG)
biomG %>% distinct(present)

#each sequence is only in one cluster
head(biomG)
biomG %>% filter(present == 1) %>% group_by(sample) %>% summarize(n = n()) %>% filter(n != 1)
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
head(seqsNoProSyn)
seqsNoProSyn %>% distinct(V1) %>% nrow()

head(seqG)
head(seqsNoProSyn)

#excludes Pro and Syn sequences from seqG
seqG <- seqG %>% semi_join(seqsNoProSyn, by = c("sequence" = "V1"))

#correct number of sequences and samples
#seqG no longer includes pro and syn
seqG %>% distinct(sequence) %>% nrow()
seqG %>% distinct(sample) %>% nrow()

#biomG includes Pro and Syn sequences in sean samples
head(biomG)
biomG %>% filter(str_detect(sample, "sean")) %>% nrow()

##makes the sean sequence IDs in biomG match seqG
biomG %>% filter(str_detect(sample, "sean")) %>% head()
biomG <- biomG %>% mutate(sample = ifelse(str_detect(sample, "sean"), str_replace(sample, "X", ""), sample))

biomG %>% filter(str_detect(sample, "sean")) %>% head()
biomG <- biomG %>% mutate(sample = ifelse(str_detect(sample, "sean"), str_c("contig_", sample), sample))

biomG %>% filter(str_detect(sample, "sean")) %>% head()
biomG <- biomG %>% mutate(sample = ifelse(str_detect(sample, "sean"), str_replace(sample, "sean", ""), sample))

#these are the sean sequences
biomG %>% filter(str_detect(sample, "contig")) %>% head()
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

#there are clusters in sean samples that have multiple sequences across 
#sean samples
head(biomG)
biomG %>% filter(str_detect(sample, "contig"))
biomG %>% filter(str_detect(sample, "contig")) %>% group_by(OTU_ID) %>% summarize(n = n()) %>% arrange(desc(n))

#gets just the clusters that sean sequences are in 
head(biomG)
seanBiom <- biomG %>% filter(str_detect(sample, "contig"))
head(seanBiom)
seanBiom %>% distinct(sample) %>% nrow()

#gets just the clusters that biogeotraces sequences are in
head(biomG)
biogeoBiom <- biomG %>% filter(!str_detect(sample, "contig"))

nrow(biomG)
nrow(seanBiom)
nrow(biogeoBiom)
nrow(seanBiom)+nrow(biogeoBiom)==nrow(biomG)

head(seqG)
head(seanBiom)

#doesn't include pro and syn sequences
seqG %>% distinct(sequence) %>% nrow()

#includes pro and syn sequences
seanBiom %>% distinct(sample) %>% nrow()

head(seqG)
head(seanBiom)
seanBiom %>% distinct(present)

#adds the cluster each sean sequence is assigned to to the abundance of
#sean sequences across sean samples
#with Pro and Syn
nrow(seqG)
merged <- seqG %>% left_join(seanBiom, by = c("sequence" = "sample"))
nrow(merged)

merged %>% filter(is.na(present))

#gets rid of unnecessary variable
merged <- merged %>% select(-present)

head(merged)

#the one OTU missing from merged I believe is Pro, Syn
merged %>% distinct(OTU_ID) %>% nrow()
seanBiom %>% distinct(OTU_ID) %>% nrow()

#doesn't include pro and syn
merged %>% distinct(sequence) %>% nrow()

#correct number of samples
merged %>% distinct(sample) %>% nrow()

head(merged)

head(biogeoBiom)
biogeoBiom %>% distinct(present)

#gets rid of unnecessary variable
biogeoBiom <- biogeoBiom %>% select(-present)

#there are 28008 OTUs present in biogeotraces
biogeoBiom %>% distinct(OTU_ID) %>% nrow()

#there are 27925 OTUs present in biogeotraces that are not present in sean samples
head(biogeoBiom)
head(merged)
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

#renames variables
colnames(biogeoBiom)[2]
colnames(biogeoBiom)[2] <- "sequence"

#extracts the biogeotraces sample from the biogeotraces sequence ID
head(biogeoBiom)
biogeoBiom <- biogeoBiom %>% mutate(biogeotracesSample = str_extract(sequence, "[a-zA-Z0-9]*"))
head(biogeoBiom)
biogeoBiom %>% filter(is.na(biogeotracesSample))
biogeoBiom %>% filter(!str_detect(biogeotracesSample, "SRR"))

#all of the sequences are only found in one OTU, biogeotraces sample
biogeoBiom %>% group_by(sequence) %>% distinct() %>% summarize(n = n()) %>% filter(n != 1)

#gets the number of sequences of each OTU across all samples for biogeotraces
head(biogeoBiom)
biogeo_summary <- biogeoBiom %>% group_by(OTU_ID) %>% distinct(sequence) %>% summarize(numSeqs = n())
head(biogeo_summary)
nrow(biogeo_summary)
biogeo_summary %>% summarize(totalSeqs = sum(numSeqs))

#loads taxonomies of OTUs
tax <- read_tsv("sequencesForTree_noContamination_noMock_withProSyn_no1223_biogeotracesTaxonomy.tsv")

#gets rid of first row because it just has comments
tax[1,]
tax <- tax[-1,]

#gets rid of unnecessary variable
head(tax)
tax <- tax %>% select(1,2)
nrow(tax)

head(biogeo_summary)
head(tax)

#adds taxonomy to number of sequences in each OTU in biogeotraces
nrow(biogeo_summary)
biogeo_summary <- biogeo_summary %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(biogeo_summary)

head(biogeo_summary)

biogeo_summary %>% filter(is.na(Taxon))

#I haven't excluded the pro or syn sequences in biogeo_summary yet
biogeo_summary %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)
biogeo_summary %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(Taxon)

#excludes Pro and Syn OTUs from biogeo_summary
biogeo_summary %>% distinct(Taxon) %>% nrow()
biogeo_summary <- biogeo_summary %>% filter(!str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc"))
biogeo_summary %>% distinct(Taxon) %>% nrow()

biogeo_summary %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)

#biogeo_summary %>% arrange(desc(numSeqs)) %>% View()
str(biogeo_summary)
str(biogeo_summary$numSeqs)
biogeo_summary %>% arrange(desc(numSeqs)) %>% head()

#gets the OTUs that have the most sequences across all biogeotraces samples
mostCommonBiogeoOTUs <- biogeo_summary %>% arrange(desc(numSeqs)) %>% slice(1:10)

mostCommonBiogeoOTUs
nrow(mostCommonBiogeoOTUs)

#merged has the OTUs, sequences, and their abundances across each sample 
#doesn't include Pro or Syn
head(merged)
merged %>% distinct(sample) %>% nrow()
merged %>% distinct(sequence) %>% nrow()

head(mostCommonBiogeoOTUs)
head(merged)

#none of the most common biogeotraces OTUs are found in sean cultures
mostCommonBiogeoOTUs %>% semi_join(merged, by = c("OTU_ID"))

#these are the biogeotraces clusters that are also in sean cultures arranged 
#by their total number of sequences across biogeotraces samples
biogeo_summary %>% semi_join(merged, by = c("OTU_ID")) %>% nrow()
biogeo_summary %>% semi_join(merged, by = c("OTU_ID")) %>% str()
biogeo_summary %>% semi_join(merged, by = c("OTU_ID")) %>% arrange(desc(numSeqs))

#gets the OTUs with the top 10 highest total number of sequences across biogeotraces sequences
#that are also in sean culture samples
biogeo_summary %>% semi_join(merged, by = c("OTU_ID")) %>% str()
mostCommonBiogeoOTUs_alsoInCulture <- biogeo_summary %>% semi_join(merged, by = c("OTU_ID")) %>% 
  arrange(desc(numSeqs)) %>% slice(1:10)

mostCommonBiogeoOTUs_alsoInCulture
nrow(mostCommonBiogeoOTUs_alsoInCulture)

mostCommonBiogeoOTUs %>% semi_join(mostCommonBiogeoOTUs_alsoInCulture, by = c("OTU_ID"))

#makes variable for whether OTU is in sean culture or not
mostCommonBiogeoOTUs <- mostCommonBiogeoOTUs %>% mutate(source = "justBiogeotraces")
mostCommonBiogeoOTUs_alsoInCulture <- mostCommonBiogeoOTUs_alsoInCulture %>% mutate(source = "biogeotracesAndCultures")

mostCommonBiogeoOTUs
mostCommonBiogeoOTUs_alsoInCulture

#combines the two most common biogeotraces OTU dataframes
comb <- bind_rows(mostCommonBiogeoOTUs, mostCommonBiogeoOTUs_alsoInCulture) 

head(comb)
comb %>% nrow()
comb %>% distinct(OTU_ID) %>% nrow()

head(comb)
head(biogeoBiom)

#biogeoBiom still includes pro and syn sequences
biogeoBiom %>% distinct(biogeotracesSample) %>% nrow()
biogeoBiom %>% distinct(sequence) %>% nrow()
biogeoBiom %>% distinct(OTU_ID) %>% nrow()

head(comb)
head(biogeoBiom)

#adds which samples the biogeotraces sequences of each OTU are present in to comb
nrow(comb)
comb <- comb %>% left_join(biogeoBiom, by = c("OTU_ID"))
nrow(comb)

head(comb)
comb %>% filter(is.na(biogeotracesSample))

#correct number of OTUs
comb %>% distinct(OTU_ID) %>% nrow()
comb %>% distinct(OTU_ID, numSeqs) %>% nrow()

#numSeqs is the total number of sequences across all biogeotraces samples that the OTU
#is found in
#so makes this variable name more explanatory
head(comb)
colnames(comb)[2]
colnames(comb)[2] <- "totalNumSeqsAcrossBiogeotraces"

comb %>% head()

#gets the total number of sequences of each OTU found in each biogeotraces sample
comb <- comb %>% ungroup()
comb <- comb %>% group_by(OTU_ID, Taxon, totalNumSeqsAcrossBiogeotraces, source, biogeotracesSample) %>% 
  distinct(sequence) %>% summarize(numSeqs = n())

head(comb)
comb %>% nrow()
comb %>% distinct(OTU_ID, Taxon, totalNumSeqsAcrossBiogeotraces, source) %>% nrow()

head(comb)
comb %>% filter(is.na(biogeotracesSample))


##loads metadata for biogeotraces samples
#from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP109831&o=acc_s%3Aa&s=SRR6507277,SRR6507278,SRR6507279,SRR6507280,SRR5720219,SRR5720220,SRR5720221,SRR5720222,SRR5720223,SRR5720224,SRR5720225,SRR5720226,SRR5720227,SRR5720228,SRR5720229,SRR5720230,SRR5720231,SRR5720232,SRR5720233,SRR5720234,SRR5720235,SRR5720236,SRR5720237,SRR5720238,SRR5720239,SRR5720240,SRR5720241,SRR5720242,SRR5720243,SRR5720244,SRR5720245,SRR5720246,SRR5720247,SRR5720248,SRR5720249,SRR5720250,SRR5720251,SRR5720252,SRR5720253,SRR5720254,SRR5720255,SRR5720256,SRR5720257,SRR5720258,SRR5720259,SRR5720260,SRR5720261,SRR5720262,SRR5720263,SRR5720264,SRR5720265,SRR5720266,SRR5720267,SRR5720268,SRR5720269,SRR5720270,SRR5720271,SRR5720272,SRR5720273,SRR5720274,SRR5720275,SRR5720276,SRR5720277,SRR5720278,SRR5720279,SRR5720280,SRR5720281,SRR5720282,SRR5720283,SRR5720284,SRR5720285,SRR5720286,SRR5720287,SRR5720288,SRR5720289,SRR5720290,SRR5720291,SRR5720292,SRR5720293,SRR5720294,SRR5720295,SRR5720296,SRR5720297,SRR5720298,SRR5720299,SRR5720300,SRR5720301,SRR5720302,SRR5720303,SRR5720304,SRR5720305,SRR5720306,SRR5720307,SRR5720308,SRR5720309,SRR5720310,SRR5720311,SRR5720312,SRR5720313,SRR5720314,SRR5720315,SRR5720316,SRR5720317,SRR5720318,SRR5720319,SRR5720320,SRR5720321,SRR5720322,SRR5720323,SRR5720324,SRR5720325,SRR5720326,SRR5720327,SRR5720328,SRR5720329,SRR5720330,SRR5720331,SRR5720332,SRR5720333,SRR5720334,SRR5720335,SRR5720336,SRR5720337,SRR5720338,SRR5720339,SRR5720340,SRR5720341,SRR5720342,SRR5720343,SRR5720344
meta1 <- read_csv("SraRunTable.csv")

#from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP110813
meta2 <- read_csv("SraRunTable2.csv")

str(meta1)
str(meta2)

#cleans up metadata bottle ID and collection date
meta1$bottle_id <- as.character(meta1$bottle_id)
meta1$collection_date <- dmy(meta1$collection_date)

#combines meta1 and meta2
meta <- bind_rows(meta1, meta2)

str(meta)

#there are 610 rows in meta which makes sense because this is the number of biogeotraces 
#samples
#there is only one row for each sample in the metadata
head(meta)
nrow(meta)
meta %>% group_by(Run) %>% summarize(n = n()) %>% filter(n > 1)

#all of the biogeotraces samples have corresponding metadata in meta
#all of the biogeotraces samples are already present in comb because in total, the 
#20 OTUs are present in all of the biogeotraces samples
head(comb)
head(meta)
comb %>% anti_join(meta, by = c("biogeotracesSample" = "Run")) %>% distinct(biogeotracesSample)
meta %>% anti_join(comb, by = c("Run" = "biogeotracesSample")) %>% distinct(Run)

comb <- comb %>% ungroup() 

comb %>% distinct(biogeotracesSample) %>% nrow()

#adds metadata of each biogeotraces data to comb
nrow(comb)
comb <- comb %>% full_join(meta, by = c("biogeotracesSample" = "Run"))
nrow(comb)

head(comb)
comb %>% filter(is.na(lat_lon)) %>% nrow()

#makes depth variable numeric
str(comb$Depth)
comb$Depth <- parse_number(comb$Depth)
str(comb$Depth)

comb %>% distinct(biogeotracesSample) %>% nrow()

#how many biogeotraces samples each OTU is found in
head(comb)
comb %>% group_by(OTU_ID) %>% summarize(n = n())
comb %>% group_by(OTU_ID) %>% distinct(biogeotracesSample) %>% summarize(n = n())


#loads number of sequences in each biogeotraces sample
#from ubiquityAnalysisSeanBiogeotraces_clustered.R (need to check)
numSeq <- read_csv("numSequencesInEachBiogeotraceSampleNoProSyn.csv")
head(numSeq)
colnames(numSeq)[2]
colnames(numSeq)[2] <- "numSeqBiogeotracesSample"
head(numSeq)

nrow(numSeq)
numSeq %>% filter(is.na(numSeqBiogeotracesSample))
str(numSeq)
str(numSeq$numSeqBiogeotracesSample)

head(numSeq)
colnames(comb)

#adds the total number of non-Pro, non-Syn sequences in each biogeotraces 
#sample ot comb
nrow(comb)
comb <- comb %>% left_join(numSeq, by = c("biogeotracesSample"))
nrow(comb)

comb %>% filter(is.na(numSeqBiogeotracesSample)) %>% nrow()

head(comb)

#numSeqs is the number of sequences of the OTU in the biogeotraces sample
#the number of sequences of each OTU in a biogeotraces sample is never NA because 
#I didn't join each biogeotraces sample to each OTU
comb %>% filter(is.na(numSeqs)) %>% nrow()

str(comb$numSeqs)
str(comb$numSeqBiogeotracesSample)

#makes a variable for the number of sequences of each OTU in each biogeotraces 
#sample divided by the total number of sequences, excluding Pro and Syn in 
#the biogeotraces sample
comb <- comb %>% mutate(propOfBiogeotraces = numSeqs/numSeqBiogeotracesSample)

comb %>% filter(is.na(propOfBiogeotraces)) %>% nrow()





##makes cleaned-up lat, long variables
comb %>% select(lat_lon)
comb <- comb %>% mutate(lat = str_extract(lat_lon, "[0-9]*\\.*[0-9]*\\s(N|S)"))
comb %>% select(lat_lon, lat)

comb %>% filter(is.na(lat)) %>% nrow()
comb %>% filter(!str_detect(lat, "N|S")) %>% nrow()

comb <- comb %>% mutate(long = str_extract(lat_lon, "[0-9]*\\.*[0-9]*\\s(E|W)"))
comb %>% select(lat_lon, lat, long)

comb %>% filter(is.na(long)) %>% nrow()
comb %>% filter(!str_detect(long, "E|W")) %>% nrow()

comb <- comb %>% mutate(lat_2 = parse_number(lat))
comb <- comb %>% mutate(long_2 = parse_number(long))
comb %>% select(lat_lon, lat, long, lat_2, long_2)

comb <- comb %>% mutate(lat_2 = ifelse(str_detect(lat, "S"), -1*lat_2, lat_2))
comb <- comb %>% mutate(long_2 = ifelse(str_detect(long, "W"), -1*long_2, long_2))

comb %>% select(lat_lon, lat_2, long_2)
comb %>% distinct(lat_lon, lat_2, long_2)

#gets rid of unnecessary variables
comb <- comb %>% select(-c(lat_lon, lat, long))

colnames(comb)
colnames(comb)[46:47]
colnames(comb)[46:47] <- c("lat", "long")

str(comb$lat)
str(comb$long)

head(meta)

meta %>% distinct(Run, BioSample) %>% nrow()
comb %>% distinct(biogeotracesSample) %>% nrow()

head(comb)
comb %>% filter(is.na(propOfBiogeotraces)) %>% nrow()
str(comb$propOfBiogeotraces)

comb %>% distinct(OTU_ID) %>% nrow()

#the number of biogeotraces samples each OTU is found in
head(comb)
comb %>% group_by(OTU_ID) %>% summarize(n = n())

comb %>% distinct(OTU_ID, Taxon) %>% nrow()
comb %>% distinct(Taxon) %>% nrow()

comb %>% distinct(Taxon)

#gets rid of k__, p__, c__ ... in taxonomies 
comb <- comb %>% mutate(Taxon2 = str_replace_all(Taxon, "[a-z]__", ""))
comb %>% distinct(Taxon2)

#gets rid of ; in taxonomies 
comb <- comb %>% mutate(Taxon2 = str_replace_all(Taxon2, ";", ""))
comb %>% distinct(Taxon2)

#gets rid of spaces at the end of taxonomies
comb <- comb %>% mutate(Taxon2 = str_replace_all(Taxon2, " *$", ""))
comb %>% distinct(Taxon2)

#makes a variable for the shortened taxonomy
comb <- comb %>% mutate(shortTaxa = str_extract(Taxon2, "((Candidatus )|(HTCC2188 )|(Pelagibacter ))*[A-Za-z0-9]*$"))

comb %>% distinct(OTU_ID) %>% nrow()
#comb %>% distinct(Taxon, shortTaxa) %>% View()

#makes a variable with the combined OTU and shortTaxa
head(comb)
colnames(comb)
comb <- comb %>% mutate(otuTax = str_c(OTU_ID, shortTaxa, sep = ":\n"))
comb %>% distinct(otuTax)
comb %>% distinct(otuTax) %>% nrow()

comb %>% filter(is.na(propOfBiogeotraces)) %>% nrow()
comb %>% filter(propOfBiogeotraces == 0) %>% nrow()

colnames(comb)

#makes a dataframe with the lat, long of each biogeotraces sample
allSites <- comb %>% distinct(biogeotracesSample, lat, long)
head(allSites)
nrow(allSites)
allSites <- allSites %>% mutate(otuTax = "All BioGeoTraces sites")
head(allSites)
allSites <- allSites %>% mutate(source = "allSites")
head(allSites)
nrow(allSites)

#combines dataframe with the lat, long of each biogeotraces sample to comb
#so that I can plot a point for each biogeotraces sample
nrow(comb)
nrow(allSites)
comb <- bind_rows(comb, allSites)
nrow(comb)
7866+nrow(allSites)==nrow(comb)

comb %>% filter(is.na(propOfBiogeotraces)) %>% distinct(source)

#puts the OTUs in order by their total number of sequences across all
#biogeotraces samples
colnames(comb)
str(comb$totalNumSeqsAcrossBiogeotraces)
order <- comb %>% arrange(desc(totalNumSeqsAcrossBiogeotraces)) %>% distinct(otuTax)
order
order <- order$otuTax %>% as.list()
order
length(order)
head(comb)
comb$otuTax <- factor(comb$otuTax, levels = order)
levels(comb$otuTax)
comb %>% arrange(desc(totalNumSeqsAcrossBiogeotraces)) %>% distinct(otuTax, totalNumSeqsAcrossBiogeotraces)


library(ggworldmap)

# project the data
proj="robin"
#long_0 = -90
comb_2 <- project(comb, proj)

comb_2 %>% distinct(source)
comb_2 %>% filter(source == "biogeotracesAndCultures" | source == "allSites") %>% distinct(source)
comb_2 %>% filter(source == "biogeotracesAndCultures" | source == "allSites") %>%
  filter(propOfBiogeotraces != 0 & !is.na(propOfBiogeotraces)) %>% distinct(source)

comb_2 %>% filter(source == "allSites") %>%
  distinct(propOfBiogeotraces)

comb_2 %>% filter(source == "biogeotracesAndCultures" | source == "justBiogeotraces") %>%
  filter(propOfBiogeotraces == 0 | is.na(propOfBiogeotraces)) %>% nrow()

comb_2 %>% distinct(otuTax) %>% nrow()

comb_2 %>% distinct(source)

comb_2 %>% filter(source == "biogeotracesAndCultures") %>% summarize(min = min(log10(propOfBiogeotraces)), max = max(log10(propOfBiogeotraces)))


library(viridis)

ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), comb_2 %>% filter(source == "allSites"), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), comb_2 %>% filter(source == "biogeotracesAndCultures") %>%
                filter(propOfBiogeotraces != 0 & !is.na(propOfBiogeotraces)), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_viridis(direction = -1, limits = c(-5, 0), name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuTax) + 
  labs(title = "OTUs in cultures with top 10 highest abundances across BioGeoTraces,\narranged by total number of sequences across BioGeoTraces")
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_abundanctBiogeotracesOTUs_abundantBiogeotracesOTUsAlsoInCultures.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_abundanctBiogeotracesOTUs_abundantBiogeotracesOTUsAlsoInCultures.svg", dpi = 600, height = 10, width = 20)


comb_2 %>% filter(source == "justBiogeotraces") %>% summarize(min = min(log10(propOfBiogeotraces)), max = max(log10(propOfBiogeotraces)))

ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), comb_2 %>% filter(source == "allSites"), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), comb_2 %>% filter(source == "justBiogeotraces") %>%
                filter(propOfBiogeotraces != 0 & !is.na(propOfBiogeotraces)), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_viridis(direction = -1, limits = c(-5, 0), name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuTax) + 
  labs(title = "OTUs with top 10 highest abundances across BioGeoTraces,\narranged by total number of sequences across BioGeoTraces")
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_abundanctBiogeotracesOTUs_abundantBiogeotracesOTUs.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_abundanctBiogeotracesOTUs_abundantBiogeotracesOTUs.svg", dpi = 600, height = 10, width = 20)


##different color palettes

ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), comb_2 %>% filter(source == "allSites"), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), comb_2 %>% filter(source == "biogeotracesAndCultures") %>%
                filter(propOfBiogeotraces != 0 & !is.na(propOfBiogeotraces)), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_gradient(low = "#FEF0D9", high = "#D7301F", limits = c(-5, 0), name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuTax) + 
  labs(title = "OTUs in cultures with top 10 highest abundances across BioGeoTraces,\narranged by total number of sequences across BioGeoTraces")
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_abundanctBiogeotracesOTUs_abundantBiogeotracesOTUsAlsoInCultures_2.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_abundanctBiogeotracesOTUs_abundantBiogeotracesOTUsAlsoInCultures_2.svg", dpi = 600, height = 10, width = 20)




ggplot() +
  geom_point(aes(long, lat, color=propOfBiogeotraces), comb_2 %>% filter(source == "allSites"), size=1.5, color = "black") +
  geom_jitter(aes(long, lat, color=log10(propOfBiogeotraces)), comb_2 %>% filter(source == "justBiogeotraces") %>%
                filter(propOfBiogeotraces != 0 & !is.na(propOfBiogeotraces)), size=2) +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  scale_color_gradient(low = "#FEF0D9", high = "#D7301F", limits = c(-5, 0), name = expression(paste("log"[10],"(proportion of non-Pro., non-Syn. sequences\nin BioGeoTraces sample)"))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + facet_wrap(~otuTax) + 
  labs(title = "OTUs with top 10 highest abundances across BioGeoTraces,\narranged by total number of sequences across BioGeoTraces")
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("seanSeqsNoProSyn_abundanctBiogeotracesOTUs_abundantBiogeotracesOTUs_2.png", dpi = 600, height = 10, width = 20)
ggsave("seanSeqsNoProSyn_abundanctBiogeotracesOTUs_abundantBiogeotracesOTUs_2.svg", dpi = 600, height = 10, width = 20)

