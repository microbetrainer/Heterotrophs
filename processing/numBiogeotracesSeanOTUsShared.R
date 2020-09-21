library(tidyverse)

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

#loads file with the sequences in sean culture samples, after 
#excluding jw3 and 1223 samples
#this includes pro and syn sequences
seanWithProSyn <- read.table("sequencesForTree_noContamination_noMock_withProSyn_no1223.txt")
seanWithProSyn %>% head()
seanWithProSyn %>% nrow()

biomG %>% head()

#these are the sean sequences
biomG %>% filter(!str_detect(sample, "SRR")) %>% head()

#makes the sean sequence ID format match the format of the sequences 
#in seanWithProSyn
str_extract("X240sean", "[0-9]{1,10000}")
biomG <- biomG %>% mutate(sample = ifelse(!str_detect(sample, "SRR"), str_c("contig_", str_extract(sample, "[0-9]{1,10000}")), sample))

#the new format of the sean sequences
biomG %>% filter(!str_detect(sample, "SRR")) %>% head()

biomG %>% head()
seanWithProSyn %>% head()

biomG %>% filter(str_detect(sample, "contig_")) %>% nrow()
biomG %>% semi_join(seanWithProSyn, by = c("sample" = "V1")) %>% nrow()

#makes a dataframe for the sean sequences in seanWithProSyn 
biom_sean <- biomG %>% semi_join(seanWithProSyn, by = c("sample" = "V1"))
biomG %>% filter(str_detect(sample, "contig_")) %>% nrow()
biom_sean %>% nrow()

#makes a dataframe for the biogeotraces sequences
biom_biogeo <- biomG %>% filter(str_detect(sample, "SRR"))
biom_biogeo %>% nrow()

biom_sean
biom_biogeo

###with pro and syn

biom_sean %>% distinct(sample) %>% nrow()
biom_biogeo %>% distinct(sample) %>% nrow()

#without excluding pro and syn, there are 98 sean OTUs and 28008 biogeotraces OTUs
biom_sean %>% distinct(OTU_ID) %>% nrow()
biom_biogeo %>% distinct(OTU_ID) %>% nrow()

#without excluding pro and syn, 84 OTUs are in common between sean and biogeotraces samples
biom_sean %>% semi_join(biom_biogeo, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

#without excluding pro and syn, 14 OTUs are unique to sean samples
biom_sean %>% anti_join(biom_biogeo, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

#without excluding pro and syn, 27924 OTUs are unique to biogeotraces samples
biom_biogeo %>% anti_join(biom_sean, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()


###without pro and syn

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

biomG %>% distinct(OTU_ID) %>% nrow()
tax %>% distinct(`Feature ID`) %>% nrow()

biomG %>% head()
tax %>% head()

#adds taxonomy of each OTU to numOTUs
nrow(biomG)
biomG <- biomG %>% left_join(tax, by = c("OTU_ID" = "Feature ID"))
nrow(biomG)

#biom_sean %>% semi_join(biom_biogeo, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% 
  #left_join(tax, by = c("OTU_ID" = "Feature ID")) %>% distinct(Taxon)  %>% 
  #write.table("~/Desktop/testWithProSyn.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

biomG %>% head()

biomG %>% filter(is.na(Taxon))

biomG %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)
biomG %>% filter(str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc")) %>% distinct(Taxon)

#excludes pro and syn OTUs
biomG <- biomG %>% filter(!str_detect(Taxon, "Synechococcaceae|synechococcaceae|Proc|proc"))

biomG %>% filter(str_detect(Taxon, "Synec|synec|Proc|proc")) %>% distinct(Taxon)

#makes a dataframe for sean sequences
head(biomG)
head(seanWithProSyn)
biom_sean <- biomG %>% semi_join(seanWithProSyn, by = c("sample" = "V1"))

#makes a dataframe for biogeotraces sequences
biom_biogeo <- biomG %>% filter(str_detect(sample, "SRR"))

biom_sean %>% distinct(sample) %>% nrow()
biom_biogeo %>% distinct(sample) %>% nrow()

#after excluding pro and syn sequences, there are 97 OTUs in sean samples
biom_sean %>% distinct(OTU_ID) %>% nrow()

#after excluding pro and syn sequences, there are 27966 OTUs in biogeotraces samples
biom_biogeo %>% distinct(OTU_ID) %>% nrow()

#after excluding pro and syn sequences, there are 83 OTUs in common between sean and 
#biogeotraces samples
biom_sean %>% semi_join(biom_biogeo, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

#after excluding pro and syn sequences, there are 14 OTUs unique to sean samples
biom_sean %>% anti_join(biom_biogeo, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

#after excluding pro and syn sequences, there are 27883 OTUs unique to biogeotraces samples
biom_biogeo %>% anti_join(biom_sean, by = c("OTU_ID")) %>% distinct(OTU_ID) %>% nrow()

##this is all checked 

#biom_sean %>% semi_join(biom_biogeo, by = c("OTU_ID")) %>% distinct(Taxon) %>% 
  #write.table("~/Desktop/test.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
biom_sean %>% semi_join(biom_biogeo, by = c("OTU_ID")) %>% distinct(Taxon) %>% 
  mutate(Taxon = str_replace_all(Taxon, "[a-z]__", "")) %>% 
  mutate(Taxon = str_replace_all(Taxon, ";", "")) %>% 
  mutate(Taxon = str_replace_all(Taxon, " *$", "")) %>% distinct(Taxon)

