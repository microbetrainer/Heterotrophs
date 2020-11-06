library(tidyverse)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads taxonomy of Sean features
#for this, I excluded sample 1223
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223.tsv")

#gets rid of first row that just has notes
tax <- tax[-1,]

#this is the correct number of rows
tax %>% nrow()

#extracts all text after k__ because all of sequences have at least this
tax <- tax %>% mutate(taxa = str_extract(Taxon, "k__.*"))

tax %>% filter(is.na(taxa))

#gets rid of letter and dash combinations that tell you the level of taxonomy

tax$taxa <- str_replace(tax$taxa, "k__", "")

tax$taxa <- str_replace(tax$taxa, "; p__", " ")

tax$taxa <- str_replace(tax$taxa, "; c__", " ")

tax$taxa <- str_replace(tax$taxa, "; o__", " ")

tax$taxa <- str_replace(tax$taxa, "; f__", " ")

tax$taxa <- str_replace(tax$taxa, "; g__", " ")

tax$taxa <- str_replace(tax$taxa, "; s__", " ")

tax %>% filter(is.na(taxa))


#loads representative sequences so that I know which sequence IDs correspond
#for this, I excluded sample 1223
repSeqs <- read.table("representative_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.fasta", fill = TRUE)

#gets just the rows that correspond to the IDs
repSeqs <- repSeqs %>% filter(str_detect(V1, ">"))

#gets rid of ">" in IDs
repSeqs <- repSeqs %>% mutate(V1 = str_replace(V1, ">", ""))

nrow(repSeqs)

nrow(tax)

#adds alternative sequence ID to tax
tax <- tax %>% left_join(repSeqs, by = c("Feature ID" = "V1"))

nrow(tax)

tax %>% filter(is.na(V2))



#loads the sequences that I want to be included in the tree, excluding sequences 
#just in JW3 and just in 1223
tree <- read.table("sequencesForTree_noContamination_noMock_noProSyn_no1223.txt")

nrow(tree)

#adds taxonomy to sequences in tree
tree <- tree %>% left_join(tax, by = c("V1" = "V2"))

nrow(tree)

tree %>% filter(is.na(taxa))

#gets just the feature and cleaned up taxonomy variables from tax
#tax <- tax %>% select(1,4)

#makes a variable for the shortened taxonomy
tree <- tree %>% mutate(shortTaxa = str_extract(taxa, "[A-Z][a-z\\s0-9I]*$"))

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "OD1"), "OD1", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree %>% filter(str_detect(taxa, "Deltaproteobacteria Sva0853") & str_detect(shortTaxa, "Sva0853"))

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Deltaproteobacteria Sva0853") & str_detect(shortTaxa, "Sva0853"), "Deltaproteobacteria Sva0853", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Acidimicrobiales OCS155"), "Acidimicrobiales OCS155", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Alteromonadales HTCC2188 HTCC"), "Alteromonadales HTCC2188 HTCC", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Deltaproteobacteria PB19"), "Deltaproteobacteria PB19", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Thermoplasmata E2 Marine group II"), "Thermoplasmata E2 Marine group II", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "KSA1"), "Balneolaceae KSA1", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "SBR1093 A712011"), "SBR1093 A712011", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "WGS"), "SAR406 AB16 Arctic96B-7 A714017 SargSea-WGS", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Deltaproteobacteria Sva0853 SAR324"), "Deltaproteobacteria Sva0853 SAR324", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Alteromonadales OM60"), "Alteromonadales OM60", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Gammaproteobacteria HTCC2188 HTCC2089"), "Gammaproteobacteria HTCC2188 HTCC2089", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "SGSH944"), "SAR406 AB16 Arctic96B-7 A714017 SGSH944", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Alteromonadales J115"), "Alteromonadales J115", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Acidimicrobiales C111"), "Acidimicrobiales C111", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Verrucomicrobiaceae MSBL3"), "Verrucomicrobiaceae MSBL3", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Myxococcales OM27"), "Myxococcales OM27", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Deltaproteobacteria GMD14H09"), "Deltaproteobacteria GMD14H09", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Alteromonadaceae ZD0117"), "Alteromonadaceae ZD0117", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Bacteria Proteobacteria Alphaproteobacteria Rickettsiales AEGEAN_112"), "Rickettsiales AEGEAN_112", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Bacteria SAR406 AB16 ZA3648c AEGEAN_185"), "SAR406 AB16 ZA3648c AEGEAN_185", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Bacteria SAR406 AB16 Arctic96B-7 A714017 Arctic95A-2"), "SAR406 AB16 Arctic96B-7 A714017 Arctic95A-2", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "SAR406 AB16 Arctic96B-7 Sc-NB04"), "SAR406 AB16 Arctic96B-7 Sc-NB04", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

tree %>% filter(str_detect(taxa, "Bacteria Proteobacteria Gammaproteobacteria HTCC2188\\s*$"))

tree <- tree %>% mutate(shortTaxa = ifelse(str_detect(taxa, "Bacteria Proteobacteria Gammaproteobacteria HTCC2188\\s*$"), "Gammaproteobacteria HTCC2188", shortTaxa))

#tree %>% distinct(taxa, shortTaxa) %>% View()

d <- tree %>% distinct(taxa, shortTaxa)

nrow(tree)



#tree <- tree %>% select(V1, shortTaxa)


str_extract(tree$shortTaxa, "\\s$")

#gets rid of space at end of shortTaxa
tree <- tree %>% mutate(shortTaxa = str_replace(shortTaxa, "\\s*$", ""))


write_tsv(tree, "taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv")


#combines the two sequence IDs into one variable
tree <- tree %>% mutate(V1 = str_c(`Feature ID`, V1, sep = " "))

#tree %>% filter(V1 == "contig_7")

#tree <- tree %>% mutate(shortTaxa = ifelse(V1 == "contig_7", "Alteromonas macleodii MIT1002", shortTaxa))

#adds quotes to the feature IDs and short taxonomies because 
#this is needed to replace tree node names
tree <- tree %>% mutate(V1 = str_c("'", V1, "'"))
tree <- tree %>% mutate(shortTaxa = str_c("'", shortTaxa, "'"))

#makes variable for combined feature number and short taxonomy
tree <- tree %>% mutate(replacement = str_c(V1, shortTaxa, sep = ":"))

#gets just the combined feature and short taxonomy variable
tree <- tree %>% select(replacement)

nrow(tree)

treeList <- tree %>% as.list()

treeList[[1]]

length(treeList[[1]])

paste(treeList[[1]], collapse = ",")

write.table(paste(treeList[[1]], collapse = ","), "replaceTreeNodeNames_noContamination_noMock_noProSyn_withoutSample1223.txt")





