library(tidyverse)
library(phyloseq)
library(vegan)

setwd("/Users/skearney/Documents/prochlorococcus/experiments/Culture Collection/Metadata")

#loads corrected counts of each taxa in each sample 
#this excludes sequences that were determined to be contamination from nearby wells
seq <- read.table("seqtab_corrected.txt")

numCols <- ncol(seq)

cols <- str_c("col", 1:numCols)

#gives a number to each sequence (column)
colnames(seq) <- cols

#gets list of samples
samples <- rownames(seq)

#makes a variable for the sample in seq instead of the samples just being
#the rownames
seq <- seq %>% mutate(sample = samples)

#gathers seq into long form
seqG <- seq %>% gather(1:3170, key = "sequence", value = "abundance")

#gets the total number of sequences in each sample
summ <- seqG %>% group_by(sample) %>% summarize(totalSeqs = sum(abundance))

#adds number of sequences in sample to each sequence
seqG <- seqG %>% left_join(summ, by = c("sample"))

#gets rid of the rows for the samples that are unrelated (Sean told me which ones to exclude)
seqG <- seqG %>% filter(!str_detect(sample, "Con|Epi|Pcn|Lin"))

#gets rid of Mock samples
seqG <- seqG %>% filter(!str_detect(sample, "Mock"))

#gets rid of the rows for which a sequence is not present in a sample
#this excludes the sequences that are just present in the unrelated samples (Con|Epi|Pcn|Lin)
seqG <- seqG %>% filter(abundance > 0)

#gets the proportion that each sequence is in each sample
seqG <- seqG %>% mutate(propSample = abundance/totalSeqs)

#excludes the sequences in samples that have a low relative abundances in the sample
seqG <- seqG %>% filter(propSample >= 0.002)

#gets rid of unnecessary variables
seqG <- seqG %>% select(-c(abundance, totalSeqs))

#makes sequence names in seqG match sequence names in repSeqs
seqG$sequence <- str_replace(seqG$sequence, "col", "contig_")

#loads the representative sequences fasta file so that I can get the 
#contig ID assigned to each feature ID
#this is with excluding Pro and Syn
repSeqs <- read.table("representative_sameLengthSeqs_noContamination_noMock_noProSyn_withoutSample1223.fasta", fill = TRUE)

#gets just the rows in repSeqs that correspond to IDs
repSeqs <- repSeqs %>% filter(str_detect(V1, ">"))

#adds contig IDs to feature IDs in seqG
seqG <- seqG %>% left_join(repSeqs, by = c("sequence" = "V2"))

#gets rid of ">" in feature IDs
seqG <- seqG %>% mutate(V1 = str_replace(V1, ">", ""))

#makes feature ID, contig ID combined variable that matches 
#the IDs used in the tree file
seqG <- seqG %>% mutate(featureContig = str_c(V1, sequence, sep = ""))

#loads cleaned up taxonomy of the same length features
#corresponds to representative_sameLengthSeqs.fasta
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv")

#adds taxonomy to seqG
seqG <- seqG %>% left_join(tax, by = c("V1" = "Feature ID"))

#gets rid of pro and syn sequences
seqG <- seqG %>% filter(!is.na(shortTaxa))

#JW3 only has Pro sequences so it gets excluded which is good because it only has 
#one sequence
seqG <- seqG %>% filter(sample != "JW3")

#excludes sample 1223
seqG <- seqG %>% filter(sample != "1223")

#according to Sean:
#pro: JW7, B7, DV, C12B, JW2, JW4, 1213
#syn: 8102, RS9916

#gets samples that are either Pro or Syn, which is all of them
seqG <- seqG %>% filter(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916|SYN|WH|9220|S9503|8102|RS9916|0604|9202|9215|9314|9321|9322|9401|MED1|MED4|NATL1A|NATL2A|9313|9211|9303|1214|1201|ASN9601|9123|1300|1205|1304|SS35|SS51|SS52|SS2|B5|C4|C8|9302|9201.1|9201.2|LG|C9B|0701|0702|0703|1227|0601|GP2|9301|9312|9515|SB|SS120|9107|0602|0603|1307|1341|1223|JW7|B7|DV|C12B|JW2|JW4|1213"))

#makes a variable for whether sample is Pro or Syn
seqG <- seqG %>% mutate(culture = ifelse(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916"), "Synechococcus", "Prochlorococcus"))

#gets rid of variables so I can spread data
seqG <- seqG %>% select(sample, culture, featureContig, propSample)

#spreads seqG so that the sequences that are not in a sample will have 0 abundance for that sample
seqSpread <- seqG %>% spread(key = featureContig, value = propSample, fill = 0)

#gathers seqG again, now with rows for sequences that are not present in a sample
seqG <- seqSpread %>% gather(3:237, key = "V1", value = "propSample")

rownames(seqSpread) <- seqSpread$sample

#makes column names not start with numbers
colnames(seqSpread) <- str_c("col", colnames(seqSpread))

#loads cleaned up environmental variables for each culture
#this came from cleanUpEnvVariablesOfCulturesForPCA.R
env <- read_csv("051820Metadata_SK.csv")

##from comparingSampleIDs.R
#1213 is 1313
#WH7803 is WH7803
#SYN1320 is SYN1220
#9201.2 is 9311

#AS9601 is correct ID rather than ASN9601
#SYN9503 is correct ID rather than S9503

seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "1213", "1313", colsample))
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "9201.2", "9311", colsample))
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "ASN9601", "AS9601", colsample))
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "9201.1", "9201", colsample))
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "S9503", "SYN9503", colsample))
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "SYN1320", "SYN1220", colsample))

#seqG$sample <- ifelse(seqG$sample == "8102", "SYNCLADEXVI", seqG$sample)
seqSpread <- seqSpread %>% mutate(colsample = ifelse(colsample == "8102", "SYNCLADEXVI", colsample))


#the entries that I don't know the month, day, or both still have a dayMonth value

#makes a variable for the month that cultures were isolated
# env <- env %>% mutate(month = str_extract(`DATE ISOLATED_Edited`, "(?<=-)[0-9]*(?=-)"))
# 
# #all I know is that WH6501 is from 1965
# #the culture website and excel sheet say MED4 and MED1 were isolated Jan. 1989
# env %>% filter(month == "01")
# 
# #makes the month NA for WH6501
# env <- env %>% mutate(month = ifelse(CULTURE == "WH6501", NA, month))
# 
# #makes a variable with the month name
# env <- env %>% mutate(monthNew = ifelse(month == "01", "Jan.", month))
# env <- env %>% mutate(monthNew = ifelse(month == "02", "Feb.", monthNew))
# env <- env %>% mutate(monthNew = ifelse(month == "03", "Mar.", monthNew))
# env <- env %>% mutate(monthNew = ifelse(month == "04", "Apr.", monthNew))
# env <- env %>% mutate(monthNew = ifelse(month == "05", "May.", monthNew))
# env <- env %>% mutate(monthNew = ifelse(month == "06", "Jun.", monthNew))
# env <- env %>% mutate(monthNew = ifelse(month == "07", "Jul.", monthNew))
# env <- env %>% mutate(monthNew = ifelse(month == "08", "Aug.", monthNew))
# env <- env %>% mutate(monthNew = ifelse(month == "09", "Sep.", monthNew))
# env <- env %>% mutate(monthNew = ifelse(month == "10", "Oct.", monthNew))
# env <- env %>% mutate(monthNew = ifelse(month == "11", "Nov.", monthNew))
# env <- env %>% mutate(monthNew = ifelse(month == "12", "Dec.", monthNew))

#makes year of isolation variable
#env <- env %>% mutate(year = str_extract(`DATE ISOLATED_Edited`, "[0-9]*"))
#env <- env %>% mutate(year = as.numeric(year)) 

#makes depth range variable
env <- env %>% mutate(depthRange = ifelse(Depth <= 50, "0-50m", NA))
env <- env %>% mutate(depthRange = ifelse(Depth > 50 & Depth <= 100, "50-100m", depthRange))
env <- env %>% mutate(depthRange = ifelse(Depth > 100 & Depth <= 150, "100-150m", depthRange))
env <- env %>% mutate(depthRange = ifelse(Depth > 150 & Depth <= 200, "150-200m", depthRange))

#makes a variable that has both the month and year of isolation
#env <- env %>% mutate(monthNewYear = str_c(as.character(env$monthNew), as.character(env$year), sep = " "))

#makes a variable for the age of the culture at the time they were sequenced
#env <- env %>% mutate(yearsInAge = interval(ymd(env$`DATE ISOLATED_Edited`), ymd("2019-06-01 UTC"))/years(1))

#fixes "Possibly EqPac/IRONEX" in cruise variable
#env <- env %>% mutate(CRUISE = ifelse(str_detect(CRUISE, "EqPac/IRONEX\\?"), "Possibly EqPac/IRONEX", CRUISE))


###phylosift

#https://joey711.github.io/phyloseq/import-data.html#phyloseq-ize_data_already_in_r

#phyloseq - Takes as argument an otu_table and any unordered list of valid phyloseq components: 
#sample_data, tax_table, phylo, or XStringSet. The tip labels of a phylo-object (tree) must match 
#the OTU names of the otu_table, and similarly, the sequence names of an XStringSet object must match 
#the OTU names of the otu_table

phylo_seqSpread <- seqSpread[c(1:237)]

#makes the rownames in seqSpread the sample
rownames(phylo_seqSpread) <- phylo_seqSpread$colsample

#gets rid of the sample variable in seqSpread
phylo_seqSpread <- phylo_seqSpread %>% select(-c(colsample, colculture))

colnames(phylo_seqSpread) <- str_replace(colnames(phylo_seqSpread), "col", "")

#makes phyloseq object out of otu table
otu <- otu_table(phylo_seqSpread, taxa_are_rows = FALSE)

#loads cleaned up taxonomy of the same length features, no pro, no syn, no contamination, no mock samples, no 1223 or JW3
#this is from makeTreeNodes_noContamination_noMock_noProSyn_without1223.R
tax <- read_tsv("taxonomySameLength_noContamination_noMock_noProSyn_withoutSample1223_cleanedUp.tsv")

tax <- tax %>% mutate(featureContig = str_c(`Feature ID`, V1))

#there are no Pro or Syn features in tax
#makes dataframe with the taxonomy for each feature
taxPhyloseq <- tax %>% distinct(featureContig, shortTaxa)

#makes rownames of taxPhyloseq the feature ID
rownames(taxPhyloseq) <- taxPhyloseq$featureContig

#gets rid of feature ID variable
taxPhyloseq <- taxPhyloseq %>% select(-featureContig)

#makes phyloseq object out of taxPhyloseq
TAX <- tax_table(as.matrix(taxPhyloseq))

#loads rooted tree
#for this, I excluded sample 1223
tree <- read_tree("tree.nwk")


env <- env %>% as.data.frame()
row.names(env) <- env$OldName
#env <- env %>% select(-c(CULTURE, `ALTERNATIVE NAME`, DEPTH, Notes, latTotal, lonTotal, taxa, `DATE ISOLATED_Edited`, month))
colnames(env)

#rearranges columns
#env <- env %>% select(cultureType:`METHOD (for isolation)`, oceanType, monthNew, depthRange, year, monthNewYear, yearsInAge)

x <- env$`Date Isolated`
x.new <- sapply(x,function(h) format(as.Date(h),"%Y/%m/%d"))
x.num <- sapply(x.new,function(x) as.numeric(as.Date(x)))

date.measured <- as.Date("2019/05/01")
dm.num <- as.numeric(x.num)



# x[which(x=="Jan.")] <- 0
# x[which(x=="Feb.")] <- 1
# x[which(x=="Mar.")] <- 2
# x[which(x=="Apr.")] <- 3
# x[which(x=="May.")] <- 4
# x[which(x=="Jun.")] <- 5
# x[which(x=="Jul.")] <- 6
# x[which(x=="Aug.")] <- 7
# x[which(x=="Sep.")] <- 8
# x[which(x=="Oct.")] <- 9
# x[which(x=="Nov.")] <- 10
# x[which(x=="Dec.")] <- 11
#
x <- as.numeric(x)
y <- as.numeric(env$year)
yx <- y + x/12
z <- 2019+5/12
age <- z-yx
#
# env$age <- age

# add in Longhurst Data
#lhrst <- read.csv(file = "longhurst.csv")
#crs <- read.csv(file = "cruise.csv")
#env$CRUISE <- sapply(crs[,1],as.character)
#env$LonghurstCode <- sapply(lhrst[,1],as.character)
#env$LonghurstProvince <- sapply(lhrst[,2],as.character)
#env$Biome <- sapply(lhrst[,3],as.character)


#str(env)
#row.names(env) <- env$CULTURE_originalName
sampleData <- sample_data(env)

#makes combined phyloseq object out of otu, TAX, tree, and sample data
physeq <- phyloseq(otu, TAX, tree, sampleData)

set.seed(7)

#there is a warning
unifracRes <- UniFrac(physeq, weighted = FALSE)

##sample IDs in unifracRes and env need to be in the same order for permanova
order <- labels(unifracRes) %>% as.list()
order
env$OldName <- factor(env$OldName, levels = order)

levels(env$OldName) %>% length()
labels(unifracRes) %>% length()

levels(env$OldName)
labels(unifracRes)

unique(levels(env$OldName) == labels(unifracRes))

env <- env %>% arrange(OldName)

env$OldName
labels(unifracRes)
unique(env$OldName == labels(unifracRes))

#example of permanova when there are no NAs in metadata
#set.seed(7)

#adonis2(unifracRes~Genus, data=env, permutations = 999, method="unifrac", strata="Genus")

#example of permanova when there are NAs in metadata
#set.seed(7)
#adonis2(UniFrac(subset_samples(physeq, (!is.na(Ecotype) & !is.na(Genus))), weighted = FALSE)~Ecotype,
#data=env %>% filter(!is.na(Ecotype) & !is.na(Genus)), permutations = 999, method="unifrac")

#example with Jaccard distance
phyloseq::distance(physeq,method="jaccard") -> jac
adonis2(jac~Genus, data=env, permutations = 999, method="jaccard")


#permanova each metaadata variable
set.seed(1234)
#cultureType
M01 <- UniFrac(subset_samples(physeq, !is.na(Genus)), weighted = FALSE)
S01 <- sample_data(physeq) %>% filter(!is.na(Genus))
PA01 <- adonis2(M01 ~ Genus, data=S01, permutations = 999, method="unifrac")

#Clade
M02 <- UniFrac(subset_samples(physeq, !is.na(Clade)), weighted=FALSE)
S02 <- sample_data(physeq) %>% filter(!is.na(Clade))
PA02 <- adonis2(M02 ~ Clade, data=S02, permutations = 999, method="unifrac")

#depthRange
M03 <- UniFrac(subset_samples(physeq, !is.na(depthRange)), weighted=FALSE)
S03 <- sample_data(physeq) %>% filter(!is.na(depthRange))
PA03 <- adonis2(M03 ~ depthRange, data=S03, permutations = 999, method="unifrac")

#Location
M04 <- UniFrac(subset_samples(physeq, !is.na(Location)), weighted =FALSE)
S04 <- sample_data(physeq) %>% filter(!is.na(Location))
A04 <- adonis2(M04 ~ Location, data=sample_data(physeq) %>% filter(!is.na(Location)), permutations = 999, method="unifrac")

#Isolation Method
M05 <- UniFrac(subset_samples(physeq, !is.na(Isolation.Method)), weighted = FALSE)
S05 <- sample_data(physeq) %>% filter(!is.na(Isolation.Method))
PA05 <- adonis2(M05 ~ Isolation.Method, 
               data=sample_data(physeq) %>% filter(!is.na(Isolation.Method)), permutations = 999, method="unifrac")

#Ecotype
M06 <- UniFrac(subset_samples(physeq, !is.na(Ecotype)), weighted = FALSE)
S06 <- sample_data(physeq) %>% filter(!is.na(Ecotype))
PA06 <-  adonis2(M06 ~ Ecotype, 
                data=sample_data(physeq) %>% filter(!is.na(Ecotype)), permutations = 999, method="unifrac")

#Cruise Name
M07 <- UniFrac(subset_samples(physeq, !is.na(Cruise.Name)), weighted = FALSE)
S07 <- sample_data(physeq) %>% filter(!is.na(Cruise.Name))

PA07 <- adonis2(M07 ~ Cruise.Name, 
               data=sample_data(physeq) %>% filter(!is.na(Cruise.Name)), permutations = 999, method="unifrac")

#monthNewYear
M08 <- UniFrac(subset_samples(physeq, !is.na(monthNewYear)), weighted = FALSE)
S08 <- sample_data(physeq) %>% filter(!is.na(monthNewYear))
PA08 <- adonis2(M08 ~ monthNewYear, data=sample_data(physeq) %>% filter(!is.na(monthNewYear)), permutations = 999, method="unifrac")

PA.pvalues <- c(PA01[[5]][1],PA02[[5]][1],PA03[[5]][1],PA04[[5]][1],PA05[[5]][1],PA06[[5]][1],PA07[[5]][1],PA08[[5]][1]) 

#GKTau for environmental variable associations
library(GoodmanKruskal)
mdata <- sample_data(physeq)
mdata.clean <- mdata[,c(-2,-7,-9,-10,-11,-12)]
names(mdata.clean) <- c("Genus","Ecotype","Clade","Location","Cruise","Isolation.Method","Depth","Date")
psenv.l <- length(mdata.clean)
gkmat <- matrix(NA,nrow = psenv.l,ncol=psenv.l)
mdata.names <- names(mdata.clean)
colnames(gkmat) <- rownames(gkmat) <- mdata.names
for (i in 1:psenv.l){
  for (j in 1:psenv.l){
    if (i <= j){
      tmp <- GKtau(mdata.clean[[i]],mdata.clean[[j]])
      a <- tmp$tauxy
      b <- tmp$tauyx
      gkmat[i,j] <- a
      gkmat[j,i] <- b
    }
  }
}

corrplot::corrplot(gkmat,order="hclust",cl.lim = c(0,1),tl.col = "black",method="color")


#testing some NDMS
library(patchwork)
library(ggrepel)
set.seed(1234)
culture.ord <- phyloseq::ordinate(physeq,"NMDS","unifrac")

p1 = plot_ordination(physeq, culture.ord, color="Genus") + geom_point(size=3) + theme_linedraw() #+ theme(legend.position= "none")
p2 = plot_ordination(physeq, culture.ord, color="depthRange") + geom_point(size=3) + theme_linedraw() #+ theme(legend.position= "none")
p3 = plot_ordination(physeq, culture.ord, color="Cruise.Name") + geom_point(size=3) + theme_linedraw() #+ theme(legend.position = "none")
p4 = plot_ordination(physeq, culture.ord, color="Isolation.Method") + geom_point(size=3) + theme_linedraw() #+ theme(legend.position = "none")

(p1 + p2) / (p3 + p4)

p5 = plot_ordination(physeq, culture.ord, color="Ecotype",shape="Genus") + geom_point(size=3) + theme_linedraw() #+ theme(legend.position = "none")
p6 = plot_ordination(physeq, culture.ord, color="Clade",shape="Genus") + geom_point(size=3) + theme_linedraw() #+ theme(legend.position = "none")
p7 = plot_ordination(physeq, culture.ord, color="Location") + geom_point(size=3) + theme_linedraw() #+ theme(legend.position = "none")
p8 = plot_ordination(physeq, culture.ord, color="yearsInAge") + geom_point(size=3) + theme_linedraw() + scale_color_gradient(low="#fff7bc",high = "#d95f0e") #+ theme(legend.position = "none")

# plot just the NMDS with sample labels 
p9 = plot_ordination(physeq, culture.ord) + geom_point(size=3) + theme_linedraw() + geom_label_repel(mapping = aes(label = Culture.Name), size = 6, segment.size=0.1) #+ scale_color_gradient(low="#fff7bc",high = "#d95f0e") 

(p5 + p6) / (p7 + p8)


# testing of all pairwise interaction effects
# Genus/Ecotype NULL
A01 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Genus) & !is.na(Ecotype)), weighted = FALSE) ~ Genus*Ecotype, 
               data=sample_data(physeq) %>% filter(!is.na(Genus) & !is.na(Ecotype)), permutations = 999, method="unifrac",by="margin")
# Genus/Clade NULL
A02 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Genus) & !is.na(Clade)), weighted = FALSE) ~ Genus*Clade, 
               data=sample_data(physeq) %>% filter(!is.na(Genus) & !is.na(Clade)), permutations = 999, method="unifrac",by="margin")
# Genus/Depth
A03 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Genus) & !is.na(depthRange)), weighted = FALSE) ~ Genus*depthRange, 
               data=sample_data(physeq) %>% filter(!is.na(Genus) & !is.na(depthRange)), permutations = 999, method="unifrac",by="margin")
# Genus/Cruise
A04 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Genus) & !is.na(Cruise.Name)), weighted = FALSE) ~ Genus*Cruise.Name, 
               data=sample_data(physeq) %>% filter(!is.na(Genus) & !is.na(Cruise.Name)), permutations = 999, method="unifrac",by="margin")

# Genus/Isolation
A05 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Genus) & !is.na(Isolation.Method)), weighted = FALSE) ~ Genus*Isolation.Method, 
               data=sample_data(physeq) %>% filter(!is.na(Genus) & !is.na(Isolation.Method)), permutations = 999, method="unifrac",by="margin")
# Genus/Location
A06 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Genus) & !is.na(Location)), weighted = FALSE) ~ Genus*Location, 
               data=sample_data(physeq) %>% filter(!is.na(Genus) & !is.na(Location)), permutations = 999, method="unifrac",by="margin")

# Ecotype/Clade NULL
A07 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Clade) & !is.na(Ecotype)), weighted = FALSE) ~ Clade*Ecotype, 
               data=sample_data(physeq) %>% filter(!is.na(Clade) & !is.na(Ecotype)), permutations = 999, method="unifrac",by="margin")

# Ecotype/Depth
A08 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Ecotype) & !is.na(depthRange)), weighted = FALSE) ~ Ecotype*depthRange, 
               data=sample_data(physeq) %>% filter(!is.na(Ecotype) & !is.na(depthRange)), permutations = 999, method="unifrac",by="margin")

# Ecotype/Cruise
A09 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Ecotype) & !is.na(Cruise.Name)), weighted = FALSE) ~ Ecotype*Cruise.Name, 
               data=sample_data(physeq) %>% filter(!is.na(Ecotype) & !is.na(Cruise.Name)), permutations = 999, method="unifrac",by="margin")

# Ecotype/Isolation
A10 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Ecotype) & !is.na(Isolation.Method)), weighted = FALSE) ~ Ecotype*Isolation.Method, 
                      data=sample_data(physeq) %>% filter(!is.na(Ecotype) & !is.na(Isolation.Method)), permutations = 999, method="unifrac",by="margin")

# Ecotype/Location
A11 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Ecotype) & !is.na(Location)), weighted = FALSE) ~ Ecotype*Location, 
               data=sample_data(physeq) %>% filter(!is.na(Genus) & !is.na(Location)), permutations = 999, method="unifrac",by="margin")

# Clade/Depth
A12 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Clade) & !is.na(depthRange)), weighted = FALSE) ~ Clade*depthRange, 
               data=sample_data(physeq) %>% filter(!is.na(Clade) & !is.na(depthRange)), permutations = 999, method="unifrac",by="margin")

# Clade/Cruise NULL
A13 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Clade) & !is.na(Cruise.Name)), weighted = FALSE) ~ Clade*Cruise.Name, 
        data=sample_data(physeq) %>% filter(!is.na(Clade) & !is.na(Cruise.Name)), permutations = 999, method="unifrac",by="margin")

# Clade/Isolation
A14 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Clade) & !is.na(Isolation.Method)), weighted = FALSE) ~ Clade*Isolation.Method, 
               data=sample_data(physeq) %>% filter(!is.na(Clade) & !is.na(Isolation.Method)), permutations = 999, method="unifrac",by="margin")

# Clade/Location
A15 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Clade) & !is.na(Location)), weighted = FALSE) ~ Clade*Location, 
        data=sample_data(physeq) %>% filter(!is.na(Clade) & !is.na(Location)), permutations = 999, method="unifrac",by="margin")

# Depth/Cruise
A16 <- adonis2(UniFrac(subset_samples(physeq, !is.na(depthRange) & !is.na(Cruise.Name)), weighted = FALSE) ~ depthRange*Cruise.Name, 
        data=sample_data(physeq) %>% filter(!is.na(depthRange) & !is.na(Cruise.Name)), permutations = 999, method="unifrac",by="margin")

# Depth/Isolation
A17 <- adonis2(UniFrac(subset_samples(physeq, !is.na(depthRange) & !is.na(Isolation.Method)), weighted = FALSE) ~ depthRange*Isolation.Method, 
        data=sample_data(physeq) %>% filter(!is.na(depthRange) & !is.na(Isolation.Method)), permutations = 999, method="unifrac",by="margin")

# Depth/Location
A18 <- adonis2(UniFrac(subset_samples(physeq, !is.na(depthRange) & !is.na(Location)), weighted = FALSE) ~ depthRange*Location, 
               data=sample_data(physeq) %>% filter(!is.na(depthRange) & !is.na(Location)), permutations = 999, method="unifrac",by="margin")

# Cruise/Isolation
A19 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Cruise.Name) & !is.na(Isolation.Method)), weighted = FALSE) ~ Cruise.Name*Isolation.Method, 
               data=sample_data(physeq) %>% filter(!is.na(Cruise.Name) & !is.na(Isolation.Method)), permutations = 999, method="unifrac",by="margin")

# Cruise/Location
A20 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Cruise.Name) & !is.na(Location)), weighted = FALSE) ~ Cruise.Name*Location, 
               data=sample_data(physeq) %>% filter(!is.na(Cruise.Name) & !is.na(Location)), permutations = 999, method="unifrac",by="margin")

# Isolation/Location
A21 <- adonis2(UniFrac(subset_samples(physeq, !is.na(Location) & !is.na(Isolation.Method)), weighted = FALSE) ~ Location*Isolation.Method, 
               data=sample_data(physeq) %>% filter(!is.na(Location) & !is.na(Isolation.Method)), permutations = 999, method="unifrac",by="margin")


### Additional Post-hoc tests requested by the editor
## Beta dispersion 
# for cruise names
Depth.PostHoc <- pairwise.adonis(M03,S03$depthRange,p.adjust='holm')
Cruise.PostHoc <- pairwise.adonis(M07,S07$Cruise.Name,p.adjust='holm')


# Additional Age Regression Plot

df1 <- as.data.frame(cbind(sample_data(physeq)$yearsInAge,rowSums(otu_table(physeq) > 0)))
colnames(df1) <- c("Age","Richness")

# regression ggplot

ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + geom_point() + stat_smooth(method = "lm", col = "red") + labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
  "Intercept =",signif(fit$coef[[1]],5 ),
  " Slope =",signif(fit$coef[[2]], 5),
  " P =",signif(summary(fit)$coef[2,4], 5)))
}


ggplotRegression(lm(Richness~Age,data=df1)) + theme_bw()


