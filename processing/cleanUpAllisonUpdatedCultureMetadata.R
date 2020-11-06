library(tidyverse)
library(readxl)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads updated metadata for cultures that Allison sent over email on 5/15/20
meta <- read_excel("Supplementary Tables_AC.xlsx", sheet = 1)

#makes first row the variable names
colnames(meta)
meta[1,]
colnames(meta) <- meta[1,]
colnames(meta)
meta[1,]
meta <- meta[-1,]

#gets rid of rows with all NAs
#meta %>% filter(is.na(`Culture Name`)) %>% View()
meta <- meta %>% filter(!is.na(`Culture Name`))

#gets rid of row that just has title
#meta %>% filter(`Culture Name` == "Enrichment (not unialgal) cultures") %>% View()
meta <- meta %>% filter(`Culture Name` != "Enrichment (not unialgal) cultures")

str(meta)

#replaces all instances of "NA" with NA
meta[meta == "NA"]
meta[meta == "NA"] <- NA

#sean says LG is a LLII/III rather than a HLII
#I think it was wrong in the Culture Collection.xlsx in Dropbox that I originally used
meta %>% filter(`Culture Name` == "LG")
meta %>% filter(`Culture Name` == "LG") %>% select(Ecotype, Clade)
meta <- meta %>% mutate(Clade = ifelse(`Culture Name` == "LG", "LLII/III", Clade))
meta %>% filter(`Culture Name` == "LG") %>% select(Ecotype, Clade)
meta <- meta %>% mutate(Ecotype = ifelse(`Culture Name` == "LG", "Low Light", Ecotype))
meta %>% filter(`Culture Name` == "LG") %>% select(Ecotype, Clade)
#meta %>% filter(`Culture Name` == "LG") %>% View()

meta %>% nrow()
meta %>% distinct(`Culture Name`) %>% nrow()
meta %>% group_by(`Culture Name`) %>% summarize(n = n()) %>% filter(n > 1)

colnames(meta)
meta %>% distinct(Genus)

#correct number of Pro and Syn cultures
meta %>% group_by(Genus) %>% summarize(n = n())

meta %>% distinct(Ecotype)
meta %>% distinct(Clade) %>% arrange(Clade)
meta %>% distinct(`Place of origin`) %>% arrange(`Place of origin`)
meta %>% distinct(`Cruise Name`) %>% arrange(`Cruise Name`)

#makes depth numeric
str(meta$`Depth (m)`)
meta$`Depth (m)` <- as.numeric(meta$`Depth (m)`)
str(meta$`Depth (m)`)
meta %>% distinct(`Depth (m)`) %>% arrange(`Depth (m)`)

#meta %>% distinct(Isolator) %>% arrange(Isolator) %>% View()

meta %>% distinct(`Isolation Method`) %>% arrange(`Isolation Method`)
meta %>% distinct(Notes) %>% arrange(Notes)

#makes latitude numeric
str(meta$Latitude)
meta$Latitude <- as.numeric(meta$Latitude)
str(meta$Latitude)

#makes longitude numeric
str(meta$Longitude)
meta$Longitude <- as.numeric(meta$Longitude)
str(meta$Longitude)

meta %>% distinct(`Cruise ID`)

meta %>% distinct(Taxa)
meta %>% group_by(Taxa) %>% summarize(n = n())

#makes isolation date a date variable
str(meta$`Date Isolated`)
meta$`Date Isolated`
as.Date(meta$`Date Isolated`)
meta$`Date Isolated` <- as.Date(meta$`Date Isolated`)

meta %>% distinct(LonghurstCode) %>% arrange(LonghurstCode)
meta %>% distinct(LonghurstProvince) %>% arrange(LonghurstProvince)
meta %>% distinct(Biome) %>% arrange(Biome)

colnames(meta)

#makes a variable for the age of the culture at the time they were sequenced
meta <- meta %>% mutate(yearsInAge = interval(ymd(`Date Isolated`), ymd("2019-06-01 UTC"))/years(1))

#makes depth range variable
meta <- meta %>% mutate(depthRange = ifelse(`Depth (m)` <= 50, "0-50m", NA))
meta <- meta %>% mutate(depthRange = ifelse(`Depth (m)` > 50 & `Depth (m)` <= 100, "50-100m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(`Depth (m)` > 100 & `Depth (m)` <= 150, "100-150m", depthRange))
meta <- meta %>% mutate(depthRange = ifelse(`Depth (m)` > 150 & `Depth (m)` <= 200, "150-200m", depthRange))

meta %>% filter(is.na(yearsInAge)) %>% distinct(`Date Isolated`)
meta %>% filter(is.na(depthRange)) %>% distinct(`Depth (m)`)

colnames(meta)

write_csv(meta, "cleanedUpCultureEnvVariablesForPCA_ACUpdates.csv")


###make map of culture isolation locations

meta %>% filter(is.na(Latitude))
meta %>% filter(is.na(Longitude))

colnames(meta)[11]
colnames(meta)[11] <- "lat"

colnames(meta)[12]
colnames(meta)[12] <- "long"


#rounds lat and long values
meta <- meta %>% mutate(lat=floor(lat), long=floor(long)) 

meta %>% filter(is.na(lat))
meta %>% filter(is.na(long))

meta %>% group_by(`Culture Name`) %>% summarize(n = n()) %>% filter(n > 1)

#calculates the number number of cultures at each lat, long
meta <- meta %>% group_by(lat, long) %>% summarize(n = n())

meta %>% filter(is.na(lat))

#excludes NA lat, long row
meta <- meta %>% filter(!is.na(lat))

# project the data
proj="robin"
#long_0 = -90
meta_2 <- project(meta, proj)

ggplot() +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + 
  geom_point(aes(long, lat, size = n), meta_2, color = "black") + 
  labs(size = "Number of cultures")
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

ggsave("mapOfCultures_ACUpdates.png", dpi = 600, height = 10, width = 12)
ggsave("mapOfCultures_ACUpdates.svg", dpi = 600, height = 10, width = 12)


