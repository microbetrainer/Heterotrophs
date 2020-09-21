library(tidyverse)
library(readxl)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads cultures that Sean sampled 
#he sent this excel sheet to me over slack
env <- read_excel("2019_06_13_CC_DNAYields.xlsx", col_types = "text") %>% select(CULTURE)


##from comparingSampleIDs.R
#1213 is 1313
#WH7803 is WH7801
#SYN1320 is SYN1220
#9201.2 is 9311

#AS9601 is correct ID rather than ASN9601
#SYN9503 is correct ID rather than S9503

#seqG$sample <- ifelse(seqG$sample == "8102", "SYNCLADEXVI", seqG$sample)


#fixes sample IDs that were mismatched
env <- env %>% mutate(CULTURE = ifelse(CULTURE == "1213", "1313", CULTURE))
env <- env %>% mutate(CULTURE = ifelse(CULTURE == "WH7801", "WH7803", CULTURE))
env <- env %>% mutate(CULTURE = ifelse(CULTURE == "SYN1320", "SYN1220", CULTURE))
env <- env %>% mutate(CULTURE = ifelse(CULTURE == "9201.2", "9311", CULTURE))
env <- env %>% mutate(CULTURE = ifelse(CULTURE == "9201.1", "9201", CULTURE))
env <- env %>% mutate(CULTURE = ifelse(CULTURE == "ASN9601", "AS9601", CULTURE))
env <- env %>% mutate(CULTURE = ifelse(CULTURE == "S9503", "SYN9503", CULTURE))
env <- env %>% mutate(CULTURE = ifelse(CULTURE == "8102", "SYNCLADEXVI", CULTURE))

#excludes JW3 and 1223
env <- env %>% filter(CULTURE != "JW3")
env <- env %>% filter(CULTURE != "1223")

#correct number of samples
env %>% nrow()

#makes a variable for whether sample is a Pro or Syn culture
env <- env %>% mutate(cultureType = ifelse(str_detect(CULTURE, "SYN|WH|9220|S9503|8102|RS9916"), "Synechococcus", "Prochlorococcus"))

env %>% group_by(cultureType) %>% summarize(n = n())
env %>% filter(is.na(cultureType))


##########Pro############


#loads metadata for Pro cultures
#from cleanUpCulture.R (checked)
meta <- read_csv("cultureCollectionCleanedUp.csv")

head(meta)
meta %>% nrow()
meta %>% filter(!is.na(NAME)) %>% nrow()
meta %>% filter(is.na(NAME))

#gets rid of unnecessary variables
colnames(meta)
meta <- meta %>% select(NAME:`METHOD (for isolation)`, Notes:cruise_id)

#makes a variable with the culture name without "Ax" or "ax"0
meta <- meta %>% mutate(nameNoAXLabel = str_replace(NAME, "Ax|ax", ""))

meta %>% filter(!is.na(NAME)) %>% nrow()
meta %>% filter(!is.na(nameNoAXLabel)) %>% nrow()
meta %>% filter(is.na(NAME))
meta %>% filter(is.na(nameNoAXLabel))

#MIT0917ax was accidentally preserved again
#meta %>% distinct(NAME, nameNoAXLabel) %>% View()
meta %>% group_by(NAME) %>% summarize(n = n()) %>% filter(n > 1)

#these are the rows that have differences in variables between the nonaxenic and axenic cultures
samplesWithDiffs <- meta %>% select(-NAME) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1)
#meta %>% semi_join(samplesWithDiffs, by = c("nameNoAXLabel")) %>% View()

###the differences in notes between axenic and nonaxenic strains don't matter
samplesWithDiffs <- meta %>% select(-c(NAME, Notes)) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1)
#meta %>% semi_join(samplesWithDiffs, by = c("nameNoAXLabel")) %>% View()


###when there are differences in method of isolation between nonaxenic and axenic strains, 
###I can't tell that the method of isolation should be the same

meta %>% select(-c(NAME, Notes, `METHOD (for isolation)`)) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1)

#MED4 and MIT0917 have two different isolators
meta %>% select(-c(NAME, Notes, `METHOD (for isolation)`, ISOLATOR)) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1)

#MIT0801 has different dates of isolation
meta %>% select(-c(NAME, Notes, `METHOD (for isolation)`, `DATE ISOLATED`)) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1)

#PAC1 has different isolation locations and different dates of isolation and depths and latitude and longitude
meta %>% select(-c(NAME, Notes, `METHOD (for isolation)`, `PLACE OF ORIGIN`, `DATE ISOLATED`, DEPTH, latTotal, lonTotal)) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1)

###also there is more data for PAC1 on the website then there is in Culture Collection.xlsx
###https://chisholmlab.mit.edu/cultures/prochlorococcus-stock-cultures
###also why is the PAC1 date so messed up?????
###is pac1ax actually eqpac1????


###these are the problem ones: see above
samplesWithDiffs <- meta %>% select(-c(NAME, Notes, `METHOD (for isolation)`)) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1)
#meta %>% semi_join(samplesWithDiffs, by = c("nameNoAXLabel")) %>% View()



###excluding the ones that are not in Sean samples

meta %>% select(-NAME) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1)

env %>% filter(str_detect(CULTURE, "9601"))
env %>% filter(str_detect(CULTURE, "MED4"))
#no
env %>% filter(str_detect(CULTURE, "0801"))
#no
env %>% filter(str_detect(CULTURE, "0917"))
env %>% filter(str_detect(CULTURE, "1205"))
env %>% filter(str_detect(CULTURE, "1214"))
#no
env %>% filter(str_detect(CULTURE, "1223"))
env %>% filter(str_detect(CULTURE, "1304"))
#no
env %>% filter(str_detect(CULTURE, "1314"))
env %>% filter(str_detect(CULTURE, "9215"))
env %>% filter(str_detect(CULTURE, "9515"))
#no
env %>% filter(str_detect(CULTURE, "PAC"))


nrow(meta)
meta %>% filter(is.na(nameNoAXLabel))
meta %>% filter(nameNoAXLabel == "MIT0801" | nameNoAXLabel == "MIT0917" | nameNoAXLabel == "MIT1223" | nameNoAXLabel == "MIT1314" | nameNoAXLabel == "PAC1")
meta %>% filter(nameNoAXLabel == "MIT0801" | nameNoAXLabel == "MIT0917" | nameNoAXLabel == "MIT1223" | nameNoAXLabel == "MIT1314" | nameNoAXLabel == "PAC1") %>% distinct(NAME)
meta <- meta %>% filter(nameNoAXLabel != "MIT0801" & nameNoAXLabel != "MIT0917" & nameNoAXLabel != "MIT1223" & nameNoAXLabel != "MIT1314" & nameNoAXLabel != "PAC1")
nrow(meta)



###these are the real problem ones: see above
samplesWithDiffs <- meta %>% select(-c(NAME, Notes, `METHOD (for isolation)`)) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1)
#meta %>% semi_join(samplesWithDiffs, by = c("nameNoAXLabel")) %>% View()



meta %>% filter(is.na(ISOLATOR)) %>% nrow()

#for MED4, the lab website (https://chisholmlab.mit.edu/cultures/prochlorococcus-stock-cultures) says the isolator was D. Vaulot, F. Partensky	
#so I am going to trust this
meta <- meta %>% mutate(ISOLATOR = ifelse(nameNoAXLabel == "MED4", "D. Vaulot, F. Partensky", ISOLATOR))

meta %>% filter(is.na(ISOLATOR)) %>% nrow()
meta %>% filter(is.na(nameNoAXLabel))

#now, none of the cultures that Sean sampled have problems
meta %>% select(-c(NAME, Notes, `METHOD (for isolation)`)) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1)



###make culture names match between env and meta

env %>% filter(is.na(CULTURE))

env <- env %>% mutate(CULTURE_originalName = CULTURE)

env %>% filter(is.na(CULTURE_originalName))

#env %>% filter(cultureType == "Prochlorococcus") %>% anti_join(meta, by = c("CULTURE" = "NAME")) %>% View()
env %>% filter(str_detect(CULTURE, "^9"))
env <- env %>% mutate(CULTURE = ifelse(str_detect(CULTURE, "^9") & cultureType == "Prochlorococcus", str_c("MIT", CULTURE), CULTURE))
#env %>% filter(cultureType == "Prochlorococcus") %>% anti_join(meta, by = c("CULTURE" = "NAME")) %>% View()
env %>% filter(str_detect(CULTURE, "^0"))
env <- env %>% mutate(CULTURE = ifelse(str_detect(CULTURE, "^0"), str_c("MIT", CULTURE), CULTURE))
env %>% filter(cultureType == "Prochlorococcus") %>% anti_join(meta, by = c("CULTURE" = "NAME"))
env %>% filter(str_detect(CULTURE, "^1"))
env <- env %>% mutate(CULTURE = ifelse(str_detect(CULTURE, "^1"), str_c("MIT", CULTURE), CULTURE))
env %>% filter(cultureType == "Prochlorococcus") %>% anti_join(meta, by = c("CULTURE" = "NAME"))
#from Culture Collection.xlsx:
#MIT0604	C12B
#MIT0604ax	C12B
env <- env %>% mutate(CULTURE = ifelse(CULTURE == "C12B", "MIT0604", CULTURE))
env %>% filter(cultureType == "Prochlorococcus") %>% anti_join(meta, by = c("CULTURE" = "NAME"))


env %>% filter(is.na(CULTURE))
env %>% filter(is.na(CULTURE_originalName))
env %>% filter(is.na(cultureType))

#gets just rows in meta that match to Pro cultures that Sean sampled
meta <- meta %>% semi_join(env, by = c("NAME" = "CULTURE"))


#loads axenic Pro tab of Culture Collection.xlsx
metaPro_axenic <- read_excel("~/Dropbox (MIT)/Lab Files/Cultures/Culture Collection.xlsx", sheet = "Axenic Pro only", col_types = "text")

#gets rid of unnecessary variables
colnames(metaPro_axenic)
metaPro_axenic <- metaPro_axenic %>% select(-c(Reference:`Axenicity confirmed`))

#makes a variable with the culture name without "Ax" or "ax"
metaPro_axenic <- metaPro_axenic %>% mutate(nameNoAXLabel = str_replace(NAME, "Ax|ax", ""))

metaPro_axenic %>% filter(is.na(NAME))
metaPro_axenic %>% filter(is.na(nameNoAXLabel))

metaPro_axenic %>% filter(NAME %in% c("B7", "DV", "JW2", "JW4", "JW7"))
metaPro_axenic <- metaPro_axenic %>% semi_join(meta, by = c("nameNoAXLabel"))

#AS9601 is the same
#C9B is the same
#MED4 in meta is correct 
#MIT9202 is the same
#MIT9215 in meta is correct
#MIT9301 is the same
#MIT9303 is the same
#MIT9312 is the same
#MIT9313 is the same 
#MIT9515 in meta is correct
#NATL1A is the same 
#SB in meta is correct 
#MIT0604 in meta is correct
#MIT1205 in meta is correct
#MIT1214 in meta is correct 
#MIT1304: the method should be flow sorted then dilution to extinction I think 
#change for the other sort 4 as well

meta %>% filter(is.na(latTotal)) %>% nrow()
meta %>% filter(is.na(lonTotal)) %>% nrow()


##########Syn############


#loads syn culture metadata
metaSyn <- read_excel("~/Dropbox (MIT)/Lab Files/Cultures/Culture Collection.xlsx", sheet = "Syn")

#laods axenic syn culture metadata
metaSyn_axenic <- read_excel("~/Dropbox (MIT)/Lab Files/Cultures/Culture Collection.xlsx", sheet = "Axenic Syn only", col_types = "text")

#gets rid of unnecessary variables
colnames(metaSyn)
metaSyn <- metaSyn %>% select(-c(Reference:need_pacbio))

colnames(metaSyn_axenic)
metaSyn_axenic <- metaSyn_axenic %>% select(-c(Reference:`Axenicity confirmed`))

metaSyn_axenic %>% filter(str_detect(NAME, "9916"))

#binds the two syn excel sheets
metaSyn <- bind_rows(metaSyn, metaSyn_axenic)

metaSyn %>% filter(is.na(NAME))

#makes a variable with the culture name without "Ax" or "ax"
metaSyn <- metaSyn %>% mutate(nameNoAXLabel = str_replace(NAME, "Ax|ax", ""))

metaSyn %>% filter(is.na(NAME))
metaSyn %>% filter(is.na(nameNoAXLabel))

#these are the rows that have differences in variables between the nonaxenic and axenic cultures
samplesWithDiffs <- metaSyn %>% select(-NAME) %>% group_by(nameNoAXLabel) %>% distinct() %>% summarize(n = n()) %>% filter(n > 1) 
#metaSyn %>% semi_join(samplesWithDiffs, by = c("nameNoAXLabel")) %>% View()

samplesWithDiffs

#Sean did not sample WH7803
env %>% filter(cultureType == "Synechococcus") %>% filter(str_detect(CULTURE, "7803"))

#the differences don't matter because I have to check all of the dates anyway
#metaSyn %>% semi_join(samplesWithDiffs, by = c("nameNoAXLabel")) %>% filter(nameNoAXLabel != "WH7803") %>% View()


###make culture names match between env and metaSyn

env %>% filter(cultureType == "Synechococcus") %>% anti_join(metaSyn, by = c("CULTURE" = "NAME"))
env %>% filter(cultureType == "Synechococcus") %>% semi_join(metaSyn, by = c("CULTURE" = "NAME"))

env %>% filter(is.na(CULTURE))
env %>% filter(is.na(cultureType))

env <- env %>% mutate(CULTURE = ifelse(cultureType == "Synechococcus", str_replace(CULTURE, "SYN", "MIT S"), CULTURE))

env %>% filter(is.na(CULTURE))
env %>% filter(is.na(cultureType))

env %>% filter(cultureType == "Synechococcus") %>% anti_join(metaSyn, by = c("CULTURE" = "NAME"))

env <- env %>% mutate(CULTURE = ifelse(CULTURE == "MITS9220", "MIT S9220", CULTURE))

env %>% filter(is.na(CULTURE))
env %>% filter(is.na(cultureType))

env %>% filter(cultureType == "Synechococcus") %>% anti_join(metaSyn, by = c("CULTURE" = "NAME"))
env <- env %>% mutate(CULTURE = ifelse(CULTURE == "MIT S1220", "MIT1220", CULTURE))

env %>% filter(is.na(CULTURE))
env %>% filter(is.na(cultureType))

env %>% filter(cultureType == "Synechococcus") %>% anti_join(metaSyn, by = c("CULTURE" = "NAME"))

metaSyn <- metaSyn %>% semi_join(env, by = c("NAME" = "CULTURE"))


#replaces NA with "" so I can perform stringr commands
metaSyn[is.na(metaSyn)] <- ""





###fixes ecotypes

#fixes Syn 5.1 ID
metaSyn %>% distinct(`ECOTYPE (subcluster)`)

metaSyn %>% filter(is.na(`ECOTYPE (subcluster)`))
metaSyn %>% filter(`ECOTYPE (subcluster)` == "") %>% nrow()

metaSyn$`ECOTYPE (subcluster)` <- str_replace(metaSyn$`ECOTYPE (subcluster)`, "5.0999999999999996", "5.1")

metaSyn %>% distinct(`ECOTYPE (subcluster)`)
metaSyn %>% filter(is.na(`ECOTYPE (subcluster)`))
metaSyn %>% filter(`ECOTYPE (subcluster)` == "") %>% nrow()

metaSyn %>% distinct(ECOTYPE)
metaSyn <- metaSyn %>% select(-ECOTYPE)

#switches ecotypes and clades
metaSyn %>% filter(NAME == "MIT1220")

metaSyn %>% filter(is.na(`ECOTYPE (subcluster)`))
metaSyn %>% filter(is.na(NAME))
metaSyn %>% filter(is.na(CLADE))
metaSyn %>% filter(`ECOTYPE (subcluster)` == "") %>% nrow()
metaSyn %>% filter(CLADE == "") %>% nrow()
metaSyn %>% filter(NAME == "") %>% nrow()

metaSyn %>% filter(NAME == "MIT1220")
metaSyn <- metaSyn %>% mutate(`ECOTYPE (subcluster)` = ifelse(NAME == "MIT1220", "5.1", `ECOTYPE (subcluster)`))
metaSyn <- metaSyn %>% mutate(CLADE = ifelse(NAME == "MIT1220", "IIh", CLADE))
metaSyn %>% filter(NAME == "MIT1220")

metaSyn %>% filter(is.na(`ECOTYPE (subcluster)`))
metaSyn %>% filter(is.na(NAME))
metaSyn %>% filter(is.na(CLADE))
metaSyn %>% filter(`ECOTYPE (subcluster)` == "") %>% nrow()
metaSyn %>% filter(CLADE == "") %>% nrow()
metaSyn %>% filter(NAME == "") %>% nrow()


###fixes depth 

metaSyn %>% distinct(DEPTH)
metaSyn %>% filter(is.na(DEPTH))
metaSyn %>% filter(DEPTH == "") %>% nrow()

metaSyn <- metaSyn %>% mutate(DEPTH = ifelse(DEPTH == "surface", "0", DEPTH))
metaSyn %>% distinct(DEPTH)
metaSyn$DEPTH <- parse_number(metaSyn$DEPTH)
metaSyn %>% distinct(DEPTH)
metaSyn %>% filter(is.na(DEPTH)) %>% nrow()
metaSyn %>% filter(DEPTH == "") %>% nrow()



###fixes notes

metaSyn %>% distinct(Notes)
metaSyn %>% filter(Notes == "MIT S1220 ?")
metaSyn$Notes <- NA


###fixes isolators

metaSyn %>% filter(is.na(ISOLATOR))
metaSyn %>% filter(ISOLATOR == "")


metaSyn %>% distinct(ISOLATOR)
metaSyn <- metaSyn %>% mutate(ISOLATOR = ifelse(ISOLATOR == "Liz Mann", "L. Mann", ISOLATOR))
metaSyn %>% distinct(ISOLATOR)
metaSyn <- metaSyn %>% mutate(ISOLATOR = ifelse(ISOLATOR == "P. Berube, J. Becker", "P.M. Berube, J. Becker", ISOLATOR))
metaSyn %>% distinct(ISOLATOR)
metaSyn <- metaSyn %>% mutate(ISOLATOR = ifelse(ISOLATOR == "L.Brand", "L. Brand", ISOLATOR))
metaSyn %>% distinct(ISOLATOR)
metaSyn <- metaSyn %>% mutate(ISOLATOR = ifelse(ISOLATOR == "J.Waterbury", "J. Waterbury", ISOLATOR))
metaSyn %>% distinct(ISOLATOR)
metaSyn <- metaSyn %>% mutate(ISOLATOR = ifelse(ISOLATOR == "F.Valois", "F. Valois", ISOLATOR))
metaSyn %>% distinct(ISOLATOR)
#metaSyn %>% select(ISOLATOR) %>% bind_rows(meta %>% select(ISOLATOR)) %>% distinct(ISOLATOR) %>% arrange(ISOLATOR) %>% View()

metaSyn %>% filter(is.na(ISOLATOR))
metaSyn %>% filter(ISOLATOR == "")



###fixes lat, longs 

metaSyn %>% distinct(COORDINATES)
metaSyn %>% filter(is.na(COORDINATES))
metaSyn %>% filter(COORDINATES == "") %>% nrow()

#woods hole = 41.5265° N, 70.6731° W
metaSyn <- metaSyn %>% mutate(COORDINATES = ifelse(COORDINATES == "Woods Hole", "41.5265° N, 70.6731° W", COORDINATES))
metaSyn %>% distinct(COORDINATES)
metaSyn %>% filter(COORDINATES == "") %>% nrow()

#makes lat and lon variables from the Coordinates character strings
metaSyn <- separate(metaSyn, COORDINATES, c("lat", "lon"), sep = "[;,]", remove = FALSE)
metaSyn %>% distinct(COORDINATES, lat, lon)
metaSyn <- metaSyn %>% mutate(lat = ifelse(is.na(lon), NA, lat))
metaSyn %>% distinct(COORDINATES, lat, lon)
metaSyn <- metaSyn %>% mutate(lat = ifelse(is.na(lon), str_extract(COORDINATES, ".*[NS]"), lat))
metaSyn %>% distinct(COORDINATES, lat, lon)
metaSyn <- metaSyn %>% mutate(lon = ifelse(is.na(lon), str_extract(COORDINATES, "(?<=[NS]).*"), lon))
metaSyn %>% distinct(COORDINATES, lat, lon)
metaSyn <- metaSyn %>% mutate(latDirection = str_extract(lat, "S|N"))
metaSyn <- metaSyn %>% mutate(lonDirection = str_extract(lon, "E|W"))
metaSyn %>% distinct(COORDINATES, lat, lon, latDirection, lonDirection)
metaSyn <- metaSyn %>% mutate(lat = str_replace(lat, "\\s*[NS]", ""))
metaSyn %>% distinct(COORDINATES, lat, lon, latDirection, lonDirection)
metaSyn <- metaSyn %>% mutate(lon = str_replace(lon, "\\s*[EW]", ""))
metaSyn %>% distinct(COORDINATES, lat, lon, latDirection, lonDirection)

metaSyn <- metaSyn %>% mutate(latMin = str_extract(lat, "[\\.0-9]*(?=’)"))
metaSyn <- metaSyn %>% mutate(lonMin = str_extract(lon, "[\\.0-9]*(?=’)"))
metaSyn %>% distinct(COORDINATES, lat, latMin, lon, lonMin, latDirection, lonDirection)
metaSyn <- metaSyn %>% mutate(latMin = as.numeric(latMin)/60)
metaSyn <- metaSyn %>% mutate(lonMin = as.numeric(lonMin)/60)
metaSyn %>% distinct(COORDINATES, lat, latMin, lon, lonMin, latDirection, lonDirection)

metaSyn$lat <- parse_number(metaSyn$lat)
metaSyn$lon <- parse_number(metaSyn$lon)
metaSyn %>% distinct(COORDINATES, lat, latMin, lon, lonMin, latDirection, lonDirection)

metaSyn$latMin <- ifelse(is.na(metaSyn$latMin), 0, metaSyn$latMin)
metaSyn$lonMin <- ifelse(is.na(metaSyn$lonMin), 0, metaSyn$lonMin)
metaSyn %>% distinct(COORDINATES, lat, latMin, lon, lonMin, latDirection, lonDirection)

metaSyn$latTotal = metaSyn$lat + metaSyn$latMin
metaSyn$lonTotal = metaSyn$lon + metaSyn$lonMin
metaSyn %>% distinct(COORDINATES, lat, latMin, latTotal, lon, lonMin, lonTotal, latDirection, lonDirection)

metaSyn$latDirection <- ifelse(is.na(metaSyn$latDirection), "NA", metaSyn$latDirection)
metaSyn$lonDirection <- ifelse(is.na(metaSyn$lonDirection), "NA", metaSyn$lonDirection)
metaSyn %>% distinct(COORDINATES, lat, latMin, latTotal, lon, lonMin, lonTotal, latDirection, lonDirection)

metaSyn <- metaSyn %>% mutate(latTotal = ifelse(latDirection == "S", (-1)*latTotal, latTotal))
metaSyn <- metaSyn %>% mutate(lonTotal = ifelse(lonDirection == "W", (-1)*lonTotal, lonTotal))
metaSyn %>% distinct(COORDINATES, latTotal, lonTotal)

metaSyn %>% filter(COORDINATES == "") %>% nrow()

#gets rid of unnecessary variables
metaSyn <- metaSyn %>% select(-c(COORDINATES, latMin, lonMin, lat, lon, latDirection, lonDirection))

metaSyn %>% filter(latTotal == "") %>% nrow()
metaSyn %>% filter(lonTotal == "") %>% nrow()

metaSyn %>% filter(is.na(latTotal)) %>% nrow()
metaSyn %>% filter(is.na(lonTotal)) %>% nrow()




###fixes place of origin

metaSyn %>% distinct(`PLACE OF ORIGIN`)

#if the x coordinate is less than 0, and the place of origin is Atlantic, it should be South Atlantic
metaSyn %>% filter(str_detect(`PLACE OF ORIGIN`, "Atlantic")) %>% filter(latTotal < 0) %>% distinct(`PLACE OF ORIGIN`)

#if the x coordinate is greater than 0, and the place of origin is Atlantic, it should be North Atlantic
metaSyn %>% filter(str_detect(`PLACE OF ORIGIN`, "Atlantic")) %>% filter(latTotal > 0) %>% distinct(`PLACE OF ORIGIN`)

metaSyn %>% filter(is.na(`PLACE OF ORIGIN`))
metaSyn %>% filter(`PLACE OF ORIGIN` ==  "") %>% nrow()
metaSyn %>% filter(is.na(latTotal)) %>% nrow()
metaSyn %>% filter(is.na(lonTotal)) %>% nrow()

metaSyn$`PLACE OF ORIGIN` <- ifelse((str_detect(metaSyn$`PLACE OF ORIGIN`, "Atlantic") & (metaSyn$latTotal > 0)), "North Atlantic", metaSyn$`PLACE OF ORIGIN`)

metaSyn %>% distinct(`PLACE OF ORIGIN`)
metaSyn %>% distinct(latTotal, lonTotal, `PLACE OF ORIGIN`)

metaSyn %>% filter(is.na(`PLACE OF ORIGIN`)) %>% nrow()
metaSyn %>% filter(`PLACE OF ORIGIN` ==  "") %>% nrow()
metaSyn %>% filter(is.na(latTotal)) %>% nrow()
metaSyn %>% filter(is.na(lonTotal)) %>% nrow()



###merges pro and syn datasets 

meta <- meta %>% mutate(taxa = "Pro")
metaSyn <- metaSyn %>% mutate(taxa = "Syn")

colnames(meta)
colnames(metaSyn)
meta <- meta %>% select(-c(`DATE ISOLATED`, nameNoAXLabel))
metaSyn %>% distinct(`METHOD (axenic method?)`)
metaSyn <- metaSyn %>% select(-c(`DATE ISOLATED`, nameNoAXLabel, `METHOD (axenic method?)`))
colnames(metaSyn)[3]
colnames(metaSyn)[3] <- "ECOTYPE"

meta <- bind_rows(meta, metaSyn)



###fixes variables

meta %>% distinct(`PLACE OF ORIGIN`) %>% arrange(`PLACE OF ORIGIN`)

#Sean says "dilution-to-extinction" is equivalent to "filtered (dilution to extinction)"
meta %>% distinct(`METHOD (for isolation)`)

meta %>% filter(`METHOD (for isolation)` == "" | is.na(`METHOD (for isolation)`)) %>% nrow()

meta$`METHOD (for isolation)` <- ifelse(meta$`METHOD (for isolation)` == "Dilution-to-extinction", "Filtered (dilution to extinction)", meta$`METHOD (for isolation)`)
meta %>% distinct(`METHOD (for isolation)`)

meta %>% filter(`METHOD (for isolation)` == "" | is.na(`METHOD (for isolation)`)) %>% nrow()

#gets rid of this note because it isn't important
meta %>% distinct(Notes)
meta %>% filter(is.na(Notes)) %>% nrow()
meta$Notes <- ifelse(meta$Notes == "Axenic LLII/III strains are rare in culture", "", meta$Notes)
meta %>% distinct(Notes)
meta %>% filter(is.na(Notes)|Notes == "") %>% nrow()

meta %>% distinct(Notes)
meta %>% filter(Notes == "Sequenced using PacBio")
meta <- meta %>% mutate(Notes = ifelse(Notes == "Sequenced using PacBio", NA, Notes))
meta %>% filter(is.na(Notes)|Notes == "") %>% nrow()

meta[meta == ""] <- NA
meta %>% filter(is.na(Notes)) %>% nrow()

meta %>% distinct(Notes)



#loads dates cultures were isolated 
#I edited this so that if the day was missing, I made it the first of the month
#if the month was also missing, I made it the first day of the year
##I have checked all of these excel documents multiple times
dateEdited <- read_excel("Culture Collection_EditedForDate.xlsx", sheet = 1, col_types = c("text", "text", "date")) %>% select(1,3)
dateEdited2 <- read_excel("Culture Collection_EditedForDate.xlsx", sheet = 2, col_types = c("text", "text", "date")) %>% select(1,3)
dateEdited3 <- read_excel("Culture Collection_EditedForDate.xlsx", sheet = 3, col_types = c("text", "text", "date")) %>% select(1,3)
dateEdited4 <- read_excel("Culture Collection_EditedForDate.xlsx", sheet = 4, col_types = c("text", "text", "date")) %>% select(1,3)

dateEdited <- bind_rows(dateEdited, dateEdited2, dateEdited3, dateEdited4) 

head(dateEdited)

dateEdited %>% filter(is.na(NAME))

dateEdited <- dateEdited %>% mutate(nameNoAXLabel = str_replace(NAME, "Ax|ax", ""))

dateEdited %>% filter(is.na(NAME))
dateEdited %>% filter(is.na(nameNoAXLabel))

probs <- dateEdited %>% group_by(nameNoAXLabel) %>% distinct(`DATE ISOLATED_Edited`) %>% summarize(n = n()) %>% filter(n > 1)
dateEdited %>% semi_join(probs, by = c("nameNoAXLabel")) %>% arrange(nameNoAXLabel)
dateEdited %>% semi_join(probs, by = c("nameNoAXLabel")) %>% semi_join(env, by = c("nameNoAXLabel" = "CULTURE")) %>% arrange(nameNoAXLabel)
meta %>% filter(NAME == "MIT1304")
dateEdited %>% filter(NAME == "MIT1304")

#I made WH7803 have 1978-07-01 because this is what WH7803ax has
#I made WH8102 have 1981-03-15 because this is what WH8102ax has
dateEdited %>% filter(NAME %in% c("WH7803", "WH8102"))

meta %>% anti_join(dateEdited, by = c("NAME"))

meta %>% group_by(NAME) %>% summarize(n = n()) %>% filter(n > 1)
meta %>% left_join(dateEdited, by = c("NAME")) %>% group_by(NAME) %>% summarize(n = n()) %>% filter(n > 1)
dateEdited %>% filter(NAME == "MIT1220")
dateEdited %>% filter(NAME == "MIT1220") %>% distinct()
dateEdited %>% nrow()
dateEdited %>% distinct() %>% nrow()
dateEdited <- dateEdited %>% distinct()

head(meta)

meta %>% filter(NAME == "B7")

b7 <- data.frame(NAME = "B7")

nrow(meta)
meta <- bind_rows(meta, b7)
nrow(meta)

meta %>% filter(NAME == "B7")

nrow(meta)
meta <- meta %>% left_join(dateEdited, by = c("NAME"))
nrow(meta)

meta <- meta %>% select(-nameNoAXLabel)




#from https://docs.google.com/spreadsheets/d/1IhEKygYDi8xrM2J7AcBmYSqkw3SxeH6F-9_aqWxJ5dI/edit#gid=1447452535:
#150NLHA	Prochlorococcus sp. 150NLHA				LLIV	Isolate	150	22.75	-158
#150SLHA	Prochlorococcus sp. 150SLHA				LLIV	Isolate	150	22.75	-158
#150SLHB	Prochlorococcus sp. 150SLHB				LLIV	Isolate	150	22.75	-158


#from Jessie thesis table 2.1:
#150NLHA LLIV LL to HL transition
#150SLHA LLIV LL to HL transition
#150SLHB LLIV LL to HL transition


#from image of tubes:
#JW2 = 150NLHB
#JW4 = 150SLHB
#JW7 = 150NLLBHL sort 4


###from getSequencesForTree_noContamination_noProSyn_withoutSample1223.R:

#col2 is in JW7
#0701, 0702, 0703, 1205, 9303, b5, c4, c8 are LLIV
#1300 is LLVII
#I'm not sure what clade b7 is

#col20 is in JW4
#1201, 1227, 9313 are LLIV
#1304 is LLII/III

#col42 is in JW2
#1300 is LLVII


#JW2 = 150NLHB
#JW4 = 150SLHB is definitely LLIV (from Jessie and Thomas tables)
#JW7 = 150NLLBHL sort 4


#JW4 = 150SLHB is definitely LL to HL transition from Jessie tables
#JW2 = 150NLHB should have the same method of isolation as 150NLHA (LL to HL transition)

#JW7 = 150NLLBHL sort 4
#according to Allison: LLB to HL (so it started at LL and then transitioned to HL)
#then it says sort 4 so this culture would have been sorted
#flow sorted then LL to HL transition



#from Jessie's thesis: 
#Samples were obtained during the HOE-PhoR cruise, which took place over May 22-June 5 (e.g. del Valle
#and Karl, 2014). Samples from 150m resulting in successful isolations were taken 6-2-2013 at Station Aloha.
#Full information on enrichment conditions prepared at sea are in Supplementary Table S2.2.

jessieUnpublished <- read_excel("jessieUnpublishedStrains.xlsx", col_types = c(replicate(8,"text"), "numeric", "text", "text", "text", "date"))
jessieUnpublished

colnames(jessieUnpublished)[1:2] <- c("NAME", "ALTERNATIVE NAME")
colnames(jessieUnpublished)[6:7] <- c("latTotal", "lonTotal")
colnames(jessieUnpublished)[13] <- "DATE ISOLATED_Edited"
jessieUnpublished <- jessieUnpublished %>% mutate(taxa = "Pro")
colnames(jessieUnpublished)
colnames(meta)

jessieUnpublished %>% select(CRUISE)
meta %>% filter(CRUISE == "HOE-PhoR")
meta %>% filter(CRUISE == "HOE-PhoR") %>% distinct(cruise_id)

str(meta)
str(jessieUnpublished)

jessieUnpublished$latTotal <- as.numeric(jessieUnpublished$latTotal)
jessieUnpublished$lonTotal <- as.numeric(jessieUnpublished$lonTotal)

jessieUnpublished %>% filter(`PLACE OF ORIGIN` == "Station ALOHA/North Pacific") %>% distinct(latTotal, lonTotal)
meta %>% filter(`PLACE OF ORIGIN` == "Station ALOHA/North Pacific") %>% distinct(latTotal, lonTotal)

jessieUnpublished

meta <- bind_rows(meta, jessieUnpublished)

env %>% anti_join(meta, by = c("CULTURE" = "NAME"))

nrow(env)
env <- env %>% left_join(meta, by = c("CULTURE" = "NAME"))
nrow(env)

env %>% filter(is.na(taxa))

env %>% filter(is.na(CULTURE))

env <- env %>% mutate(CULTURE = ifelse(CULTURE == "MIT SCLADEXVI", "SYNCLADEXVI", CULTURE))

env %>% filter(is.na(taxa))

#seqG <- seqG %>% mutate(culture = ifelse(str_detect(sample, "SYN|WH|9220|S9503|8102|RS9916"), "Syn", "Pro"))
env <- env %>% mutate(taxa = ifelse(is.na(taxa) & CULTURE %in% c("SYNCLADEXVI", "RS9916"), "Syn", taxa))
env %>% filter(is.na(taxa))
env <- env %>% mutate(taxa = ifelse(is.na(taxa), "Pro", taxa))

env %>% filter(is.na(taxa))
env %>% group_by(taxa) %>% summarize(n = n())
env %>% filter(is.na(CULTURE))


###this additional data is from https://chisholmlab.mit.edu/cultures/synechococcus-stock-cultures

env %>% filter(taxa == "Syn")

#all of the data that is on https://chisholmlab.mit.edu/cultures/synechococcus-stock-cultures
#is already here
env %>% filter(CULTURE == "MIT S9220")

#MITS9502 is not in Sean's data
env %>% filter(str_detect(CULTURE, "9502"))

env %>% filter(is.na(CULTURE))
env %>% filter(CULTURE == "")
env %>% filter(is.na(`PLACE OF ORIGIN`)) %>% nrow()
env %>% filter(`PLACE OF ORIGIN` == "") %>% nrow()

#this additional data is on the website
env %>% filter(CULTURE == "MIT S9504")
env <- env %>% mutate(`PLACE OF ORIGIN` = ifelse(CULTURE == "MIT S9504", "Equatorial Pacific", `PLACE OF ORIGIN`))
env %>% filter(CULTURE == "MIT S9504")

env %>% filter(is.na(CULTURE))
env %>% filter(is.na(`PLACE OF ORIGIN`)) %>% nrow()

#MITS9505 is not in Sean's data
env %>% filter(str_detect(CULTURE, "9505"))

#this additional data is on the website
env %>% filter(CULTURE == "MIT S9508")
env <- env %>% mutate(`PLACE OF ORIGIN` = ifelse(CULTURE == "MIT S9508", "Equatorial Pacific", `PLACE OF ORIGIN`))
env %>% filter(is.na(CULTURE))
env %>% filter(is.na(`PLACE OF ORIGIN`)) %>% nrow()

env %>% filter(is.na(DEPTH)) %>% nrow()
env <- env %>% mutate(DEPTH = ifelse(CULTURE == "MIT S9508", 0, DEPTH))
env %>% filter(CULTURE == "MIT S9508")
env %>% filter(is.na(DEPTH)) %>% nrow()

env %>% filter(is.na(CULTURE))


##compares the metadata for Pro cultures to the Pro lab website info
##from https://chisholmlab.mit.edu/cultures/prochlorococcus-stock-cultures

#SS120 is the same
#SS35 is the same
#SS51	is the same
#SS52 is the same 
#SS2 is the same
#LG	is the same
#I don't think Sean sampled PAC1
#MIT9301 is the same
#MIT9302 is the same
#MIT9303 is the same
#the culture website says MIT9401 was isolated by L.R. Moore
env %>% filter(is.na(ISOLATOR)) %>% nrow()
env <- env %>% mutate(ISOLATOR = ifelse(CULTURE == "MIT9401", "L.R. Moore", ISOLATOR))
env %>% filter(is.na(ISOLATOR)) %>% nrow()
env %>% filter(is.na(CULTURE)) %>% nrow()
#MIT9311 is the same
#MIT9312 is the same
#MIT9313 is the same
#MIT9314 is the same
#NATL1A	is the same
#NATL2A	is the same
#I don't think Sean sampled MED
#MED1	is the same
#MED4	is the same
#MIT9107 is the same
#I don't think Sean sampled MIT9116
#MIT9123 is the same
#MIT9201 is the same
#MIT9202 is the same
#MIT9211 is the same
#MIT9215 is the same
#MIT9321 is the same





#the culture website says MIT9322 was isolated at 0 degrees, 16'N; 93 degrees W
env %>% filter(is.na(latTotal)) %>% nrow()
env %>% filter(is.na(lonTotal)) %>% nrow()
env <- env %>% mutate(latTotal = ifelse(CULTURE == "MIT9322", 0.26666666666, latTotal))
env <- env %>% mutate(lonTotal = ifelse(CULTURE == "MIT9322", -93, lonTotal))
env %>% filter(is.na(latTotal)) %>% nrow()
env %>% filter(is.na(lonTotal)) %>% nrow()

#MIT9515 is the same
#GP2 is the same
#SB	is the same
#AS9601	is the same
#MIT0601 dates are different
#MIT0602 dates and method of isolation (I don't think this should be changed though) are different
#MIT0603 dates are different
#MIT0604 is the same

###I asked Allison about MIT0601-4 and she pointed me to https://www.nature.com/articles/sdata201434/tables/2
###from https://www.nature.com/articles/sdata201434/tables/2:

#MIT0601	 	eMIT9211/LLII,III	Station ALOHA/North Pacific	22.75°N 158°W	125	17-Nov-2006	This work (same as culture website)
#MIT0602	 	eSS120/LLII,III	Station ALOHA/North Pacific	22.75°N 158°W	125	17-Nov-2006	This work (same as culture website)
#MIT0603	 	eSS120/LLII,III	Station ALOHA/North Pacific	22.75°N 158°W	125	17-Nov-2006	This work (same as culture website)
#MIT0604	 	eMIT9312/HLII	Station ALOHA/North Pacific	22.75°N 158°W	175	May-2006	This work (same as culture website)

#I fixed the dates for MIT0601-3 in Culture Collection_EditedForDate.xlsx


#http://hahana.soest.hawaii.edu/hot/csreports/cs181.html
#HOT-181 Chief Scientist's Cruise Report
#R/V Kilo Moana
#May 25-29, 2006

#http://hahana.soest.hawaii.edu/hot/cruises.html
#186          18 Oct - 24 Oct 06      R/V Kilo-Moana       Grabowski               186
#187           7 Nov - 11 Nov 06      R/V Kilo-Moana       Mandujano               187
#188           8 Dec - 12 Dec 06      R/V Kilo-Moana       Grabowski               188

env %>% filter(str_detect(CULTURE, "MIT060")) %>% select(CULTURE, CRUISE, cruise_id, `DATE ISOLATED_Edited`) %>% arrange(CULTURE)

#I'm not confident in the cruises or cruise IDs for these strains
env %>% filter(is.na(CRUISE)) %>% nrow()
env %>% filter(is.na(cruise_id)) %>% nrow()
env %>% filter(CULTURE %in% c("MIT0603", "MIT0601", "MIT0602", "MIT0604")) %>% nrow()
env <- env %>% mutate(CRUISE = ifelse(CULTURE %in% c("MIT0603", "MIT0601", "MIT0602", "MIT0604"), NA, CRUISE))
env <- env %>% mutate(cruise_id = ifelse(CULTURE %in% c("MIT0603", "MIT0601", "MIT0602", "MIT0604"), NA, cruise_id))
env %>% filter(is.na(CRUISE)) %>% nrow()
env %>% filter(is.na(cruise_id)) %>% nrow()

env %>% filter(is.na(CULTURE)) %>% nrow()


env %>% filter(str_detect(CULTURE, "MIT060")) %>% select(CULTURE, CRUISE, cruise_id, `DATE ISOLATED_Edited`) %>% arrange(CULTURE)

#MIT0701 is the same
#MIT0702 is the same
#MIT0703 is the same
#I don't think Sean sampled MIT0801





#from https://aem.asm.org/content/69/5/2430/figures-only table 1
#Clade no.	RCC no.h	Strain	Site of isolation	Depth (m)	Isolation datel	PUB/PEB ratio	Motil- ityi	N require- ment(s)	Obligately marinej	Isolating scientist(s)
#IX	555	RS9916a	Gulf of Aqaba 29°28′N, 34°55′E	10	22-11-99	0.7	−	NO3−, NH4+	+	N. Fuller, this study

#I added the above date for RS9916 to Culture Collection_EditedForDate.xlsx
env %>% filter(CULTURE == "RS9916")
RS9916_date <- dateEdited %>% filter(NAME == "RS9916") %>% select(NAME, `DATE ISOLATED_Edited`)
RS9916_date
colnames(RS9916_date)[1] <- "CULTURE"
nrow(env)
env <- bind_rows(env, RS9916_date)
nrow(env)
env %>% filter(CULTURE == "RS9916")
env <- env %>% mutate(cultureType = ifelse(CULTURE == "RS9916", "Synechococcus", cultureType))
env <- env %>% mutate(CULTURE_originalName = ifelse(CULTURE == "RS9916", "RS9916", CULTURE_originalName))
env <- env %>% mutate(taxa = ifelse(CULTURE == "RS9916", "Syn", taxa))
env %>% filter(CULTURE == "RS9916")
env %>% filter(is.na(cultureType))
env %>% filter(is.na(CULTURE))
env %>% filter(is.na(CULTURE_originalName))
env %>% filter(is.na(taxa))

env %>% filter(`DATE ISOLATED_Edited` == "")
env %>% filter(is.na(`DATE ISOLATED_Edited`))

nrow(env)
env <- env %>% filter(CULTURE != "RS9916" | !is.na(`DATE ISOLATED_Edited`))
nrow(env)
env %>% filter(CULTURE == "RS9916")
env %>% filter(is.na(`DATE ISOLATED_Edited`))
env %>% filter(is.na(CULTURE))


#RS9916: https://www.ncbi.nlm.nih.gov/bioproject/13557
#Synechococcus sp. RS9916. This cyanobacterium was collected at Eliat, Red Sea at a depth of 10 meters, and was isolated by enrichment
env %>% filter(CULTURE == "RS9916")
env %>% distinct(`PLACE OF ORIGIN`) %>% arrange(`PLACE OF ORIGIN`)


env %>% filter(is.na(`PLACE OF ORIGIN`)) %>% nrow()
env %>% filter(is.na(DEPTH)) %>% nrow()

env <- env %>% mutate(`PLACE OF ORIGIN` = ifelse(CULTURE == "RS9916", "Red Sea", `PLACE OF ORIGIN`))
env <- env %>% mutate(DEPTH = ifelse(CULTURE == "RS9916", 10, DEPTH))

env %>% filter(is.na(`PLACE OF ORIGIN`)) %>% nrow()
env %>% filter(is.na(DEPTH)) %>% nrow()
env %>% filter(is.na(CULTURE))


#https://gold.jgi.doe.gov/project?id=Gp0002094
#Member of marine Synechococcus cluster 5.1; member of clade IX (16S, Fuller et al., 2003). 
env %>% distinct(ECOTYPE, CLADE)

env %>% filter(is.na(ECOTYPE)) %>% nrow()
env %>% filter(is.na(CLADE)) %>% nrow()

env <- env %>% mutate(ECOTYPE = ifelse(CULTURE == "RS9916", "5.1", ECOTYPE))
env <- env %>% mutate(CLADE = ifelse(CULTURE == "RS9916", "IX", CLADE))

env %>% filter(is.na(ECOTYPE)) %>% nrow()
env %>% filter(is.na(CLADE)) %>% nrow()

env %>% filter(CULTURE == "RS9916")
env %>% filter(is.na(CULTURE))


env %>% distinct(ISOLATOR) %>% arrange(ISOLATOR)
env %>% distinct(ISOLATOR) %>% filter(str_detect(ISOLATOR, "Fuller"))

env %>% filter(is.na(ISOLATOR)) %>% nrow()

env <- env %>% mutate(ISOLATOR = ifelse(CULTURE == "RS9916", "N. Fuller", ISOLATOR))

env %>% filter(is.na(ISOLATOR)) %>% nrow()
env %>% filter(CULTURE == "RS9916")
env %>% filter(is.na(CULTURE))


env %>% filter(is.na(latTotal)) %>% nrow()
env %>% filter(is.na(lonTotal)) %>% nrow()

#29°28′N, 34°55′E
env <- env %>% mutate(latTotal = ifelse(CULTURE == "RS9916", 29 + 28/60, latTotal))
env <- env %>% mutate(lonTotal = ifelse(CULTURE == "RS9916", 34 + 55/60, lonTotal))
env %>% filter(CULTURE == "RS9916")

env %>% filter(is.na(latTotal)) %>% nrow()
env %>% filter(is.na(lonTotal)) %>% nrow()
env %>% filter(is.na(CULTURE))


env %>% filter(is.na(`METHOD (for isolation)`)) %>% nrow()

#from https://aem.asm.org/content/69/5/2430:
#Culture isolation and growth. Seawater samples were filtered through 25-mm diameter, 0.8-μm-pore-size cellulose acetate filters (Whatman) without refrigeration. Nutrients based on SN medium (74) were added at a 1:10 dilution, containing either 100 μM NaNO3 or 100 μM NH4Cl as the sole N source. Cycloheximide (final concentration, 0.5 mg ml−1) was then added to inhibit the growth of eukaryotic algae, and the samples were incubated at 25°C and 10 microeinsteins m−2 s−1 until Synechococcus cell pellets were observed (typically 1 to 2 weeks). Cell pellets were then transferred to full-strength SN medium (1 mM N source) with cycloheximide for maintenance. Clonal Synechococcus cultures were obtained by successively plating isolates on solid SN medium 3 times, as described by Brahamsha (9, 11). Diluted cells were mixed at 37°C with SN medium containing 0.6% (wt/vol) washed Bacto agar (Difco) and containing 1 mM of the appropriate N source. The mixture was then poured immediately into petri dishes before incubation at 25°C and 5 microeinsteins m−2 s−1 for 24 h then at 10 microeinsteins m−2 s−1 until colonies appeared (typically 1 to 4 weeks). Subsequent clonal Synechococcus sp. strains were maintained in SN (with 1 mM N source) at 25°C and 10 microeinsteins m−2 s−1 and transferred to fresh medium every 2 to 4 weeks.
#over slack, Sean said to classify this as filtered
env <- env %>% mutate(`METHOD (for isolation)` = ifelse(CULTURE == "RS9916", "Filtered", `METHOD (for isolation)`))
env %>% filter(CULTURE == "RS9916")

env %>% filter(is.na(`METHOD (for isolation)`)) %>% nrow()
env %>% filter(is.na(CULTURE))


##https://www.pnas.org/content/pnas/suppl/2017/06/16/1700990114.DCSupplemental/pnas.1700990114.sapp.pdf:
#Synechococcus MIT9509 5.1B CRD1 3.09 Equatorial Pacific 80 29 Yes (MIT9509) This study
#Synechococcus MIT9508 5.1B CRD1 2.50 Equatorial Pacific 8 0 No Moore et al., 2002 and this study
#Synechococcus MIT9504 5.1B CRD1 3.09 Equatorial Pacific 78 29 Yes (MIT9509) Moore et al., 2002 and this study

env %>% distinct(ECOTYPE, CLADE)

env %>% filter(CULTURE %in% c("MIT S9509", "MIT S9508", "MIT S9504"))

env %>% filter(is.na(ECOTYPE)) %>% nrow()
env %>% filter(is.na(CLADE)) %>% nrow()
env %>% filter(is.na(`PLACE OF ORIGIN`)) %>% nrow()

env <- env %>% mutate(ECOTYPE = ifelse(CULTURE %in% c("MIT S9509", "MIT S9508", "MIT S9504"), "5.1b", ECOTYPE))
env <- env %>% mutate(CLADE = ifelse(CULTURE %in% c("MIT S9509", "MIT S9508", "MIT S9504"), "CRD1", CLADE))
env <- env %>% mutate(`PLACE OF ORIGIN` = ifelse(CULTURE %in% c("MIT S9509", "MIT S9508", "MIT S9504"), "Equatorial Pacific", `PLACE OF ORIGIN`))

env %>% filter(is.na(ECOTYPE)) %>% nrow()
env %>% filter(is.na(CLADE)) %>% nrow()
env %>% filter(is.na(`PLACE OF ORIGIN`)) %>% nrow()

env %>% filter(is.na(CULTURE))

env %>% filter(CULTURE %in% c("MIT S9509", "MIT S9508", "MIT S9504"))

##https://aslopubs.onlinelibrary.wiley.com/doi/epdf/10.4319/lo.2002.47.4.0989
#MIT S9220 is the same
#MIT S9504 is the same
#MIT S9508 is the same



env %>% filter(is.na(ECOTYPE)) %>% nrow()
env %>% filter(is.na(CLADE)) %>% nrow()

#from https://www.nature.com/articles/ismej2015115:
#Clade CRD1 is a member of subcluster 5.1B (Ahlgren and Rocap, 2012; Mazard et al., 2012a), while ITS phylogeny suggests clade CRD2 is closely related to 5.1 A clades (Huang et al., 2011; Ahlgren and Rocap, 2012).
env %>% filter(CLADE == "CRD1")
env <- env %>% mutate(ECOTYPE = ifelse(CLADE == "CRD1", "5.1b", ECOTYPE))
env %>% filter(CLADE == "CRD1")

env %>% filter(is.na(ECOTYPE)) %>% nrow()
env %>% filter(is.na(CLADE)) %>% nrow()

env %>% distinct(ECOTYPE, CLADE)

#from https://www.nature.com/articles/ismej2015115:
#figure 7
env <- env %>% mutate(ECOTYPE = ifelse(CLADE == "IX", "5.1b", ECOTYPE))

env %>% distinct(ECOTYPE, CLADE) %>% arrange(ECOTYPE)

env <- env %>% mutate(ECOTYPE = ifelse(CLADE == "IIh", "5.1a", ECOTYPE))

env %>% distinct(ECOTYPE, CLADE) %>% arrange(ECOTYPE)

env %>% filter(is.na(ECOTYPE)) %>% nrow()
env %>% filter(is.na(CLADE)) %>% nrow()


env %>% filter(str_detect(CULTURE, "JW"))

env %>% filter(CRUISE == "HOE-PhoR" & !str_detect(CULTURE, "JW")) 

env %>% filter(CRUISE == "HOE-PhoR")


#from Jessie's thesis: 
#Samples were obtained during the HOE-PhoR cruise, which took place over May 22-June 5 (e.g. del Valle
#and Karl, 2014). Samples from 150m resulting in successful isolations were taken 6-2-2013 at Station Aloha.
#I replaced 6/1/13 with 6/2/13 for the HOE-PhoR cultures in Culture Collection_EditedForDate
env %>% filter(CRUISE == "HOE-PhoR" & !str_detect(CULTURE, "JW")) %>% distinct(CULTURE, DEPTH, `DATE ISOLATED_Edited`)
env %>% filter(CRUISE == "HOE-PhoR" & !str_detect(CULTURE, "JW")) %>% distinct(`DATE ISOLATED_Edited`) 

env %>% filter(is.na(CULTURE))



#for the LLIV strains that Jessie collected during the HOE-phoR cruise except for the 150S and 150N ones
#(which are not selected here), I made the isolation method filtered (dilution to extinction) because 
#this is what it said in her thesis

#from Jessie's thesis Table 2.1:
#1313 LLIV Dilution to extinction

env %>% filter(CRUISE == "HOE-PhoR")
env %>% filter(CRUISE == "HOE-PhoR" & !str_detect(CULTURE, "JW"))
env %>% filter(is.na(`METHOD (for isolation)`)) %>% nrow()
str(env$`METHOD (for isolation)`)

env %>% filter(is.na(CRUISE)) %>% nrow()
env <- env %>% mutate(CRUISE = ifelse(is.na(CRUISE), "", CRUISE))
env %>% filter(CRUISE == "") %>% nrow()

env %>% mutate(`METHOD (for isolation)` = ifelse(CRUISE == "HOE-PhoR" & !str_detect(CULTURE, "JW"), "Filtered (dilution to extinction)", `METHOD (for isolation)`)) %>% filter(is.na(`METHOD (for isolation)`)) %>% nrow()
env <- env %>% mutate(`METHOD (for isolation)` = ifelse(CRUISE == "HOE-PhoR" & !str_detect(CULTURE, "JW"), "Filtered (dilution to extinction)", `METHOD (for isolation)`))
env %>% filter(is.na(`METHOD (for isolation)`)) %>% nrow()

env %>% filter(CRUISE == "HOE-PhoR")
env <- env %>% mutate(CRUISE = ifelse(CRUISE == "", NA, CRUISE))
env %>% filter(is.na(CRUISE)) %>% nrow()

env %>% distinct(`METHOD (for isolation)`)
env %>% group_by(`METHOD (for isolation)`) %>% summarize(n = n())

env %>% filter(is.na(CULTURE))



#MIT1304's alternative name is 150NLLB sort #4 and in metaPro_axenic, it said 
#the method was flow sorted then dilution to extinction
#therefore, the method should be flow sorted then dilution to extinction
env %>% filter(CULTURE == "MIT1304")
metaPro_axenic %>% filter(NAME == "MIT1304ax")
env %>% filter(is.na(`METHOD (for isolation)`)) %>% nrow()
env$`METHOD (for isolation)` <- ifelse(env$CULTURE == "MIT1304", "Flow sorted then filtered (dilution to extinction)", env$`METHOD (for isolation)`)
env %>% filter(CULTURE == "MIT1304")
env %>% filter(is.na(`METHOD (for isolation)`)) %>% nrow()

env %>% filter(CRUISE == "HOE-PhoR")

env %>% filter(is.na(CULTURE)) %>% nrow()




####

#CoFeMUG cruise KN192-05, station 13 cruise logs:
#http://data.bco-dmo.org/jg/serv/BCO/CoFeMUG/KN192-5/event_log.html0%7Bdir=data.bco-dmo.org/jg/dir/BCO/CoFeMUG/KN192-5/,info=data.bco-dmo.org/jg/info/BCO/CoFeMUG/KN192-5/event_log%7D

env %>% filter(is.na(CRUISE)) %>% nrow()


##makes C8, B5, and C4 have CoFeMUG cruise KN192-05, station 13 for the cruise variable 
##because they were isolated at the same depth, lat, lon, and time as the cultures
##from CoFeMUG cruise KN192-05, station 13

env %>% filter(`PLACE OF ORIGIN` == "South Atlantic") 
env %>% filter(`PLACE OF ORIGIN` == "South Atlantic") %>% select(CULTURE, CRUISE, DEPTH, ISOLATOR, latTotal, lonTotal, `DATE ISOLATED_Edited`)
env %>% filter(`PLACE OF ORIGIN` == "South Atlantic") %>% distinct(DEPTH, ISOLATOR, latTotal, lonTotal, `DATE ISOLATED_Edited`)

env <- env %>% mutate(CRUISE = ifelse(CULTURE %in% c("C8", "B5", "C4"), "CoFeMUG cruise KN192-05, station 13", CRUISE))

env %>% filter(`PLACE OF ORIGIN` == "South Atlantic") 
env %>% filter(`PLACE OF ORIGIN` == "South Atlantic") %>% select(CULTURE, `PLACE OF ORIGIN`, CRUISE, DEPTH, ISOLATOR, latTotal, lonTotal, `DATE ISOLATED_Edited`)
env %>% filter(`PLACE OF ORIGIN` == "South Atlantic") %>% distinct(`PLACE OF ORIGIN`, CRUISE, DEPTH, ISOLATOR, latTotal, lonTotal, `DATE ISOLATED_Edited`)

env %>% filter(CRUISE == "CoFeMUG cruise KN192-05, station 13")
env %>% filter(CRUISE == "CoFeMUG cruise KN192-05, station 13") %>% distinct(`PLACE OF ORIGIN`, CRUISE, DEPTH, ISOLATOR, latTotal, lonTotal, `DATE ISOLATED_Edited`)

env %>% filter(is.na(CRUISE)) %>% nrow()
env %>% filter(is.na(CULTURE)) %>% nrow()



###make variable for whether sample is from inland sea, open ocean, or coastal

env %>% filter(is.na(`PLACE OF ORIGIN`)) %>% nrow()
env %>% filter(is.na(latTotal)) %>% nrow()
env %>% filter(is.na(lonTotal)) %>% nrow()
env %>% filter(is.na(CULTURE)) %>% nrow()


env %>% distinct(`PLACE OF ORIGIN`) %>% arrange(`PLACE OF ORIGIN`)
env <- env %>% mutate(oceanType = ifelse(`PLACE OF ORIGIN` == "Red Sea", "Inland sea", NA))
env <- env %>% mutate(oceanType = ifelse(`PLACE OF ORIGIN` == "Mediterranean Sea", "Inland sea", oceanType))


env <- env %>% mutate(oceanType = ifelse(`PLACE OF ORIGIN` == "Station ALOHA/North Pacific", "Open ocean", oceanType))

#woods hole is at 41.5265° N, 70.6731° W
env %>% filter(CULTURE == "WH8017") %>% select(CULTURE, latTotal, lonTotal)
env <- env %>% mutate(oceanType = ifelse(CULTURE == "WH8017", "Coastal", oceanType))

env %>% filter(`PLACE OF ORIGIN` == "Sargasso Sea") %>% distinct(latTotal, lonTotal)
env <- env %>% mutate(oceanType = ifelse(`PLACE OF ORIGIN` == "Sargasso Sea", "Open ocean", oceanType))

#env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% View()

env %>% filter(`PLACE OF ORIGIN` == "Equatorial Pacific") %>% distinct(latTotal, lonTotal)
env %>% filter(`PLACE OF ORIGIN` == "Equatorial Pacific") %>% filter(!is.na(latTotal)) %>% select(CULTURE, latTotal, lonTotal, `PLACE OF ORIGIN`, oceanType)

env <- env %>% mutate(oceanType = ifelse(`PLACE OF ORIGIN` == "Equatorial Pacific" & !is.na(latTotal), "Open ocean", oceanType))
env %>% filter(`PLACE OF ORIGIN` == "Equatorial Pacific") %>% select(CULTURE, latTotal, lonTotal, `PLACE OF ORIGIN`, oceanType)

#env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% View()

env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% distinct(`PLACE OF ORIGIN`, latTotal, lonTotal) %>% arrange(`PLACE OF ORIGIN`)

env %>% filter(`PLACE OF ORIGIN` == "Western Pacific") %>% select(CULTURE, `PLACE OF ORIGIN`, latTotal, lonTotal)

env %>% filter(CULTURE == "SB") %>% select(CULTURE, `PLACE OF ORIGIN`, latTotal, lonTotal)
env <- env %>% mutate(oceanType = ifelse(CULTURE == "SB", "Coastal", oceanType))

env %>% filter(`PLACE OF ORIGIN` == "Western Pacific") %>% select(CULTURE, `PLACE OF ORIGIN`, latTotal, lonTotal, oceanType)

env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% distinct(`PLACE OF ORIGIN`, latTotal, lonTotal) %>% arrange(`PLACE OF ORIGIN`)
env %>% filter(`PLACE OF ORIGIN` == "North Atlantic") %>% distinct(latTotal, lonTotal, oceanType)
env %>% filter(`PLACE OF ORIGIN` == "North Atlantic") %>% select(CULTURE, latTotal, lonTotal, oceanType)

env %>% filter(CULTURE == "WH7803") %>% select(CULTURE, latTotal, lonTotal, oceanType)
env <- env %>% mutate(oceanType = ifelse(CULTURE == "WH7803", "Open ocean", oceanType))

env %>% filter(CULTURE == "WH8012") %>% select(CULTURE, latTotal, lonTotal, oceanType)
env <- env %>% mutate(oceanType = ifelse(CULTURE == "WH8012", "Open ocean", oceanType))

env %>% filter(CULTURE == "NATL1A") %>% select(CULTURE, latTotal, lonTotal, oceanType)
env <- env %>% mutate(oceanType = ifelse(CULTURE == "NATL1A", "Open ocean", oceanType))

env %>% filter(CULTURE == "NATL2A") %>% select(CULTURE, latTotal, lonTotal, oceanType)
env <- env %>% mutate(oceanType = ifelse(CULTURE == "NATL2A", "Open ocean", oceanType))

env %>% filter(is.na(oceanType)) %>% filter(`PLACE OF ORIGIN` == "North Atlantic") %>% select(CULTURE, latTotal, lonTotal, oceanType)

env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% distinct(`PLACE OF ORIGIN`, latTotal, lonTotal) %>% arrange(`PLACE OF ORIGIN`)

env %>% filter(`PLACE OF ORIGIN` == "South Atlantic") %>% distinct(latTotal, lonTotal)
env <- env %>% mutate(oceanType = ifelse(`PLACE OF ORIGIN` == "South Atlantic", "Open ocean", oceanType))
env %>% filter(`PLACE OF ORIGIN` == "South Atlantic") %>% distinct(latTotal, lonTotal, oceanType)

env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% distinct(`PLACE OF ORIGIN`, latTotal, lonTotal) %>% arrange(`PLACE OF ORIGIN`)

env %>% filter(`PLACE OF ORIGIN` == "Tropical Pacific") %>% distinct(latTotal, lonTotal)
env <- env %>% mutate(oceanType = ifelse(`PLACE OF ORIGIN` == "Tropical Pacific", "Open ocean", oceanType))

env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% distinct(`PLACE OF ORIGIN`, latTotal, lonTotal) %>% arrange(`PLACE OF ORIGIN`)

env %>% filter(`PLACE OF ORIGIN` == "Arabian Sea") %>% distinct(latTotal, lonTotal)
env <- env %>% mutate(oceanType = ifelse(`PLACE OF ORIGIN` == "Arabian Sea", "Open ocean", oceanType))

env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% distinct(`PLACE OF ORIGIN`, latTotal, lonTotal) %>% arrange(`PLACE OF ORIGIN`)
env %>% filter(`PLACE OF ORIGIN` == "Gulf Stream") %>% distinct(latTotal, lonTotal)
env <- env %>% mutate(oceanType = ifelse(`PLACE OF ORIGIN` == "Gulf Stream", "Open ocean", oceanType))

env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% distinct(CULTURE, `PLACE OF ORIGIN`, latTotal, lonTotal) %>% arrange(`PLACE OF ORIGIN`)

env %>% filter(CULTURE == "WH8102") %>% select(CULTURE, `PLACE OF ORIGIN`, latTotal, lonTotal, oceanType)
env <- env %>% mutate(oceanType = ifelse(CULTURE == "WH8102", "Open ocean", oceanType))

env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% distinct(CULTURE, `PLACE OF ORIGIN`, latTotal, lonTotal) %>% arrange(`PLACE OF ORIGIN`)

env <- env %>% mutate(oceanType = ifelse(CULTURE %in% c("WH8109", "WH6501", "GP2"), "Open ocean", oceanType))

env %>% filter(is.na(oceanType) & !is.na(latTotal)) %>% nrow()
env %>% filter(!is.na(oceanType) & is.na(latTotal)) %>% nrow()

env %>% filter(oceanType == "Coastal")
env %>% filter(oceanType == "Coastal") %>% select(latTotal, lonTotal)

env %>% filter(is.na(`PLACE OF ORIGIN`)) %>% nrow()
env %>% filter(is.na(latTotal)) %>% nrow()
env %>% filter(is.na(lonTotal)) %>% nrow()
env %>% filter(is.na(CULTURE)) %>% nrow()

env %>% nrow()

env %>% distinct(oceanType)
env %>% group_by(oceanType) %>% summarize(n = n())
env %>% filter(is.na(oceanType)) %>% distinct(latTotal, lonTotal)

env %>% filter(oceanType == "Inland sea")
env %>% filter(oceanType == "Inland sea") %>% distinct(`PLACE OF ORIGIN`)

env %>% filter(is.na(oceanType))

head(env)

##fills in metadata for B7 that Sean sent me over Slack 4/10/20

#env %>% filter(CULTURE %in% c("B7", "C8", "B5")) %>% View()

env <- env %>% mutate(`PLACE OF ORIGIN` = ifelse(CULTURE == "B7", "South Atlantic", `PLACE OF ORIGIN`))
env <- env %>% mutate(CRUISE = ifelse(CULTURE == "B7", "CoFeMUG cruise KN192-05, station 13", CRUISE))
env <- env %>% mutate(CLADE = ifelse(CULTURE == "B7", "LLIV", CLADE))
env <- env %>% mutate(ECOTYPE = ifelse(CULTURE == "B7", "Low Light", ECOTYPE))
env <- env %>% mutate(DEPTH = ifelse(CULTURE == "B7", 150, DEPTH))
env <- env %>% mutate(ISOLATOR = ifelse(CULTURE == "B7", "L.R. Moore", ISOLATOR))
env <- env %>% mutate(`METHOD (for isolation)` = ifelse(CULTURE == "B7", "Filtered (dilution to extinction)", `METHOD (for isolation)`))

env %>% filter(CULTURE %in% c("B7", "C8", "B5")) %>% select(latTotal, lonTotal)

lat <- env %>% filter(CULTURE == "C8") %>% select(latTotal)
lat
lat <- lat$latTotal
lat

lon <- env %>% filter(CULTURE == "C8") %>% select(lonTotal)
lon
lon <- lon$lonTotal
lon

env <- env %>% mutate(latTotal = ifelse(CULTURE == "B7", lat, latTotal))
env <- env %>% mutate(lonTotal = ifelse(CULTURE == "B7", lon, lonTotal))
env <- env %>% mutate(oceanType = ifelse(CULTURE == "B7", "Open ocean", oceanType))

#env %>% filter(CULTURE %in% c("B7", "C8", "B5")) %>% View()

str(env$`DATE ISOLATED_Edited`)

env %>% filter(CULTURE %in% c("B7", "C8", "B5")) %>% 
  select(-c(CULTURE, CULTURE_originalName, `ALTERNATIVE NAME`)) %>% distinct()
  

nrow(env)

#sean says LG is a LLII/III rather than a HLII
#I think it was wrong in the Culture Collection.xlsx in Dropbox that I originally used
env %>% filter(CULTURE == "LG")
env %>% filter(CULTURE == "LG") %>% select(ECOTYPE, CLADE)
env <- env %>% mutate(CLADE = ifelse(CULTURE == "LG", "LLII/III", CLADE))
env %>% filter(CULTURE == "LG") %>% select(ECOTYPE, CLADE)
#env %>% filter(CULTURE == "LG") %>% View()
env %>% distinct(ECOTYPE)
env <- env %>% mutate(ECOTYPE = ifelse(CULTURE == "LG", "Low Light", ECOTYPE))
env %>% filter(CULTURE == "LG") %>% select(ECOTYPE, CLADE)
#env %>% filter(CULTURE == "LG") %>% View()

write_csv(env, "cleanedUpCultureEnvVariablesForPCA.csv")




env %>% filter(CULTURE == "JW3")

forSean <- env

forSean <- forSean %>% mutate(`X Coordinate` = ifelse(CULTURE == "JW3", "22.75", `X Coordinate`))
forSean <- forSean %>% mutate(`Y Coordinate` = ifelse(CULTURE == "JW3", "-158.0", `Y Coordinate`))

forSean <- forSean %>% select(CULTURE, `X Coordinate`, `Y Coordinate`)

write_csv(forSean, "2019_06_13_CC_DNAYieldsLatLong_withJW.csv")






