library(tidyverse)
library(readxl)
library(ggworldmap)

setwd("~/Dropbox (MIT)/Sean16SSep2019/")

#loads cleaned up metadata for cultures
#from cleanUpEnvVariablesOfCulturesForPCA.R
meta <- read_csv("cleanedUpCultureEnvVariablesForPCA.csv")

#JW3 and 1223 have already been excluded
meta %>% filter(CULTURE == "JW3")
meta %>% filter(str_detect(CULTURE, "1223"))

colnames(meta)[13:14]
colnames(meta)[13:14] <- c("lat", "long")

# project the data
proj="robin"
#long_0 = -90
meta_2 <- project(meta, proj)

###is it okay to get rid of the long_0 in the projection??
# plot the world
ggplot() +
  geom_worldmap(proj = proj, color = "gray90", fill = "gray90") +
  geom_graticule(proj = proj, color = "gray80") +
  geom_gratframe(proj = proj, color = "gray80") +
  geom_degree(proj = proj, color = "gray80", long_by = 80) +
  theme_worldmap_light() +
  coord_equal() +
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 10)) + 
  theme(legend.key.size = unit(5, "mm")) + 
  geom_point(aes(long, lat), meta_2, size=1.5, color = "black") + 
  geom_text(aes(label = CULTURE, long, lat), meta_2, position = position_jitter(width = 4, height = 4))
#+ ggtitle("PO4 (mmol/m^3), time: 2019-03-16, depth: 92.33 m") 

str(meta$lat)
str(meta$long)

meta %>% filter(is.na(lat)) %>% nrow()
meta %>% filter(is.na(long)) %>% nrow()

#rounds lat and long values
meta <- meta %>% mutate(lat=floor(lat), long=floor(long)) 

meta %>% filter(is.na(lat)) %>% nrow()
meta %>% filter(is.na(long)) %>% nrow()

meta %>% group_by(CULTURE) %>% summarize(n = n()) %>% filter(n > 1)
meta %>% group_by(CULTURE_originalName) %>% summarize(n = n()) %>% filter(n > 1)

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

ggsave("mapOfCultures.png", dpi = 600, height = 10, width = 12)
ggsave("mapOfCultures.svg", dpi = 600, height = 10, width = 12)

