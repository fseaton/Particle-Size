## The code to read in data from CEH computer
## Saved separately to make it easier to move from computer to computer

## The particle size distributions
GMEP <- read.csv("***REMOVED***")

library(RODBC)

channel1 <- odbcConnect("***REMOVED***")

## list tables under GMEP_DERIVED
sqlTables(channel1, schema="GMEP_DERIVED")

Soil <- sqlFetch(channel1, "GMEP_DERIVED.WP8_SOILMETRICS", na.strings = "-9999.999")
summary(Soil)

Ancil <- sqlFetch(channel1, "GMEP_DERIVED.WP8_SOIL_ANCILLARY_Y1_4", na.strings = "-9999.999")
summary(Ancil)
Ancil <- Ancil[!is.na(Ancil$REP_ID),]
Ancil <- dplyr::select(Ancil, -X_COORD, -Y_COORD) #removing sensitive information

## Microbial diversity not in sql database yet so fetching from server
MIC_RICH <- read.csv("***REMOVED***")
FUNG_RICH_BLAST <- read.csv("***REMOVED***", row.names = 1)
colnames(FUNG_RICH_BLAST) <- c("REP_ID", paste0(colnames(FUNG_RICH_BLAST)[2:12], "_BLAST"))

MIC_RICH <- merge(MIC_RICH, FUNG_RICH_BLAST, by="REP_ID", all=T)
summary(MIC_RICH)


## Bacteria data 

BACT <- read.csv("***REMOVED***")


## Fungi data
BLAST_m <- read.csv("***REMOVED***")
BLAST_um <- read.csv("***REMOVED***")

# save to file just in case
rm(channel1)
save.image("Data_files.RData")
