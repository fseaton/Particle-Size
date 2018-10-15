## Calculating Hill numbers for microbial data 
# mainly for interest
library(tidyverse)
library(vegan)
library(qgraph)
library(psych)

#### bacteria data prep ####
## get data 
dim(BACT)
colnames(BACT)

## quick check to see if there are weird dominant taxa
summary(BACT$taxonomy)

BACT <- separate(BACT, col=taxonomy, sep=";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                 remove = FALSE, convert = T)
summary(as.factor(BACT$Kingdom))
summary(as.factor(BACT$Family))

BACT_ORD <- BACT

## filter to bacteria only and remove mitochondrial and chloroplast DNA
BACT_ORD <- filter(BACT_ORD, Kingdom == "k__Bacteria")
BACT_ORD <- filter(BACT_ORD, Family != " f__mitochondria")
BACT_ORD <- filter(BACT_ORD, Class != " c__Chloroplast")

## check if this worked
length(grep("mitochondria", BACT$taxonomy, ignore.case = TRUE))# 2449
length(grep("mitochondria", BACT_ORD$taxonomy, ignore.case = TRUE)) # 0
length(grep("chloroplast", BACT$taxonomy, ignore.case = TRUE))# 245
length(grep("chloroplast", BACT_ORD$taxonomy, ignore.case = TRUE)) # 0

## remove OTU_ID
rownames(BACT_ORD) <- BACT_ORD$OTU_ID
rownames(BACT_ORD)[1:10]
colnames(BACT_ORD)
BACT_ORD <- BACT_ORD[,2:439]
colnames(BACT_ORD)

## get rid of 'empty' OTUs and transpose data for vegan
BACT_ORD <- BACT_ORD[rowSums(BACT_ORD)>0,]
BACT_ORD <- t(BACT_ORD)



## removing reps with read depth less than 2000 (there's two, and the next smallest is 18873)
sort(rowSums(BACT_ORD), decreasing = T)
dim(BACT_ORD) #438 58165
BACT_ORD <- BACT_ORD[rowSums(BACT_ORD)>2000,]
dim(BACT_ORD) #436 58165



## bacteria ####
sort(rowSums(BACT_ORD))
## try rarefying to 18800 reads
dim(BACT_ORD)
min(colSums(BACT_ORD)) #1
BACT_ORD <- BACT_ORD[,colSums(BACT_ORD)>1]
dim(BACT_ORD)

## rarefaction
BACT_RARE_ARR <- array(dim=c(436,58164,10))
set.seed(250618)

for (i in 1:10){
  BACT_RARE_ARR[,,i] <- rrarefy(BACT_ORD, 18800)
}
saveRDS(BACT_RARE_ARR, file="./Bacterial rarefaction array 18800 reads.rds")

## try median first
BACT_RARE_MEDIAN <- apply(BACT_RARE_ARR, c(1,2), median)
dim(BACT_RARE_MEDIAN)
rowSums(BACT_RARE_MEDIAN) # row sums of median are not consistent


BACT_RARE_MEAN <- apply(BACT_RARE_ARR, c(1,2), mean)
dim(BACT_RARE_MEAN)
rowSums(BACT_RARE_MEAN)

## convert to proportion
BACT_RARE_MEAN <- BACT_RARE_MEAN/18800
rowSums(BACT_RARE_MEAN)

## Dq calculations
q <- -5:5
BACT_Dq <- matrix(nrow = 436, ncol = length(q))

for (i in 1:ncol(BACT_Dq)){
  for (j in 1:nrow(BACT_Dq)){
    
    Temp <- BACT_RARE_MEAN[j,] # take obs per site
    Temp <- Temp[Temp>0] # remove 0 values
    ifelse(q[i] == 1,
           BACT_Dq[j,i] <- exp(-sum(Temp*log(Temp))),
           BACT_Dq[j,i] <- 1/sum(Temp^q[i])^(1/(q[i] - 1))) # Calculate Hill's number
    
    
  }
  
}

head(BACT_Dq)
summary(BACT_Dq)

colnames(BACT_Dq) <- paste("BACTq",-5:5, sep="_")
BACT_Dq <- as.data.frame(BACT_Dq)
BACT_Dq$REP_ID <- str_sub(rownames(BACT_ORD), start=2)

summary(BACT_Dq)
multi.hist(BACT_Dq[,1:11])
pairs.panels(BACT_Dq[,1:11], method="spearman", rug=FALSE)

## Fungi ####
## data prep

FUNGI <- rbind(BLAST_m, BLAST_um)
colnames(FUNGI)

FUNGI_ORD <- FUNGI

rownames(FUNGI_ORD) <- FUNGI_ORD$OTU_ID
rownames(FUNGI_ORD)[1:10]

## limit to plot data
FUNGI_ORD <- FUNGI_ORD[,2:438]
colnames(FUNGI_ORD)

## remove OTUs with no obs and transpose matrix
FUNGI_ORD <- FUNGI_ORD[rowSums(FUNGI_ORD)>0,]
FUNGI_ORD <- t(FUNGI_ORD)

## rarefaction to create random community assembly
set.seed(260618)
sort(rowSums(FUNGI_ORD))
## try rarefying to 1750 reads
dim(FUNGI_ORD)
FUNG_RARE_ARR <- array(dim=c(437,8407,10))

for (i in 1:10){
  FUNG_RARE_ARR[,,i] <- rrarefy(FUNGI_ORD, 1750)
}

FUNG_RARE_MEAN <- apply(FUNG_RARE_ARR, c(1,2), mean)
dim(FUNG_RARE_MEAN)
rowSums(FUNG_RARE_MEAN)

## convert to proportion
FUNG_RARE_MEAN <- FUNG_RARE_MEAN/1750
rowSums(FUNG_RARE_MEAN)

## Dq calculations
q <- -5:5
FUNG_Dq <- matrix(nrow = 437, ncol = length(q))

for (i in 1:ncol(FUNG_Dq)){
  for (j in 1:nrow(FUNG_Dq)){
    
    Temp <- FUNG_RARE_MEAN[j,] # take obs per site
    Temp <- Temp[Temp>0] # remove 0 values
    ifelse(q[i] == 1,
           FUNG_Dq[j,i] <- exp(-sum(Temp*log(Temp))),
           FUNG_Dq[j,i] <- 1/sum(Temp^q[i])^(1/(q[i] - 1))) # Calculate Hill's number
    
    
  }
  
}

head(FUNG_Dq)
summary(FUNG_Dq)

colnames(FUNG_Dq) <- paste("FUNGq",-5:5, sep="_")
FUNG_Dq <- as.data.frame(FUNG_Dq)
FUNG_Dq$REP_ID <- str_sub(rownames(FUNGI_ORD), start=2)

summary(FUNG_Dq)
multi.hist(FUNG_Dq[,1:11])
pairs.panels(FUNG_Dq[,1:11], method="spearman", rug=FALSE)

## Compare to Dq from texture
# This data from Multifractals v3.R
summary(Res_fil_wide)
Dq_mer <- Res_fil_wide
colnames(Dq_mer) <- c("REP_ID", paste("TEXT",-5:5, sep="_"))

Dq_mer <- merge(BACT_Dq, Dq_mer, by="REP_ID", all=TRUE) %>% merge(FUNG_Dq, by="REP_ID", all=TRUE)
summary(Dq_mer)
multi.hist(Dq_mer[,12:23])
par(mfrow=c(1,1))

Dq_cor <- cor(na.omit(Dq_mer[,2:34]))
Dq_cor

qgraph(Dq_cor, layout = "spring")

Dq_cor <- cor(na.omit(Dq_mer[,2:34]), method="spearman")
Dq_cor

qgraph(Dq_cor, layout = "spring", posCol = "#56B4E9", negCol = "#D55E00",
       color = c(rep("#FFB5C5",11), rep("#EEC900",11), rep("aquamarine2",11)),
              filetype = "png", filename = "Dq texture bacteria fungi",
              width=8, height=6)

Dq_mer_SSC <- merge(Dq_mer, Res_fil_wide_cor_dt[,1:120], by="REP_ID") 

DqSSC_cor <- cor(na.omit(Dq_mer_SSC[,2:153]), method="spearman")
DqSSC_cor

palette <- colorRampPalette(colors = c("firebrick1","lightgoldenrod1","dodgerblue1"))
col <- c(rep("darkorchid1",11), rep("hotpink",11), rep("springgreen2",11),
         palette(116),"firebrick1", "lightgoldenrod1","dodgerblue1")
qgraph(DqSSC_cor, layout = "spring", posCol = "#56B4E9", negCol = "#D55E00",
       color = col, filetype = "png", filename = "Dq texture bacteria fungi SSC",
       width=8, height=6)
