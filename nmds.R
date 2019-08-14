## Doing NMDS for bacterial data
# mainly for interest
library(tidyverse)
library(vegan)
library(qgraph)
library(psych)
library(magrittr)

#### bacteria data prep ####
## get data 
dim(BACT)
colnames(BACT)


## quick check to see if there are weird dominant taxa
summary(BACT$taxonomy)

BACT <- separate(BACT, col=taxonomy, sep=";", 
                 into = c("Kingdom", "Phylum", "Class", "Order", 
                          "Family", "Genus", "Species"),
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


# limit to sites within particle size data
BACT_ORD <- BACT_ORD[,colnames(BACT_ORD) %in% paste0("X",GMEP$REP_ID)]
# 361 variables

## get rid of 'empty' OTUs and transpose data for vegan
BACT_ORD <- BACT_ORD[rowSums(BACT_ORD)>0,]
BACT_ORD <- t(BACT_ORD)



## removing reps with read depth less than 2000 (there's two, and the next
## smallest is 18873)
sort(rowSums(BACT_ORD), decreasing = T)
dim(BACT_ORD) #361 57477
BACT_ORD <- BACT_ORD[rowSums(BACT_ORD)>2000,]
dim(BACT_ORD) # 359 57477



## bacteria ####
sort(rowSums(BACT_ORD))
## try rarefying to 18800 reads
dim(BACT_ORD)
min(colSums(BACT_ORD)) #1
BACT_ORD <- BACT_ORD[,colSums(BACT_ORD)>1]
dim(BACT_ORD) #[1]   359 56733

## rarefaction
BACT_RARE_ARR <- array(dim=c(359,56733,10))
set.seed(010819)

for (i in 1:10){
  BACT_RARE_ARR[,,i] <- rrarefy(BACT_ORD, 18800)
}

BACT_RARE_MEAN <- apply(BACT_RARE_ARR, c(1,2), mean)
dim(BACT_RARE_MEAN)
rowSums(BACT_RARE_MEAN)

# correlation network ####
dimnames(BACT_RARE_MEAN) <- dimnames(BACT_ORD)
rownames(BACT_RARE_MEAN) <- substring(rownames(BACT_RARE_MEAN), 2)
rownames(PSD)
PSD <- PSD[rownames(BACT_RARE_MEAN),]

COR_BACT_MER <- cbind(BACT_RARE_MEAN, PSD)
# 56849 variables
# need to appear in 50% of sites
COR_BACT_MER <- na.omit(COR_BACT_MER) # remove 30 sites
COR_BACT_MER <- COR_BACT_MER[,colSums(COR_BACT_MER>0)>(0.5*359)]
dim(COR_BACT_MER) #329 3778

cor_bact <- Hmisc::rcorr(as.matrix(COR_BACT_MER), type = "spearman")

cor_bact2 <- cor(COR_BACT_MER, method = "spearman")


summary(cor_bact$r[1:3664,3665:3778])

# limit to correlations with p value below bonferroni point
n <- ncol(COR_BACT_MER)
corsel <- cor_bact$P < 0.05/((n*(n+1))/2)
diag(corsel) <- FALSE

cor2 <- matrix(ncol = n, nrow = n)
for(i in 1:n){
  for(j in 1:n){
    cor2[i,j] <- ifelse(corsel[i,j], cor_bact$r[i,j], 0)
  }
}

apply(cor2[,3665:3778], 2, function(x) length(x[x > 0]))
# [1]  31  31  28  26  29  31  32  28  29  30  32  51 212 579 809 878 806 759
# 689 639 610 594 574 562 556 557 550 561 566 575 584 613 645 628 579 536 533
# 532 500 471 439 375 323 282 252 238 222 214 203 199 188 183 179 166 164 156
# 148 140 123 111 102  88  75  63  58  55  52  43  36  36  36  36  36 36  35  35
# 36  36  35  34  31  35  37  38  34  31  30 29  29  30  30 28  28  29  29  29
# 30  32  32  32  33  34  32  31  28  22  19  17  16  16  16 18  19  20

cor_check <- cor2
dimnames(cor_check) <- list(colnames(COR_BACT_MER),colnames(COR_BACT_MER))

bact_wo_match <- rowSums(cor_check[,3665:3778]) > 0

cor_check <- cor_check[bact_wo_match, bact_wo_match]
dim(cor_check) #[1] 1147 1147

colnames(cor_check)[1066:1147]

# remove all bacteria to bacteria correlations and particle-size to particle
# size correlations
cor_check[1:1065,1:1065] <- 0
cor_check[1066:1147,1066:1147] <- 0

palette <- colorRampPalette(colors = c("firebrick1","lightgoldenrod1","dodgerblue1"))
cols <- c(rep("black", 1065),palette(82))
qgraph::qgraph(cor_check, layout = "spring", color = cols)

# positive only
corsel <- cor_check > 0
cor3 <- matrix(ncol = ncol(corsel), nrow = nrow(corsel))
for(i in 1:ncol(corsel)){
  for(j in 1:ncol(corsel)){
    cor3[i,j] <- ifelse(corsel[i,j], cor_check[i,j], 0)
  }
}
qgraph::qgraph(cor3, layout = "spring", color = cols, labels = FALSE)

# get names of bacteria
bact_otu <- colnames(cor_check)[1:1065]
bact_otu_id <- BACT[BACT$OTU_ID %in% bact_otu, 445:452]
table(bact_otu_id$Phylum)

bact_otu_id <- droplevels(bact_otu_id)
length(unique(bact_otu_id$taxonomy)) #219

bact_tab <- count(bact_otu_id, taxonomy)

write.csv(bact_tab, "../bact_tab.csv", row.names = FALSE)

# get proportion of phyla across whole dataset
overall_bact_phyla <- count(BACT[rowSums(BACT[,2:439]>0)>218,], Phylum)
cor_bact_phyla <- count(bact_otu_id, Phylum)
bact_phyla <- full_join(overall_bact_phyla, cor_bact_phyla, by = "Phylum") %>%
  filter(n.x > 10) %>%
  mutate(Phylum = substring(Phylum, 5)) %>%
  tidyr::gather(key = "key", value = "value", n.x,n.y) %>%
  mutate(key = recode_factor(key, n.x = "Whole dataset", 
                             n.y = "Correlated with particle size bin"))

ggplot(bact_phyla, aes(x = reorder(Phylum, -value), y = value)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~key, ncol = 1, scales = "free_y") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Count", x = "Phylum")
ggsave("Bact phyla comparison overall and psd correlated.png", path = "./Graphs/",
       width = 12, height = 16, units = "cm")

## nmds ####
bact.nmds <- metaMDS(BACT_RARE_MEAN, trymax = 1000, k = 3)
bact.nmds
# Call:
#   metaMDS(comm = BACT_RARE_MEAN, k = 3, trymax = 1000) 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     wisconsin(sqrt(BACT_RARE_MEAN)) 
# Distance: bray 
# 
# Dimensions: 3 
# Stress:     0.06625096 
# Stress type 1, weak ties
# Two convergent solutions found after 61 tries
# Scaling: centring, PC rotation, halfchange scaling 
# Species: expanded scores based on ‘wisconsin(sqrt(BACT_RARE_MEAN))’ 

plot(bact.nmds, display = "sites")
plot(bact.nmds, display = "sites", choices = c(1,3))

D1 <- Res_fil_mer_wideD$qD1
names(D1) <- paste0("X", Res_fil_mer_wideD$REP_ID)
D1 <- D1[rownames(BACT_ORD)]

bact.D1.ordi <- ordisurf(bact.nmds, D1) #linear 
bact_env_dat <- Res_fil_mer_wideD %>% 
  select(REP_ID, CTOT = C_FE_CTOTAL, pH = C_B_PH_CACL2, qD0:D2_D1) %>%
  left_join(select(GMEP, REP_ID, CLAY:SAND)) %>%
  set_rownames(Res_fil_mer_wideD$REP_ID) %>% select(-REP_ID)
bact_env_dat <- bact_env_dat[substring(rownames(BACT_ORD),2),]
summary(bact_env_dat)

# envfit
bact.envfit23 <- envfit(bact.nmds,bact_env_dat, na.rm = TRUE, choices = 1:3)

par(mfrow=c(1,3))
plot(bact.nmds, display = "sites")
plot(bact.envfit23)
plot(bact.nmds, display = "sites", choices = c(1,3))
plot(bact.envfit23, choices = c(1,3))
plot(bact.nmds, display = "sites", choices = 2:3)
plot(bact.envfit23, choices = 2:3)
par(mfrow=c(1,1))


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

# limit to data with textural data
FUNGI_ORD <- FUNGI_ORD[,colnames(FUNGI_ORD) %in% paste0("X",GMEP$REP_ID)]
# 361 variables

## remove OTUs with no obs and transpose matrix
FUNGI_ORD <- FUNGI_ORD[rowSums(FUNGI_ORD)>0,]
FUNGI_ORD <- t(FUNGI_ORD)

## rarefaction to create random community assembly
set.seed(010819)
sort(rowSums(FUNGI_ORD))
## try rarefying to 1750 reads
dim(FUNGI_ORD)
FUNG_RARE_ARR <- array(dim=c(361,8158,10))

for (i in 1:10){
  FUNG_RARE_ARR[,,i] <- rrarefy(FUNGI_ORD, 1750)
}

FUNG_RARE_MEAN <- apply(FUNG_RARE_ARR, c(1,2), mean)
dim(FUNG_RARE_MEAN)
rowSums(FUNG_RARE_MEAN)

# correlation network ####
dimnames(FUNG_RARE_MEAN) <- dimnames(FUNGI_ORD)
rownames(FUNG_RARE_MEAN) <- substring(rownames(FUNG_RARE_MEAN), 2)
rownames(PSD)
F_PSD <- PSD[rownames(FUNG_RARE_MEAN),]

COR_FUNG_MER <- cbind(FUNG_RARE_MEAN, F_PSD)
# 56849 variables
# need to appear in 50% of sites
COR_FUNG_MER <- na.omit(COR_FUNG_MER) # remove 30 sites
COR_FUNG_MER <- COR_FUNG_MER[,colSums(COR_FUNG_MER>0)>(0.25*331)]
dim(COR_FUNG_MER) #331 297

cor_fung <- Hmisc::rcorr(as.matrix(COR_FUNG_MER), type = "spearman")


summary(cor_fung$r[1:184,185:297])

# limit to correlations with p value below bonferroni point
n <- ncol(COR_FUNG_MER)
corsel <- cor_fung$P < 0.05/((n*(n+1))/2)
diag(corsel) <- FALSE

cor2 <- matrix(ncol = n, nrow = n)
for(i in 1:n){
  for(j in 1:n){
    cor2[i,j] <- ifelse(corsel[i,j], cor_fung$r[i,j], 0)
  }
}

apply(cor2[,185:297], 2, function(x) length(x[x > 0]))
# [1]  35  35  34  34  34  35  34  34  33  36  40  70 100  95  90  83  79  76
# 71  68  67  67  67  67  67  67 67  68  69  69  70  72  76  76  75  72
# 72  68  66  65  62  62  61  57  57  57  57  57  57  57  57  59  59  59
# 59  61  61  63  62  63  62  64  66  61  63  59  46  41  38  37  38  37  37  37
# 38  38  38  37 37  41  42  41  41  39  38  33  29  29  29  29  30  32
# 33  32  33  32  33  33  33  33  32  32  30  28 23  20  18  16  16  18
# 20  21  22

cor_check <- cor2
dimnames(cor_check) <- list(colnames(COR_FUNG_MER),colnames(COR_FUNG_MER))

fung_wo_match <- rowSums(cor_check[,185:297]) > 0

cor_check <- cor_check[fung_wo_match, fung_wo_match]
dim(cor_check) #[1] 132 132

colnames(cor_check)[54:132]

# remove all fungi to fungi correlations and particle-size to particle
# size correlations
cor_check[1:53,1:53] <- 0
cor_check[54:132,54:132] <- 0

palette <- colorRampPalette(colors = c("firebrick1","lightgoldenrod1","dodgerblue1"))
cols <- c(rep("black", 53),palette(79))
qgraph::qgraph(cor_check, layout = "spring", color = cols)

# positive only
corsel <- cor_check > 0
cor3 <- matrix(ncol = ncol(corsel), nrow = nrow(corsel))
for(i in 1:ncol(corsel)){
  for(j in 1:ncol(corsel)){
    cor3[i,j] <- ifelse(corsel[i,j], cor_check[i,j], 0)
  }
}
qgraph::qgraph(cor3, layout = "spring", color = cols, labels = FALSE)

# get names of bacteria
fung_otu <- colnames(cor_check)[1:53]
fung_otu_id <- FUNGI[FUNGI$OTU_ID %in% fung_otu, 444:453]
table(fung_otu_id$Taxon)

fung_otu_id <- droplevels(fung_otu_id)
length(unique(fung_otu_id$taxonomy)) #31

fung_tab <- count(fung_otu_id, taxonomy)

write.csv(fung_tab, "../fung_tab.csv", row.names = FALSE)

# get proportion of phyla across whole dataset
overall_fung_phyla <- count(FUNGI[rowSums(FUNGI[,2:439]>0)>110,], taxonomy)
cor_fung_phyla <- count(fung_otu_id, taxonomy)
fung_phyla <- full_join(overall_fung_phyla, cor_fung_phyla, by = "taxonomy") %>%
  filter(n.x > 1) %>%
  # mutate(Phylum = substring(Phylum, 5)) %>%
  tidyr::gather(key = "key", value = "value", n.x,n.y, na.rm = TRUE) %>%
  mutate(key = recode_factor(key, n.x = "Whole dataset", 
                             n.y = "Correlated with particle size bin")) %>%
  separate(taxonomy, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), 
           sep ="; ", remove = FALSE) %>%
  mutate_at(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), 
            function(x) substring(x, 4)) %>%
  group_by(key, Class) %>% summarise(value = sum(value))


ggplot(fung_phyla, aes(x = reorder(Class, value, FUN = function(x){1/(max(x))}), 
                       y = value)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~key, ncol = 1, scales = "free_y") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Count", x = "Class")
ggsave("Fungi class comparison overall and psd correlated.png", path = "./Graphs/",
       width = 12, height = 16, units = "cm")

## nmds ####
fung.nmds <- metaMDS(FUNG_RARE_MEAN, trymax = 1000, k = 2)
fung.nmds
# Call:
#   metaMDS(comm = FUNG_RARE_MEAN, k = 2, trymax = 1000) 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     wisconsin(sqrt(FUNG_RARE_MEAN)) 
# Distance: bray 
# 
# Dimensions: 2 
# Stress:     0.1629087 
# Stress type 1, weak ties
# Two convergent solutions found after 72 tries
# Scaling: centring, PC rotation, halfchange scaling 
# Species: expanded scores based on ‘wisconsin(sqrt(FUNG_RARE_MEAN))’ 

plot(fung.nmds, display = "sites")

D1 <- Res_fil_mer_wideD$qD1
names(D1) <- paste0("X", Res_fil_mer_wideD$REP_ID)
D1 <- D1[rownames(FUNGI_ORD)]

fung.D1.ordi <- ordisurf(fung.nmds, D1) #linear 

fung_env_dat <- Res_fil_mer_wideD %>% 
  select(REP_ID, CTOT = C_FE_CTOTAL, PH = C_B_PH_CACL2, qD0:D2_D1) %>%
  left_join(select(GMEP, REP_ID, CLAY:SAND)) %>%
  set_rownames(Res_fil_mer_wideD$REP_ID) %>% select(-REP_ID)
fung_env_dat <- fung_env_dat[substring(rownames(FUNGI_ORD),2),]
summary(fung_env_dat)


fung.envfit <- envfit(fung.nmds,fung_env_dat, na.rm = TRUE)

plot(fung.nmds, display = "sites")
plot(fung.envfit)
