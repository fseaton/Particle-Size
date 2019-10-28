## Doing NMDS for bacterial data
# mainly for interest
library(tidyverse)
library(vegan)
library(qgraph)
library(psych)
library(magrittr)
library(vegan3d)
library(rgl)

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
BACT_RARE_MEAN2 <- BACT_RARE_MEAN[rownames(BACT_RARE_MEAN) %in% rownames(PSD),]
rownames(PSD)
PSD <- PSD[rownames(BACT_RARE_MEAN2),]

COR_BACT_MER <- cbind(BACT_RARE_MEAN2, PSD)
# 56849 variables
# need to appear in 50% of sites
COR_BACT_MER <- na.omit(COR_BACT_MER) # remove 30 sites
COR_BACT_MER <- COR_BACT_MER[,colSums(COR_BACT_MER>0)>(0.5*329)]
dim(COR_BACT_MER) #329 4393

cor_bact <- Hmisc::rcorr(as.matrix(COR_BACT_MER), type = "spearman")

summary(cor_bact$r[1:4279,4280:4393])

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

apply(cor2[,4280:4393], 2, function(x) length(x[x > 0]))
# [1]  31  31  28  26  26  31  32  28  29  30  32  46 215 593 832 891 824
# [18] 766 708 639 613 585 569 559 551 553 548 557 571 573 583 607 649 636
# [35] 588 547 546 540 504 473 442 380 329 286 251 232 218 209 202 204 193
# [52] 185 174 162 155 150 145 138 128 119 104  86  76  61  57  55  51  43
# [69]  36  36  36  36  36  36  35  35  36  36  36  35  32  37  37  38  36
# [86]  32  30  29  29  29  29  28  29  29  29  29  29  31  33  33  33  34
# [103]  33  31  26  21  19  17  16  15  16  18  19  20

cor_check <- cor2
dimnames(cor_check) <- list(colnames(COR_BACT_MER),colnames(COR_BACT_MER))

bact_w_match <- rowSums(cor_check[,4280:4393]) > 0

cor_check <- cor_check[bact_w_match, bact_w_match]
dim(cor_check) #[1] 1188 1188

colnames(cor_check)[1107:1188]

# remove all bacteria to bacteria correlations and particle-size to particle
# size correlations
cor_check[1:1106,1:1106] <- 0
cor_check[1107:1188,1107:1188] <- 0

palette <- colorRampPalette(colors = c("firebrick1","lightgoldenrod1","dodgerblue1"))
cols <- c(rep("black", 1106),palette(116)[c(1:79,108:110)])
qgraph::qgraph(cor_check, layout = "spring", color = cols)

# positive only
corsel <- cor_check > 0
cor3 <- matrix(ncol = ncol(corsel), nrow = nrow(corsel))
for(i in 1:ncol(corsel)){
  for(j in 1:ncol(corsel)){
    cor3[i,j] <- ifelse(corsel[i,j], cor_check[i,j], 0)
  }
}
colnames(cor3) <- colnames(cor_check)
colSums(cor3[1:1106,1107:1188]>0)
x <- colSums(cor3[1:1106,1107:1188]>0)
dat <- data.frame(Sizes = as.numeric(substring(names(x),2)), Correlations = x)
ggplot(dat, aes(x = Sizes, y = Correlations)) +  scale_x_log10(limits = c(0.04,2000)) + 
  geom_vline(xintercept  = 2.2, linetype = "dashed") + geom_vline(xintercept = 63, linetype = "dashed")+
  geom_line() +
  theme_bw() + theme(panel.grid.minor = element_blank()) 
ggsave("Number of bacterial OTUs correlated with each particle size bin.png", 
       path = "./Graphs/", width = 10, height = 10, units = "cm")

qgraph::qgraph(cor3, layout = "spring", color = cols, labels = FALSE,
               posCol = "grey")

# get names of bacteria
bact_otu <- colnames(cor_check)[1:1106]
bact_otu_id <- BACT[BACT$OTU_ID %in% bact_otu, c(1,445:452)]
table(bact_otu_id$Phylum)

bact_otu_id <- droplevels(bact_otu_id)
length(unique(bact_otu_id$taxonomy)) #219

bact_tab <- count(bact_otu_id, taxonomy)

write.csv(bact_tab, "../bact_tab.csv", row.names = FALSE)

# now plot correlation network with bacterial nodes coloured by class
rownames(bact_otu_id) <- bact_otu_id$OTU_ID
identical(rownames(bact_otu_id), bact_otu)
bact_otu_id <- bact_otu_id[bact_otu,]

bact_phylum_col <- substring(bact_otu_id$Phylum, 5)
summary(as.factor(bact_phylum_col))

bact_phylum_col <- recode(bact_phylum_col, 
                         "Actinobacteria" = "#E69F00",
                         "Planctomycetes" = "#56B4E9",
                         "Proteobacteria" = "#009E73",
                         "Chloroflexi" = "#F0E442",
                         "Verrucomicrobia" = "#0072B2",
                         "Bacteroidetes" = "#D55E00",
                         "Acidobacteria" = "#CC79A7")

bact_phylum_col[!grepl("#",bact_phylum_col)] <- "#000000"
summary(as.factor(bact_phylum_col))

cols <- c(bact_phylum_col,gray.colors(n=116)[c(1:79,108:110)])
shp <- c(rep("circle",1106),rep("square",82))
siz <- c(rep(1,1106),rep(0.5,82))
qgraph::qgraph(cor3, layout = "spring", color = cols, labels = FALSE,
               posCol = "grey", shape = shp)


# get proportion of phyla across whole dataset
BACT_TAX <- select(BACT, OTU_ID, Phylum:Genus) %>% column_to_rownames("OTU_ID")
BACT_TAX <- BACT_TAX[colnames(BACT_RARE_MEAN2),]
dim(BACT_TAX)
BACT_TAX <- cbind(BACT_TAX, t(BACT_RARE_MEAN2))
str(BACT_TAX)
overall_bact_phyla <- count(BACT_TAX[rowSums(BACT_TAX[,6:ncol(BACT_TAX)]>0)>(0.5*329),], Phylum)
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


# Procrustes test ####

# calculate distance objects using bray for bacteria and euclidean for texture
# and pH
psd_dist <- dist(PSD)
bact_bray <- vegdist(BACT_RARE_MEAN2)
ph <- Soil %>% filter(REP_ID %in% rownames(BACT_RARE_MEAN2)) %>%
  column_to_rownames("REP_ID") %>% select(C_B_PH_CACL2)
ph <- ph[rownames(BACT_RARE_MEAN2),]
ph_dist <- dist(ph)

# run procrustes analysis on bacteria directly
btxpr <- protest(psd_dist, bact_bray)
btxpr
# Call:
#   protest(X = psd_dist, Y = bact_bray) 
# 
# Procrustes Sum of Squares (m12 squared):        0.9096 
# Correlation in a symmetric Procrustes rotation: 0.3007 
# Significance:  0.001 
# 
# Permutation: free
# Number of permutations: 999

bphpr <- protest(ph_dist, bact_bray)
bphpr
# Call:
#   protest(X = ph_dist, Y = bact_bray) 
# 
# Procrustes Sum of Squares (m12 squared):        0.5401 
# Correlation in a symmetric Procrustes rotation: 0.6782 
# Significance:  0.001 
# 
# Permutation: free
# Number of permutations: 999

txphpr <- protest(psd_dist, ph_dist)
txphpr

# trying procrustes on residuals on the pH model. This seems wrong.
bphpr_res <- residuals(bphpr)
btxphpr <- protest(psd_dist, dist(bphpr_res))
btxphpr
# Call:
#   protest(X = psd_dist, Y = dist(bphpr_res)) 
# 
# Procrustes Sum of Squares (m12 squared):        0.9733 
# Correlation in a symmetric Procrustes rotation: 0.1634 
# Significance:  0.006 
# 
# Permutation: free
# Number of permutations: 999

# now doing procrustes on the residuals of a RDA upon pH as I think that works
# better
bphrda <- capscale(bact_bray ~ ph)
bphrda
anova(bphrda)
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = bact_bray ~ ph)
#             Df SumOfSqs      F Pr(>F)    
#   Model      1   12.754 65.844  0.001 ***
#   Residual 327   63.342               

bphrda_res <- residuals(bphrda)

bphtxrda <- protest(psd_dist, bphrda_res)
bphtxrda
# Call:
#   protest(X = psd_dist, Y = bphrda_res) 
# 
# Procrustes Sum of Squares (m12 squared):        0.8964 
# Correlation in a symmetric Procrustes rotation: 0.3219 
# Significance:  0.001 
# 
# Permutation: free
# Number of permutations: 999

PSD2 <- matrix(nrow = nrow(PSD), ncol = 29)
for (i in 1:29){
  PSD2[,i]<- rowSums(PSD[,seq(4*i-3,4*i)])
}
rownames(PSD2) <- rownames(PSD)
colnames(PSD2) <- substring(colnames(PSD)[seq(1,113,4)], 1, 5)
PSD2 <- cbind(PSD2, ph) %>% as.data.frame()
colnames(PSD2)


b_tx_ph <- dbrda(BACT_RARE_MEAN2 ~ X0.04 + X0.06 + X0.09 + X0.13 + X0.19 +
                      X0.28 + X0.41 + X0.59 + X0.86 + X1.26 + X1.83 + X2.65 +
                      X3.86 + X5.60 + X8.14 + X11.8 + X17.1 + X24.9 + X36.2 +
                      X52.6 + X76.4 + X110. + X161. + X234. + X339. + X493. +
                      X716. + X1041 + X1511 + Condition(ph), data = PSD2,
                    distance = "bray")
b_tx_ph
par(mfrow = c(1,2))
plot(b_tx_ph, choices = c(1,2))
plot(b_tx_ph, choices = c(3,4))
par(mfrow = c(1,1))
anova(b_tx_ph)
anova(b_tx_ph, by = "terms")

## nmds ####
bact.nmds <- metaMDS(BACT_RARE_MEAN2, trymax = 1000, k = 3)
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
# Stress:     0.0655898 
# Stress type 1, weak ties
# Two convergent solutions found after 43 tries
# Scaling: centring, PC rotation, halfchange scaling 
# Species: expanded scores based on ‘wisconsin(sqrt(BACT_RARE_MEAN))’ 

plot(bact.nmds, display = "sites")
plot(bact.nmds, display = "sites", choices = c(1,3))

D1 <- Res_fil_mer_wideD$qD1
names(D1) <- Res_fil_mer_wideD$REP_ID
D1 <- D1[rownames(BACT_RARE_MEAN2)]

PSD2 <- matrix(nrow = nrow(PSD), ncol = 9)
for (i in 1:8){
  PSD2[,i]<- rowSums(PSD[,seq(13*i-12,13*i)])
}
PSD2[,9] <- rowSums(PSD[,105:116])
rownames(PSD2) <- rownames(PSD)
colnames(PSD2) <- substring(colnames(PSD)[seq(7,112,13)], 2, 5)
PSD2 <- as.data.frame(PSD2) %>% rownames_to_column("REP_ID")
colnames(PSD2)
colnames(PSD2) <- c("REP_ID","FC","MC","CC","FSi","MSi","CSi","FSa","MSa","CSa")


bact_env_dat <- Res_fil_mer_wideD %>% 
  select(REP_ID, CARB = C_FE_CTOTAL, pH = C_B_PH_CACL2, qD0:qD2) %>%
  full_join(PSD2) %>%
  column_to_rownames("REP_ID")
bact_env_dat <- bact_env_dat[rownames(BACT_RARE_MEAN2),]
summary(bact_env_dat)

# envfit
bact.envfit23 <- envfit(bact.nmds,bact_env_dat, na.rm = TRUE, choices = 1:3)

par(mfrow=c(2,1), mar = c(4,4.5,1,2))
plot(bact.nmds, display = "sites", cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(bact.envfit23, cex = 1.5)
plot(bact.nmds, display = "sites", choices = c(1,3), cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(bact.envfit23, choices = c(1,3), cex = 1.5)
par(mfrow=c(1,1))

bact.envfit23
# NMDS1    NMDS2    NMDS3     r2 Pr(>r)    
# CARB  0.38118 -0.12141 -0.91649 0.4024  0.001 ***
#   pH   -0.93116 -0.24406 -0.27088 0.7090  0.001 ***
#   qD0   0.21117 -0.76843 -0.60409 0.0256  0.035 *  
#   qD1   0.53333 -0.16516 -0.82963 0.0085  0.424    
# qD2   0.10293  0.82294 -0.55872 0.0083  0.448    
# FC   -0.46268  0.81453  0.34995 0.0415  0.006 ** 
#   MC   -0.32965  0.35665  0.87414 0.2221  0.001 ***
#   CC   -0.31763  0.39478  0.86213 0.2390  0.001 ***
#   FSi  -0.17780  0.69945  0.69221 0.1710  0.001 ***
#   MSi  -0.05766  0.92923  0.36498 0.1568  0.001 ***
#   CSi  -0.10399  0.88646 -0.45098 0.0274  0.034 *  
#   FSa   0.14527 -0.24534 -0.95849 0.0882  0.001 ***
#   MSa   0.16603 -0.70895 -0.68544 0.1671  0.001 ***
#   CSa   0.25575 -0.96570 -0.04485 0.1131  0.001 ***
  ---

ordirgl(bact.nmds, display = "sites", choices=1:3, envfit=bact.envfit23, 
        ax.col = "black", arr.col = "darkorchid", col = "gray")

rgl.bg(color=c("white", "black"))
ordirgl(bact.nmds, display = "sites", choices=1:3, envfit=bact.envfit23, 
        ax.col = "black", arr.col = "darkorchid", col = "gray")

writeWebGL(dir = "Graphs/webGL/", filename = "Graphs/webGL/index.html")


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
FUNG_RARE_MEAN2 <- FUNG_RARE_MEAN[rownames(FUNG_RARE_MEAN) %in% rownames(PSD),]
rownames(PSD)
F_PSD <- PSD[rownames(FUNG_RARE_MEAN2),]

COR_FUNG_MER <- cbind(FUNG_RARE_MEAN2, F_PSD)
# 8274 variables
# need to appear in 50% of sites
COR_FUNG_MER <- na.omit(COR_FUNG_MER) # remove 30 sites
COR_FUNG_MER <- COR_FUNG_MER[,colSums(COR_FUNG_MER>0)>(0.25*329)]
dim(COR_FUNG_MER) #329 289

cor_fung <- Hmisc::rcorr(as.matrix(COR_FUNG_MER), type = "spearman")

colnames(COR_FUNG_MER)
summary(cor_fung$r[1:175,176:289])

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

apply(cor2[,176:289], 2, function(x) length(x[x > 0]))
# [1]  34  35  34  34  34  34  35  34  34  33  34  42  74 101  93  89  83
# [18]  80  72  72  69  68  68  68  67  67  67  67  67  67  70  72  72  76
# [35]  76  75  73  72  70  68  68  67  63  61  59  58  56  57  57  57  58
# [52]  59  59  60  60  60  61  62  64  63  65  63  64  66  62  63  61  47
# [69]  41  38  37  38  37  37  37  38  38  38  37  37  41  41  41  41  39
# [86]  38  32  29  29  29  29  30  31  32  32  33  32  33  33  33  33  32
# [103]  31  30  28  23  20  18  17  16  18  20  21  22

cor_check <- cor2
dimnames(cor_check) <- list(colnames(COR_FUNG_MER),colnames(COR_FUNG_MER))

fung_w_match <- rowSums(cor_check[,176:289]) > 0

cor_check <- cor_check[fung_w_match, fung_w_match]
dim(cor_check) #[1] 133 133

colnames(cor_check)[54:133]

# remove all fungi to fungi correlations and particle-size to particle
# size correlations
cor_check[1:53,1:53] <- 0
cor_check[54:133,54:133] <- 0

palette <- colorRampPalette(colors = c("firebrick1","lightgoldenrod1","dodgerblue1"))
cols <- c(rep("black", 53),palette(116)[1:80])
qgraph::qgraph(cor_check, layout = "spring", color = cols)

# positive only
corsel <- cor_check > 0
cor3 <- matrix(ncol = ncol(corsel), nrow = nrow(corsel))
for(i in 1:ncol(corsel)){
  for(j in 1:ncol(corsel)){
    cor3[i,j] <- ifelse(corsel[i,j], cor_check[i,j], 0)
  }
}
colnames(cor3) <- colnames(cor_check)
colSums(cor3[1:53,54:133]>0)
x <- colSums(cor3[1:53,54:133]>0)
dat <- data.frame(Sizes = as.numeric(substring(names(x),2)), Correlations = x)
ggplot(dat, aes(x = Sizes, y = Correlations)) +  scale_x_log10(limits = c(0.04,2000)) + 
  geom_vline(xintercept  = 2.2, linetype = "dashed") + geom_vline(xintercept = 63, linetype = "dashed")+
  geom_line() +
  theme_bw() + theme(panel.grid.minor = element_blank())
ggsave("Number of fungal OTUs correlated with each particle size bin.png", 
       path = "./Graphs/", width = 10, height = 10, units = "cm")


qgraph::qgraph(cor3, layout = "spring", color = cols, labels = FALSE,
               posCol = "grey")

# get names of fungi
fung_otu <- colnames(cor_check)[1:53]
fung_otu_id <- FUNGI[FUNGI$OTU_ID %in% fung_otu, c(1,444:453)]
table(fung_otu_id$Taxon)

fung_otu_id <- droplevels(fung_otu_id)
length(unique(fung_otu_id$taxonomy)) #31

fung_tab <- count(fung_otu_id, taxonomy)

write.csv(fung_tab, "../fung_tab.csv", row.names = FALSE)

# now plot correlation network with fungal nodes coloured by class
rownames(fung_otu_id) <- fung_otu_id$OTU_ID
identical(rownames(fung_otu_id), fung_otu)

class <- fung_otu_id %>% separate(taxonomy, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), 
                         sep ="; ", remove = FALSE) %>%
  select(Class)
fung_class_col <- substring(class$Class, 4)
fung_class_col[is.na(fung_class_col)] <- "unidentified"
fung_class_col

fung_class_col <- recode(fung_class_col, "Sordariomycetes" = "#000000",
       "Mucoromycetes" = "#E69F00",
       "Leotiomycetes" = "#56B4E9",
       "Eurotiomycetes" = "#009E73",
       "Mortierellomycetes" = "#F0E442",
       "Agaricomycetes" = "#0072B2",
       "Dothideomycetes" = "#D55E00",
       "Tremellomycetes" = "#CC79A7",
       "unidentified" = "#999999")

cols <- c(fung_class_col,gray.colors(n=116)[1:80])
shp <- c(rep("circle",53),rep("square",80))
siz <- c(rep(2,53),rep(1,80))
qgraph::qgraph(cor3, layout = "spring", color = cols, labels = FALSE,
               posCol = "grey", shape = shp, vsize = siz)


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


# Procrustes test ####
psd_dist <- dist(PSD)
fung_bray <- vegdist(FUNG_RARE_MEAN2)
ph <- Soil %>% filter(REP_ID %in% rownames(FUNG_RARE_MEAN2)) %>%
  column_to_rownames("REP_ID") %>% select(C_B_PH_CACL2)
ph <- ph[rownames(FUNG_RARE_MEAN2),]
ph_dist <- dist(ph)

ftxpr <- protest(psd_dist, fung_bray)
ftxpr
# Call:
#   protest(X = psd_dist, Y = fung_bray) 
# 
# Procrustes Sum of Squares (m12 squared):        0.9677 
# Correlation in a symmetric Procrustes rotation: 0.1798 
# Significance:  0.001 
# 
# Permutation: free
# Number of permutations: 999

fphpr <- protest(ph_dist, fung_bray)
fphpr
# Call:
#   protest(X = ph_dist, Y = fung_bray) 
# 
# Procrustes Sum of Squares (m12 squared):        0.8262 
# Correlation in a symmetric Procrustes rotation: 0.4169 
# Significance:  0.001 
# 
# Permutation: free
# Number of permutations: 999

fphpr_res <- residuals(fphpr)
ftxphpr <- protest(psd_dist, dist(fphpr_res))
ftxphpr
# Call:
#   protest(X = psd_dist, Y = dist(fphpr_res)) 
# 
# Procrustes Sum of Squares (m12 squared):        0.9752 
# Correlation in a symmetric Procrustes rotation: 0.1576 
# Significance:  0.009 
# 
# Permutation: free
# Number of permutations: 999

# now doing procrustes on the residuals of a RDA upon pH as I think that works
# better
fphrda <- capscale(fung_bray ~ ph)
fphrda
anova(fphrda)
# Model: capscale(formula = fung_bray ~ ph)
#            Df SumOfSqs      F Pr(>F)    
# Model      1    5.603 14.272  0.001 ***
# Residual 327  128.374                   

fphrda_res <- residuals(fphrda)

fphtxrda <- protest(psd_dist, fphrda_res)
fphtxrda
# Call:
#   protest(X = psd_dist, Y = fphrda_res) 
# 
# Procrustes Sum of Squares (m12 squared):        0.9722 
# Correlation in a symmetric Procrustes rotation: 0.1666 
# Significance:  0.005 
# 
# Permutation: free
# Number of permutations: 999

f_tx_ph <- dbrda(FUNG_RARE_MEAN2 ~ X0.04 + X0.06 + X0.09 + X0.13 + X0.19 +
                   X0.28 + X0.41 + X0.59 + X0.86 + X1.26 + X1.83 + X2.65 +
                   X3.86 + X5.60 + X8.14 + X11.8 + X17.1 + X24.9 + X36.2 +
                   X52.6 + X76.4 + X110. + X161. + X234. + X339. + X493. +
                   X716. + X1041 + X1511 + Condition(ph), data = PSD2,
                 distance = "bray")
f_tx_ph
par(mfrow=c(1,2))
plot(f_tx_ph, choices = c(1,2))
plot(f_tx_ph, choices = c(3,4))
par(mfrow = c(1,1))
anova(f_tx_ph, by = "terms")
# Permutation test for dbrda under reduced model Terms added sequentially (first
# to last) Permutation: free Number of permutations: 999
#
# Model: dbrda(formula = FUNG_RARE_MEAN2 ~ X0.04 + X0.06 + X0.09 + X0.13 + X0.19
# + X0.28 + X0.41 + X0.59 + X0.86 + X1.26 + X1.83 + X2.65 + X3.86 + X5.60 +
# X8.14 + X11.8 + X17.1 + X24.9 + X36.2 + X52.6 + X76.4 + X110. + X161. + X234.
# + X339. + X493. + X716. + X1041 + X1511 + Condition(ph), data = PSD2, distance
# = "bray")

# Df SumOfSqs      F Pr(>F)    
# X0.04      1    0.646 1.7603  0.009 ** 
#   X0.06      1    1.028 2.8010  0.001 ***
#   X0.09      1    1.055 2.8740  0.001 ***
#   X0.13      1    0.849 2.3127  0.001 ***
#   X0.19      1    0.319 0.8687  0.680    
# X0.28      1    0.543 1.4795  0.043 *  
#   X0.41      1    0.301 0.8200  0.800    
# X0.59      1    0.782 2.1307  0.001 ***
#   X0.86      1    0.379 1.0320  0.362    
# X1.26      1    0.617 1.6800  0.021 *  
#   X1.83      1    0.548 1.4921  0.038 *  
#   X2.65      1    0.494 1.3467  0.087 .  
# X3.86      1    0.554 1.5084  0.039 *  
#   X5.60      1    0.390 1.0631  0.320    
# X8.14      1    0.375 1.0202  0.367    
# X11.8      1    0.452 1.2320  0.132    
# X17.1      1    0.431 1.1741  0.171    
# X24.9      1    0.550 1.4988  0.042 *  
#   X36.2      1    0.491 1.3387  0.079 .  
# X52.6      1    0.298 0.8109  0.819    
# X76.4      1    1.077 2.9323  0.001 ***
#   X110.      1    0.454 1.2358  0.106    
# X161.      1    0.560 1.5265  0.024 *  
#   X234.      1    0.340 0.9248  0.590    
# X339.      1    0.346 0.9433  0.544    
# X493.      1    0.383 1.0438  0.357    
# X716.      1    0.373 1.0147  0.389    
# X1041      1    0.474 1.2915  0.116    
# X1511      1    0.322 0.8770  0.709    

## nmds ####
fung.nmds <- metaMDS(FUNG_RARE_MEAN2, trymax = 1000, k = 2)
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
# Stress:     0.1607794 
# Stress type 1, weak ties
# Two convergent solutions found after 250 tries
# Scaling: centring, PC rotation, halfchange scaling 
# Species: expanded scores based on ‘wisconsin(sqrt(FUNG_RARE_MEAN))’ 

plot(fung.nmds, display = "sites")

D1 <- Res_fil_mer_wideD$qD1
names(D1) <- Res_fil_mer_wideD$REP_ID
D1 <- D1[rownames(FUNG_RARE_MEAN2)]

fung.D1.ordi <- ordisurf(fung.nmds, D1) #linear 

fung_env_dat <- Res_fil_mer_wideD %>% 
  select(REP_ID, CTOT = C_FE_CTOTAL, PH = C_B_PH_CACL2, qD0:qD2) %>%
  full_join(PSD2) %>% column_to_rownames("REP_ID")
fung_env_dat <- fung_env_dat[rownames(FUNG_RARE_MEAN2),]
summary(fung_env_dat)


fung.envfit <- envfit(fung.nmds,fung_env_dat, na.rm = TRUE)

plot(fung.nmds, display = "sites", cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(fung.envfit, cex = 1.5)

fung.envfit
# NMDS1    NMDS2     r2 Pr(>r)    
# CTOT  0.88170  0.47182 0.3107  0.001 ***
#   PH   -0.92110 -0.38932 0.6232  0.001 ***
#   qD0   0.48294 -0.87565 0.0193  0.049 *  
#   qD1   0.89961  0.43669 0.0086  0.248    
# qD2   0.12283  0.99243 0.0111  0.159    
# FC   -0.91658  0.39985 0.0322  0.006 ** 
#   MC   -0.98231  0.18729 0.1585  0.001 ***
#   CC   -0.95320  0.30236 0.1738  0.001 ***
#   FSi  -0.58746  0.80926 0.1181  0.001 ***
#   MSi  -0.24709  0.96899 0.1170  0.001 ***
#   CSi  -0.27227  0.96222 0.0179  0.056 .  
# FSa   0.92101 -0.38955 0.0311  0.007 ** 
#   MSa   0.59634 -0.80273 0.1082  0.001 ***
#   CSa   0.52718 -0.84975 0.0965  0.001 ***