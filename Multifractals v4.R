# Multifractal parameters calculated for various q
## Based on a combination of Miranda et al Geoderma 134 (2006) 373-385 and code from
## Salat et al Physica A 473 (2017) 467–487
## Using standard grid moment method originating in Chhabra and Jensen Phys. Rev. Lett. 62 (1989) 1327–1330

## The code involves aggregating the particle size data from the 116 bin version into versions with less bins
## So e.g. take the 116 bins and cut in half, there will be two bins with the sum of all the data at < midpoint in one bin
## and > midpoint in the other.
## Here I am repeating this aggregation process 6 times (k=6) with increasing bin width and decreasing total # of bins
## each time
## Once these new datasets have been created information indices are extracted for each REP_ID. These are H0 (the number of
## non-empty bins), H1 (the Shannon index), and H2 (the Simpson index). There will be multiple versions of each index, one
## for each iteration (k), so in total the number of indices = k*num(REP_ID).
## The values for each index are then plotted against log(bin width). The slopes of the linear relationships for H0~log(e), 
## H1~log(e) and H2~ log(e) are equal to D0, D1 and D2 respectively.

## The shape of the particle size distribution below 10um is completely dependent on the optical model used to analyse the
## machine (PIDS) data. We are uncertain which model is the most accurate, or if any truly represent the actual distribution.
## Our choice of optical model is based upon similarity to %Sand, %Silt, %Clay from other methods, which cannot verify
## the distribution. Because of this, analysis of the shape of the relationship (as done in multifractal analysis) seems 
## to have inherent flaws. So we are trying out restricting our analysis to areas of the distribution that are not heavily
## impacted by optical model choice. >10um is known to be unaffected by the model, visual inspection suggests that 1 and 2um
## may also be relatively unaffected by optical model choice.

library(tidyverse)
library(Matrix)
library(gridExtra)

## get data ####
summary(GMEP)
apply(GMEP, 2, anyNA)
GMEP[is.na(GMEP)] <- 0
summary(GMEP)
## select only the particle size bins (along with identifier)
PSD <- filter(GMEP, REPEATED < 1) %>% select(REP_ID, X0.0438996:X2000)
## make identifier the rowname so it can be removed from the actual data
row.names(PSD) <- PSD$REP_ID
PSD <- select(PSD,-REP_ID)
rowSums(PSD)
PSD <- 0.01*PSD
rowSums(PSD)

## Prep data
## Setup for whole dataset analysis ####
## These are the parameters needed to break the distribution into bins
L <- log(2000)-log(0.03999)         #max(PSD$Size)-min(PSD$Size) but on a logarithmic scale (is this the right base? Not sure)
k <- 1:6                            #number of repetitions
N <- 2^k                            #number of bins per repetition. This CANNOT go above the original # bins
e <- L*(2^-k)                       #bin width
e
N

########### Salat functions: Class : Standard Grid ############

# Use for the standard moment method. 
# 'Data' is a list of vector/grids containing the data aggregated at increasing resolutions;
# 'TheSizes' is a vector containing the resolutions in increasing order.

StdGrid <- function(Data,TheSizes)
{
  SGrid <- list(
    Values = Data,
    Sizes = TheSizes
  )
  class(SGrid) <- append(class(SGrid), c("StdGrid","Moment"))
  return(SGrid)
}

MMoment.StdGrid <- function(theObject,q,dd)
{
  sizes <- theObject$Sizes
  data <- theObject$Values
  Zed <- matrix(rep(0,length(q)*length(sizes)),ncol=length(sizes),byrow=TRUE)
  Tau <- rep(0,length(q))
  Rsq <- rep(0, length(q))
  D <- rep(0,length(q))
  alpha <- rep(0,length(q))
  falpha <- rep(0,length(q))
  for(i in 1:length(q)){
    for(j in 1:length(sizes)){
      Zed[i,j] <- sum(data[[j]][which(data[[j]] != 0)]^q[i])
    }
    fit <- lm(log(Zed[i,]) ~ log(sizes))
    Tau[i] <- fit$coefficients[[2]]
    Rsq[i] <- summary(fit)$r.squared
    D[i] <- Tau[i]/(q[i]-1)
  }
  for(i in 2:(length(q)-1)){
    alpha[i] <- (Tau[i+1]-Tau[i-1])/(q[i+1]-q[i-1])
    falpha[i] <- q[i]*alpha[i]-Tau[i]
  }
  alpha[1] <- (Tau[2]-Tau[1])/(q[2]-q[1])
  alpha[length(q)] <- (Tau[length(q)]-Tau[length(q)-1])/(q[length(q)]-q[length(q)-1])
  falpha[1] <- q[1]*alpha[1]-Tau[1]
  falpha[length(q)] <- q[length(q)]*alpha[length(q)]-Tau[length(q)]
  D[match(1,q)] <- (D[match(1,q)+1]+D[match(1,q)-1])/2 #max(falpha)
  Rsq[match(1,q)] <- NA
  Result <- data.frame(q=q,Dq=D,Rsq=Rsq,alpha=alpha,falpha=falpha)
  return(Result)
}


## create blank dataframe for final result
Result <- data.frame(REP_ID = character(0),
                     q = character(0),
                     Dq = character(0),
                     Rsq = character(0),
                     alpha = character(0),
                     falpha = character(0))

## Actual calculation ####
for (i in 1:nrow(PSD)){
  ## for every rep cut data
  Data.vec <- PSD[i,] # select one rep's worth of data
  
  Data.ls <- list()
  for (j in length(N):1){
    # create an entry in the list with a list of the summed proportion of particle sizes for each size class
    Data.ls[[j]] <-  unname(sapply(split(as.numeric(Data.vec), factor(cut(1:length(Data.vec), N[j]))), sum))
    
  }
  Sizes <- e/(2*max(e))
  
  # specify q values I am interested in 
  q <- seq(-5,5,0.2)
  
  ## for every rep use code from Salat to calculate Dq, alpha and falpha
  Res.temp <- MMoment.StdGrid(StdGrid(Data.ls, Sizes), q, 1) ## Use Salat's code here
  
  Res.temp <- cbind(data.frame(REP_ID = rep(rownames(PSD)[i], length(q))),
                    Res.temp)
  ## bind every rep together into datasheet
  Result <- rbind(Result, Res.temp)
  
}
View(Result)
dim(Result)
# 17085 6

summary(Result)
# REP_ID            q              Dq              Rsq             alpha            falpha        
# 10170X1:   51   Min.   :-5.0   Min.   :0.2999   Min.   :0.5496   Min.   :0.2443   Min.   :-0.35353  
# 10170X4:   51   1st Qu.:-2.6   1st Qu.:0.8872   1st Qu.:0.9782   1st Qu.:0.8646   1st Qu.:-0.09036  
# 10170X5:   51   Median : 0.0   Median :0.9968   Median :0.9930   Median :1.1899   Median : 0.61104  
# 10412X3:   51   Mean   : 0.0   Mean   :1.6363   Mean   :0.9781   Mean   :2.1689   Mean   : 0.40298  
# 10412X4:   51   3rd Qu.: 2.6   3rd Qu.:2.2234   3rd Qu.:0.9985   3rd Qu.:3.0508   3rd Qu.: 0.84967  
# 10498X1:   51   Max.   : 5.0   Max.   :5.6908   Max.   :1.0000   Max.   :6.9279   Max.   : 1.00000  
# (Other):16779                                   NA's   :335                                        
hist(Result$Rsq)

summary(filter(Result, Rsq<0.9))
# REP_ID          q                Dq              Rsq             alpha            falpha        
# 12255X3: 24   Min.   :-5.000   Min.   :0.2999   Min.   :0.5496   Min.   :0.2443   Min.   :-0.33851  
# 22414X4: 24   1st Qu.:-4.000   1st Qu.:2.3705   1st Qu.:0.8508   1st Qu.:3.8755   1st Qu.:-0.23951  
# 23254X2: 24   Median :-2.800   Median :2.9836   Median :0.8615   Median :4.2110   Median :-0.15882  
# 25526X3: 24   Mean   :-2.716   Mean   :2.9778   Mean   :0.8551   Mean   :4.2838   Mean   :-0.13454  
# 31230X3: 24   3rd Qu.:-1.800   3rd Qu.:3.4557   3rd Qu.:0.8673   3rd Qu.:4.8547   3rd Qu.:-0.06104  
# 32721X3: 24   Max.   : 5.000   Max.   :5.6908   Max.   :0.9000   Max.   :6.9279   Max.   : 0.87557  
# (Other):934                                                                                                                                                                               
dim(filter(Result, Rsq<0.9))
# [1] 1078   6
hist(filter(Result, Rsq<0.9)[,2])
hist(filter(Result, Rsq<0.8)[,2])

# many samples at d=1 have a negligible Rsq, or none at all. D1 is not calculated by lm but seems to be
# the average of D0 and D2 which seems weird but gets rid of that issue

Result_filrsq <- filter(Result, Rsq > 0.9 | is.na(Rsq))
summary(Result_filrsq)
REPS_unique <- unique(Result_filrsq$REP_ID)

pdf("Dq against q v2.pdf", width=7, height=147)
par(mfrow=c(84,4), mar = c(4,4,4,2))
for (i in 1:length(REPS_unique)){
  plot(Dq ~ q, filter(Result_filrsq, REP_ID == REPS_unique[i]), ylim = c(0, 6), xlim=c(-5,5), type="b", 
       pch=16, cex=0.5, main = REPS_unique[i])
  abline(h=1, col="gray", lwd=0.5)
}
dev.off()

pdf("falpha vs alpha.pdf", width=7, height=147)
par(mfrow=c(84,4), mar = c(4,4,4,2))
for (i in 1:length(REPS_unique)){
  plot(falpha ~ alpha, filter(Result_filrsq, REP_ID == REPS_unique[i]), ylim = c(-0.4,1), xlim=c(0,7), 
       type="b", pch=16, cex=0.5, main=REPS_unique[i])
}
dev.off()

## cutting axes
pdf("falpha vs alpha v2.pdf", width=7, height=147)
par(mfrow=c(84,4), mar = c(4,4,4,2))
for (i in 1:length(REPS_unique)){
  plot(falpha ~ alpha, filter(Result_filrsq, REP_ID == REPS_unique[i]), ylim = c(0,1), xlim=c(0,4), 
       type="b", pch=16, cex=0.5, main=REPS_unique[i])
}
dev.off()

## correlation network with soil size bins ####
library(qgraph)
library(Hmisc)

summary(Result_filrsq)
Res_fil_wide <- Result_filrsq %>% select(REP_ID, q, Dq) %>%
  filter(q %in% -5:5) %>%
  spread(key=q, value=Dq)

Res_fil_wide_cor_dt <- select(GMEP, REP_ID, X0.0438996:SAND) %>% merge(Res_fil_wide, by="REP_ID")

cor <- rcorr(as.matrix(Res_fil_wide_cor_dt[,2:131]), type="spearman")
n <- 130
0.05/(n*n-n) #2.981515e-06
corsel <- cor$P < 2.981515e-06
cor2 <- matrix(nrow=n,ncol=n)
for (i in 1:n){
  for (j in 1:n){
    ifelse(is.na(corsel[i,j]), cor2[i,j] <- NA, 
           ifelse(corsel[i,j]==F, cor2[i,j] <- NA,cor2[i,j] <- cor$r[i,j]))
  }
}
cor2

colnames(cor2) <- colnames(Res_fil_wide_cor_dt)[2:131]
rownames(cor2) <- colnames(Res_fil_wide_cor_dt)[2:131]
colnames(cor2)
cor2[120:130,c("CLAY","SILT","SAND")] ##find correlations between MF and SSC only
#          CLAY       SILT       SAND
# -5  0.6297889         NA -0.3316516
# -4  0.6269396         NA -0.3281596
# -3  0.6216446         NA -0.3230063
# -2  0.6159378         NA -0.3134073
# -1  0.5984458         NA -0.3144124
# 0  -0.6108842 -0.4029430  0.6226100
# 1  -0.3632962 -0.4033174  0.5422590
# 2          NA -0.3406771  0.4404324
# 3          NA -0.2942029  0.3707512
# 4          NA -0.2587706  0.3210079
# 5          NA         NA  0.2850174
Size <- as.numeric(matrix(unlist(strsplit(colnames(cor2)[1:116],"X")), ncol=2, byrow=T)[,2])
Sizelab <- as.character(signif(Size,2))
qlab <- paste("q",colnames(cor2)[120:130], sep="=")
labs <- c(Sizelab, colnames(cor2)[117:119], qlab)

palette <- colorRampPalette(colors = c("firebrick1","lightgoldenrod1","dodgerblue1"))
col <- c(palette(116),"firebrick1", "lightgoldenrod1","dodgerblue1",
         rep("hotpink",11))
qgraph(cor2, layout="spring", labels=labs, color=col)
qgraph(cor2, layout="spring", labels=labs, color=col, filetype="png", 
       filename=paste("PSD network Bonferroni multifractals colour gradient",Sys.Date(),sep=" "), 
       width=20, height=16)

corsel <- abs(cor$r) > 0.5
cor2 <- matrix(nrow=n,ncol=n)
for (i in 1:n){
  for (j in 1:n){
    ifelse(is.na(corsel[i,j]), cor2[i,j] <- NA, 
           ifelse(corsel[i,j]==F, cor2[i,j] <- NA,cor2[i,j] <- cor$r[i,j]))
  }
}
cor2
qgraph(cor2, layout="spring", labels=labs, color=col)
qgraph(cor2, layout="spring", labels=labs, color=col, filetype="png", 
       filename=paste("PSD network r 0.5 multifractals colour gradient",Sys.Date(),sep=" "), 
       width=20, height=16, vsize = 2)


cor_pos <- ifelse(cor2 > 0, cor$r, 0)
cor_pos <- ifelse(cor$P < 2.981515e-06, cor_pos, 0)
cor_pos <- cor_pos[1:116,1:116] # remove multifractals and summary values
qgraph(cor_pos, layout = "spring")
graph <- graph_from_adjacency_matrix(cor_pos, mode="undirected", weighted = TRUE)
plot(graph)

clusters.eb <- cluster_edge_betweenness(graph) #warning about modularity treating links as similarities but eb treating as distances
membership(clusters.eb)

cluster.wt <- cluster_walktrap(graph)
membership(cluster.wt)

cluster.sg <- cluster_spinglass(graph)
membership(cluster.sg)

membership(cluster.wt)
# X0.0438996 X0.0481915 X0.0529029 X0.0580749 X0.0637526 X0.0699854 X0.0768275 X0.0843385 X0.0925838  X0.101635  X0.111572  X0.122479 
# 3          3          3          3          3          3          3          3          3          3          3          3 
# X0.134454  X0.147598  X0.162028  X0.177869  X0.195258  X0.214348  X0.235303  X0.258308  X0.283561  X0.311283  X0.341716  X0.375124 
# 2          2          2          2          2          2          2          2          2          2          2          2 
# X0.411798  X0.452057  X0.496252  X0.544768  X0.598027  X0.656493  X0.720675  X0.791132  X0.868477  X0.953383   X1.04659   X1.14891 
# 2          2          2          2          2          2          2          2          2          2          2          2 
# X1.26123   X1.38454    X1.5199   X1.66849   X1.83161   X2.01068   X2.20725   X2.42304   X2.65993   X2.91998   X3.20545   X3.51883 
# 2          2          2          2          2          2          2          2          2          2          2          2 
# X3.86284   X4.24049   X4.65506   X5.11017   X5.60976    X6.1582   X6.76025   X7.42117   X8.14669   X8.94315   X9.81748   X10.7773 
# 2          2          2          2          2          2          2          2          2          2          2          2 
# X11.8309   X12.9876   X14.2573   X15.6512   X17.1813    X18.861    X20.705   X22.7292   X24.9513   X27.3906   X30.0685   X33.0081 
# 2          2          2          2          2          2          2          2          3          3          3          3 
# X36.2352   X39.7777   X43.6665   X47.9356    X52.622   X57.7666   X63.4141   X69.6138   X76.4196   X83.8907   X92.0923   X101.096 
# 3          3          3          3          3          3          1          1          1          1          1          1 
# X110.979   X121.829    X133.74   X146.815   X161.168   X176.925   X194.222    X213.21   X234.054   X256.936   X282.056   X309.631 
# 1          1          1          1          1          1          1          1          1          1          1          1 
# X339.902   X373.132   X409.611   X449.657   X493.617   X541.876   X594.852   X653.008   X716.849   X786.932   X863.866   X948.322 
# 1          1          1          1          1          1          1          1          1          1          1          1 
# X1041.03   X1142.81   X1254.54   X1377.19   X1511.83   X1659.63   X1821.88      X2000 
# 1          1          1          1          1          1          1          1
modularity(cluster.wt)
# [1] 0.4138004
sizes(cluster.wt)
# Community sizes
# 1  2  3 
# 38 56 22 
## spinglass can deal with negative weights
graph.neg <- graph_from_adjacency_matrix(cor$r[1:116,1:116], mode = "undirected", weighted = TRUE)

clus2.sg <- cluster_spinglass(graph.neg, implementation = "neg", gamma.minus = 0.1)
membership(clus2.sg)

## I think it is treating negative weights as a link not as a repellent so I'm going to ignore this

## Soil texture ternary diagram ####
library(soiltexture)
TT.plot(tri.data = Res_fil_wide_cor_dt, main = "Texture classes",
        text.sum = 100, base.css.ps.lim = c(0,2.2,63,2000),
        class.sys = "UK.SSEW.TT", pch=16, frame.bg.col = "white", grid.col = "gray90",
        class.lab.col = "gray60", class.line.col	= "gray50", cex=0.7)

### Comparison to other data ####

## filtering results by habitat ####
BH <- select(Ancil, REP_ID, BH = BROAD_HABITAT_SURVEYOR)
summary(BH$BH)

BH_PSD <- BH[BH$REP_ID %in% Result$REP_ID,]
BH_PSD <- droplevels(BH_PSD)
summary(BH_PSD$BH)


## arable ####
Arable_reps <- as.character(filter(BH_PSD, BH == "Arable and horticultural")$REP_ID)


theme_set(theme_bw())

Arable_PSD <- Result_filrsq[Result_filrsq$REP_ID %in% Arable_reps,]
Arable_PSD <- select(Soil, REP_ID, C_FE_CTOTAL, C_B_PH_CACL2, TOTAL_MITES, TOTAL_INVERT) %>% 
  merge(select(MIC_RICH, REP_ID, FUNG_R_RICH_BLAST, BACT_R_RICH), by="REP_ID", all=T) %>%
  merge(Arable_PSD, by="REP_ID", all.y=T)
summary(Arable_PSD)

d1 <- ggplot(Arable_PSD, aes(x=q, y=Dq, group=REP_ID, col=C_FE_CTOTAL)) + geom_line(lwd=1.5) +
  scale_color_gradient(low="#a1dab4", high="#253494",name="Carbon") + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black"))


d2 <- ggplot(Arable_PSD, aes(x=q, y=Dq, group=REP_ID, col=C_B_PH_CACL2)) + geom_line(lwd=1.5) +
  scale_color_gradient(low="red", high="green",name="pH") + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black"))

d3 <- ggplot(Arable_PSD, aes(x=q, y=Dq, group=REP_ID, col=BACT_R_RICH)) + geom_line(lwd=1.5) +
  scale_color_gradient(low="#a1dab4", high="#253494",name="Bacterial\nRichness") + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black"))


d4 <- ggplot(Arable_PSD, aes(x=q, y=Dq, group=REP_ID, col=FUNG_R_RICH_BLAST)) + geom_line(lwd=1.5) +
  scale_color_gradient(low="#a1dab4", high="#253494",name="Fungal\nRichness") + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black"))

grid.arrange(d1,d2,d3,d4)

## falpha vs alpha
f1 <- ggplot(Arable_PSD, aes(x=alpha, y=falpha, group=REP_ID, col=C_FE_CTOTAL)) + geom_line(lwd=1.5) +
  scale_color_gradient(low="#a1dab4", high="#253494",name="Carbon") + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + scale_y_continuous(limits = c(0,1))


f2 <- ggplot(Arable_PSD, aes(x=alpha, y=falpha, group=REP_ID, col=C_B_PH_CACL2)) + geom_line(lwd=1.5) +
  scale_color_gradient(low="#b30000", high="darkolivegreen2",name="pH") + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + scale_y_continuous(limits = c(0,1))

f3 <- ggplot(Arable_PSD, aes(x=alpha, y=falpha, group=REP_ID, col=BACT_R_RICH)) + geom_line(lwd=1.5) +
  scale_color_gradient(low="#a1dab4", high="#253494",name="Bacteria\nRichness") + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + scale_y_continuous(limits = c(0,1))


f4 <- ggplot(Arable_PSD, aes(x=alpha, y=falpha, group=REP_ID, col=FUNG_R_RICH_BLAST)) + geom_line(lwd=1.5) +
  scale_color_gradient(low="#a1dab4", high="#253494",name="Fungi\nRichness") + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + scale_y_continuous(limits = c(0,1))

grid.arrange(f1,f2,f3,f4)

## plot D0, D1, D2 against carbon, pH and bacterial richness
summary(BH_PSD$BH)
BH_PSD$Habitat <- recode_factor(BH_PSD$BH,
                                "Arable and horticultural" = "Arable",
                                "Improved Grassland" = "Impr Grass",
                                "Neutral Grassland" = "Neutr Grass",
                                "Acid Grassland" = "Acid Grass",
                                "Broadleaved, mixed and yew woodland" = "Broadleaved",
                                "Coniferous Woodland" = "Conifer",
                                "Fen, Marsh and Swamp" = "Heath",
                                "Bracken" = "Heath",
                                "Dwarf Shrub Heath" = "Heath",
                                "Bog" = "Heath",
                                "Calcareous Grassland" = "Impr Grass",
                                "Supra-littoral rock" = "Other",
                                "Supra-littoral sediment" = "Other")
summary(BH_PSD$Habitat)

Res_fil_mer <- select(Soil, REP_ID, C_FE_CTOTAL, C_B_PH_CACL2, TOTAL_MITES, TOTAL_INVERT) %>% 
  merge(select(MIC_RICH, REP_ID, FUNG_R_RICH_BLAST, BACT_R_RICH), by="REP_ID", all=T) %>%
  merge(Result_filrsq, by="REP_ID", all.y=T) %>% merge(BH_PSD, by="REP_ID", all.x=T)
summary(Res_fil_mer)

# remove littoral and NA habitats
Res_fil_mer <- filter(Res_fil_mer, Habitat != "Other")
Res_fil_mer <- Res_fil_mer[!is.na(Res_fil_mer$Habitat),]

## simple plots by habitat
## Dq vs q
ggplot(Res_fil_mer, aes(x=q, y=Dq, group=REP_ID)) + geom_line(alpha=0.2) + 
  facet_wrap(~Habitat)
ggsave("Dq vs q by habitat.png", device = "png")

ggplot(Res_fil_mer, aes(x=q, y=Dq, group=REP_ID)) + geom_line(alpha=0.2) + 
  facet_wrap(~Habitat) + scale_x_continuous(limits = c(-0.5,2.5)) + scale_y_continuous(limits = c(0,2))
ggsave("Dq vs q by habitat zoomed.png", device = "png")

## falpha vs alpha
ggplot(Res_fil_mer, aes(x=alpha, y=falpha, group=REP_ID)) + geom_line(alpha=0.2) + 
  facet_wrap(~Habitat)
ggsave("falpha by alpha by habitat.png", device = "png")

ggplot(Res_fil_mer, aes(x=alpha, y=falpha, group=REP_ID)) + geom_line(alpha=0.2) + 
  facet_wrap(~Habitat) + scale_y_continuous(limits = c(0,1))
ggsave("falpha by alpha by habitat zoomed.png", device = "png")

# average by habitat
res_fil_mer_hab <- Res_fil_mer %>%
  group_by(Habitat, q) %>%
  summarise(mean_Dq = mean(Dq),
            upper_se = mean(Dq) + sd(Dq)/sqrt(length(Dq)),
            lower_se = mean(Dq) - sd(Dq)/sqrt(length(Dq)),
            upper_sd = mean(Dq) + sd(Dq),
            lower_sd = mean(Dq) - sd(Dq))
summary(res_fil_mer_hab)

ggplot(res_fil_mer_hab, aes(x = q + 0.005*as.numeric(Habitat), 
                            y = mean_Dq, colour = Habitat)) +
  geom_line(lwd = 2) +
  # geom_errorbar(aes(ymax = upper_se, ymin = lower_se)) +
  geom_linerange(aes(ymax = upper_sd, ymin = lower_sd)) +
  scale_color_brewer(palette = "Set2")

## simple plots by carbon
## Dq vs q
ggplot(Res_fil_mer, aes(x=q, y=Dq, group=REP_ID)) + geom_line(alpha=0.2) + 
  facet_wrap(~cut_number(C_FE_CTOTAL,6)) + ggtitle("Dq vs q faceted by carbon")
ggsave("Dq vs q by carbon.png", device = "png")

ggplot(Res_fil_mer, aes(x=q, y=Dq, group=REP_ID)) + geom_line(alpha=0.2) + 
  facet_wrap(~cut_number(C_B_PH_CACL2,6)) + ggtitle("Dq vs q faceted by pH")
ggsave("Dq vs q by pH.png", device = "png")

ggplot(Res_fil_mer, aes(x=q, y=Dq, group=REP_ID)) + geom_line(alpha=0.2) + 
  facet_wrap(~cut_number(BACT_R_RICH,8)) + ggtitle("Dq vs q faceted by bacteria richness")
ggsave("Dq vs q by bact.png", device = "png")

ggplot(Res_fil_mer, aes(x=q, y=Dq, group=REP_ID)) + geom_line(alpha=0.2) + 
  facet_wrap(~cut_number(FUNG_R_RICH_BLAST,8)) + ggtitle("Dq vs q faceted by fungal richness")
ggsave("Dq vs q by fungi.png", device = "png")

ggplot(Res_fil_mer, aes(x=q, y=Dq, group=REP_ID)) + geom_line(alpha=0.2) + 
  facet_wrap(~cut_number(TOTAL_INVERT,8)) + ggtitle("Dq vs q faceted by total invert catch")
ggsave("Dq vs q by invert.png", device = "png")


### Violin plots of D0, D1 and D2 by habitat
v1 <- ggplot(filter(Res_fil_mer, q==0), aes(x=Habitat, y=Dq)) + geom_violin(fill="#6CA6CD") +
  labs(y="D0", x="")
summary(lm(Dq ~ Habitat, filter(Res_fil_mer, q==0))) #R2 = 0.03981, p = 0.03994

v2 <- ggplot(filter(Res_fil_mer, q==1), aes(x=Habitat, y=Dq)) + geom_violin(fill="#FF6A6A") +
  labs(y="D1", x="")
summary(lm(Dq ~ Habitat, filter(Res_fil_mer, q==1))) #R2 = 0.02495, p = 0.2229

v3 <- ggplot(filter(Res_fil_mer, q==2), aes(x=Habitat, y=Dq)) + geom_violin(fill="#66CDAA") +
  labs(y="D2", x="")
summary(lm(Dq ~ Habitat, filter(Res_fil_mer, q==2))) #R2 = 0.003525, p = 0.9794

grid.arrange(v1,v2,v3)
ggsave("Violin plots D by habitat.png", plot = arrangeGrob(v1,v2,v3), device="png", width = 6, height=6, units="in")

## wide format data for ratio calcs ####
## limit to q = 0, 1 or 2
Res_fil_mer_wideD <- filter(Res_fil_mer, q == 0 |q == 1|q == 2) %>% select(-Rsq, -alpha, -falpha) %>% 
  spread(key=q, value=Dq, sep="D")
summary(Res_fil_mer_wideD)

Res_fil_mer_wideD$D1_D0 <- Res_fil_mer_wideD$qD1/Res_fil_mer_wideD$qD0
Res_fil_mer_wideD$D2_D1 <- Res_fil_mer_wideD$qD2/Res_fil_mer_wideD$qD1

psych::multi.hist(Res_fil_mer_wideD[,10:14])

(v4 <- ggplot(Res_fil_mer_wideD, aes(x=Habitat, y=D1_D0)) + geom_violin(fill="#FF7F24") +
    labs(y="D1/D0", x=""))
summary(lm(D1_D0 ~ Habitat, Res_fil_mer_wideD)) #R2 = 0.008456, p = 0.8383

(v5 <- ggplot(Res_fil_mer_wideD, aes(x=Habitat, y=D2_D1)) + geom_violin(fill="#7A67EE") +
    labs(y="D2/D1", x=""))
summary(lm(D2_D1 ~ Habitat, Res_fil_mer_wideD)) #R2 = 0.004011, p = 0.9713


grid.arrange(v1,v2,v3, v4, v5, layout_matrix = rbind(c(1, 2),
                                                     c(3, NA),
                                                     c(4, 5)))
ggsave("Violin plots D by habitat 2.png", plot = arrangeGrob(v1,v2,v3,v4,v5, layout_matrix = rbind(c(1, 2),
                                                                                                   c(3, NA),
                                                                                                   c(4, 5))), 
       device="png", width = 12, height=8, units="in")

## Mixed effect modelling by habitat ####
Res_fil_mer_wideD <- separate(Res_fil_mer_wideD, REP_ID, "SQNUM", sep="X", remove = FALSE, convert = TRUE)
summary(Res_fil_mer_wideD)
Res_fil_mer_wideD$SQNUM <- as.factor(Res_fil_mer_wideD$SQNUM)

summary(lm(Dq ~ Habitat, filter(Res_fil_mer, q==0))) #R2 = 0.03981, p = 0.03994


# Carbon vs D0
ggplot(filter(Res_fil_mer, q==0), aes(x=Dq, y=C_FE_CTOTAL, col=Habitat)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Carbon (%)", x="D0") +
  scale_color_brewer(palette = "Dark2")

ggplot(filter(Res_fil_mer, q==0), aes(x=Dq, y=C_FE_CTOTAL)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Carbon (%)", x="D0") +
  facet_wrap(~Habitat)

# pH vs D0
ggplot(filter(Res_fil_mer, q==0), aes(x=Dq, y=C_B_PH_CACL2, col=Habitat)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="pH", x="D0") +
  scale_color_brewer(palette = "Dark2")

ggplot(filter(Res_fil_mer, q==0), aes(x=Dq, y=C_B_PH_CACL2)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="pH", x="D0") +
  facet_wrap(~Habitat)




# Bacteria vs D0
ggplot(filter(Res_fil_mer, q==0), aes(x=Dq, y=BACT_R_RICH, col=Habitat)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D0") +
  scale_color_brewer(palette = "Dark2")

ggplot(filter(Res_fil_mer, q==0), aes(x=Dq, y=BACT_R_RICH)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D0") +
  facet_wrap(~Habitat)

(b0 <- ggplot(filter(Res_fil_mer, q==0), aes(x=Dq, y=BACT_R_RICH)) + geom_point(cex = 2) + 
  geom_smooth(method="lm") +
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D0"))


# Bacteria vs D1
ggplot(filter(Res_fil_mer, q==1), aes(x=Dq, y=BACT_R_RICH, col=Habitat)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D1") +
  scale_color_brewer(palette = "Dark2")

ggplot(filter(Res_fil_mer, q==1), aes(x=Dq, y=BACT_R_RICH)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D1") +
  facet_wrap(~Habitat)

b1 <- ggplot(filter(Res_fil_mer, q==1), aes(x=Dq, y=BACT_R_RICH)) + geom_point(cex = 2) + 
  geom_smooth(method="lm") +
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D1")


summary(lm(BACT_R_RICH ~ Dq*Habitat, filter(Res_fil_mer, q==1)))

filter(Res_fil_mer, q==1) %>% summarise(r = cor(BACT_R_RICH, Dq, use="complete.obs")) # 0.02938608
filter(Res_fil_mer, q==2) %>% summarise(r = rcorr(BACT_R_RICH, Dq, use="complete.obs")$r) #0.06332883

filter(Res_fil_mer, q==2) %>% summarise(r = rcorr(BACT_R_RICH, Dq, type="spearman")$r[2,1],
                                        P = rcorr(BACT_R_RICH, Dq, type="spearman")$P[2,1])

# Bacteria vs D2
ggplot(filter(Res_fil_mer, q==2), aes(x=Dq, y=BACT_R_RICH, col=Habitat)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D2") +
  scale_color_brewer(palette = "Dark2")

ggplot(filter(Res_fil_mer, q==2), aes(x=Dq, y=BACT_R_RICH)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D2") +
  facet_wrap(~Habitat) + geom_smooth(method="lm")

b2 <- ggplot(filter(Res_fil_mer, q==2), aes(x=Dq, y=BACT_R_RICH)) + geom_point(cex = 2) + 
  geom_smooth(method="lm") +
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D2")

summary(lm(BACT_R_RICH ~ Dq*Habitat, filter(Res_fil_mer, q==2)))

## Plots for paper ####
f0 <- ggplot(filter(Res_fil_mer, q==0), aes(x=Dq, y=FUNG_R_RICH_BLAST)) + geom_point(cex = 2) + 
  geom_smooth(method="lm") +
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Fungal Richness", x="D0")

f1 <- ggplot(filter(Res_fil_mer, q==1), aes(x=Dq, y=FUNG_R_RICH_BLAST)) + geom_point(cex = 2) + 
  geom_smooth(method="lm") +
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Fungal Richness", x="D1")

f2 <- ggplot(filter(Res_fil_mer, q==2), aes(x=Dq, y=FUNG_R_RICH_BLAST)) + geom_point(cex = 2) + 
  geom_smooth(method="lm") +
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Fungal Richness", x="D2")

library(gridExtra)
grid.arrange(b0,b1,b2,f0,f1,f2, ncol = 3)
ggsave("Bacteria Fungi D parameters lm.png", plot=arrangeGrob(b0,b1,b2,f0,f1,f2, ncol = 3),
       device = "png", path = "./Graphs/", width=30, height = 20, units="cm")
ggsave("Bacteria Fungi D parameters.pdf", plot=arrangeGrob(b0,b1,b2,f0,f1,f2, ncol = 3),
       device = "pdf", path = "./Graphs/", width=30, height = 20, units="cm")

b01 <- ggplot(Res_fil_mer_wideD, aes(x=D1_D0, y=BACT_R_RICH)) + geom_point(cex = 2) + 
  geom_smooth(method="lm") +
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D1/D0")
b12 <- ggplot(Res_fil_mer_wideD, aes(x=D2_D1, y=BACT_R_RICH)) + geom_point(cex = 2) + 
  geom_smooth(method="lm") +
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D2/D1")

f01 <- ggplot(Res_fil_mer_wideD, aes(x=D1_D0, y=FUNG_R_RICH_BLAST)) + geom_point(cex = 2) + 
  geom_smooth(method="lm") +
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Fungal Richness", x="D1/D0")

f12 <- ggplot(Res_fil_mer_wideD, aes(x=D2_D1, y=FUNG_R_RICH_BLAST)) + geom_point(cex = 2) + 
  geom_smooth(method="lm") +
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Fungal Richness", x="D2/D1")
grid.arrange(b01,b12,f01,f12, ncol = 2)
ggsave("Bacteria Fungi D ratios lm.png", plot=arrangeGrob(b01,b12,f01,f12, ncol = 2),
       device = "png", path = "./Graphs/", width=20, height = 20, units="cm")


b01 <- ggplot(Res_fil_mer_wideD, aes(x=D1_D0, y=BACT_R_RICH)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D1/D0")
b12 <- ggplot(Res_fil_mer_wideD, aes(x=D2_D1, y=BACT_R_RICH)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Bacterial Richness", x="D2/D1")

f01 <- ggplot(Res_fil_mer_wideD, aes(x=D1_D0, y=FUNG_R_RICH_BLAST)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Fungal Richness", x="D1/D0")

f12 <- ggplot(Res_fil_mer_wideD, aes(x=D2_D1, y=FUNG_R_RICH_BLAST)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Fungal Richness", x="D2/D1")
grid.arrange(b01,b12,f01,f12, ncol = 2)
ggsave("Bacteria Fungi D ratios.png", plot=arrangeGrob(b01,b12,f01,f12, ncol = 2),
       device = "png", path = "./Graphs/", width=20, height = 20, units="cm")

## Models ####
BACT_res <- resid(lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL, filter(Res_fil_mer, q==1)))

summary(lm(BACT_res ~ Dq, filter(Res_fil_mer, q==1 & BACT_R_RICH > 0)))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  -2654.1      900.9  -2.946  0.00345 **
#   Dq            2902.8      984.7   2.948  0.00344 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 557.5 on 319 degrees of freedom
# Multiple R-squared:  0.02652,	Adjusted R-squared:  0.02347 
# F-statistic:  8.69 on 1 and 319 DF,  p-value: 0.003436

ggplot(filter(Res_fil_mer, q==1 & BACT_R_RICH > 0), aes(x=Dq, y=BACT_res)) + geom_point() +
  geom_smooth(method="lm") + facet_wrap(~Habitat)

mod.bact <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL + Dq, filter(Res_fil_mer, q==1))
par(mfrow=c(2,2));plot(mod.bact);par(mfrow=c(1,1))

mod.bact <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL + qD1, Res_fil_mer_wideD)
par(mfrow=c(2,2));plot(mod.bact);par(mfrow=c(1,1))

Res_fil_mer_wideD <- na.omit(Res_fil_mer_wideD)

mod.bact <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL + qD1, Res_fil_mer_wideD)
par(mfrow=c(2,2));plot(mod.bact);par(mfrow=c(1,1))


## identify outlier
Res_fil_mer_wideD[290,]
psych::multi.hist(Res_fil_mer_wideD[,2:7])
psych::multi.hist(Res_fil_mer_wideD[,10:14])

mod.bact <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL + qD1, Res_fil_mer_wideD)
par(mfrow=c(2,2));plot(mod.bact, labels.id = Res_fil_mer_wideD$REP_ID);par(mfrow=c(1,1))
Res_fil_mer_wideD[Res_fil_mer_wideD$REP_ID=="46364X5",]
# REP_ID C_FE_CTOTAL C_B_PH_CACL2 TOTAL_MITES TOTAL_INVERT FUNG_R_RICH_BLAST BACT_R_RICH                   BH Habitat
# 290 46364X5        5.57         7.55           0            0          92.32172    2883.288 Fen, Marsh and Swamp   Heath
# qD0       qD1       qD2     D1_D0     D2_D1
# 290 0.9967543 0.9339973 0.9156323 0.9370387 0.9803372

Res_fil_mer_wideD <- filter(Res_fil_mer_wideD, REP_ID != "46364X5")


mod.bact <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL, Res_fil_mer_wideD[])
summary(mod.bact)
par(mfrow=c(2,2));plot(mod.bact);par(mfrow=c(1,1))

mod.bact1 <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL + qD1, Res_fil_mer_wideD[])
summary(mod.bact1)
par(mfrow=c(2,2));plot(mod.bact1);par(mfrow=c(1,1))

mod.bact2 <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL + qD0 + qD1 + qD2, Res_fil_mer_wideD[])
summary(mod.bact2)
par(mfrow=c(2,2));plot(mod.bact2);par(mfrow=c(1,1))

mod.bact3 <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL + D1_D0, Res_fil_mer_wideD)
summary(mod.bact3)
par(mfrow=c(2,2));plot(mod.bact3);par(mfrow=c(1,1))

library(AICcmodavg)

aictab(list(mod.bact, mod.bact1, mod.bact2, mod.bact3), modnames = c("No D", "D1", "D012","D1/D0"))
# K    AICc Delta_AICc AICcWt Cum.Wt       LL
# D012  7 4795.02       0.00   0.86   0.86 -2390.33
# D1    5 4798.92       3.89   0.12   0.98 -2394.36
# D1/D0 5 4802.42       7.40   0.02   1.00 -2396.11
# No D  4 4807.62      12.59   0.00   1.00 -2399.74

psych::pairs.panels(Res_fil_mer_wideD[,c(2:7,10:14)])

car::vif(mod.bact2)
# C_B_PH_CACL2  C_FE_CTOTAL          qD0          qD1          qD2 
# 1.171194     1.209523     2.136269    16.806995    14.732028 
ar::vif(mod.bact1)
# C_B_PH_CACL2  C_FE_CTOTAL          qD1 
# 1.169930     1.203374     1.043497 
car::vif(mod.bact3)
# C_B_PH_CACL2  C_FE_CTOTAL        D1_D0 
# 1.168666     1.194361     1.030115 

mod.bact4 <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL + qD1 + qD2, Res_fil_mer_wideD)
summary(mod.bact4)
par(mfrow=c(2,2));plot(mod.bact4);par(mfrow=c(1,1))

car::vif(mod.bact4)
# C_B_PH_CACL2  C_FE_CTOTAL          qD1          qD2 
# 1.170009     1.207907     8.018582     7.903631 

mod.bact5 <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL + qD2, Res_fil_mer_wideD)
summary(mod.bact5)
par(mfrow=c(2,2));plot(mod.bact5);par(mfrow=c(1,1))

car::vif(mod.bact5)

aictab(list(mod.bact, mod.bact1, mod.bact2, mod.bact3, mod.bact4, mod.bact5), 
       modnames = c("No D", "D1", "D012","D1/D0","D12","D2"))


## run with D1
mod.bact.h <- lm(BACT_R_RICH ~ C_B_PH_CACL2 + C_FE_CTOTAL + qD1 + Habitat, Res_fil_mer_wideD[])
summary(mod.bact.h)
par(mfrow=c(2,2));plot(mod.bact.h);par(mfrow=c(1,1))

ggplot(Res_fil_mer_wideD, aes(x=qD1, y=resid(mod.bact))) + geom_point() + 
  labs(x="D1", y="Residual bacteria richness") + 
  theme(text = element_text(size=15), axis.text = element_text(colour="black"))
ggsave("Resid bact by D1.png", device="png", width=12, height=10, units="cm")

par(mfrow=c(2,2));plot(lm(resid(mod.bact) ~ qD1, Res_fil_mer_wideD));par(mfrow=c(1,1))



# Fungi against D1
ggplot(filter(Res_fil_mer, q==1), aes(x=Dq, y=FUNG_R_RICH_BLAST, col=Habitat)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Funga; Richness", x="D1") +
  scale_color_brewer(palette = "Dark2")

ggplot(filter(Res_fil_mer, q==1), aes(x=Dq, y=FUNG_R_RICH_BLAST)) + geom_point(cex = 2) + 
  theme(text = element_text(size=15, colour = "black"),
        axis.text = element_text(colour = "black")) + labs(y="Fungal Richness", x="D1") +
  facet_wrap(~Habitat)

fung.mod <- lm(FUNG_R_RICH_BLAST ~ C_FE_CTOTAL + C_B_PH_CACL2, Res_fil_mer_wideD)
summary(fung.mod)

fung.mod1 <- lm(FUNG_R_RICH_BLAST ~ C_FE_CTOTAL + C_B_PH_CACL2 + qD1, Res_fil_mer_wideD)
summary(fung.mod1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  124.8163    91.3675   1.366   0.1729    
# C_FE_CTOTAL   -1.7145     0.5623  -3.049   0.0025 ** 
#   C_B_PH_CACL2  17.8402     3.8048   4.689 4.14e-06 ***
#   qD1           -1.0477    97.5934  -0.011   0.9914    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 53.7 on 307 degrees of freedom
# Multiple R-squared:  0.1396,	Adjusted R-squared:  0.1312 
# F-statistic:  16.6 on 3 and 307 DF,  p-value: 5.073e-10


ggplot(Res_fil_mer_wideD, aes(x=qD1, y=resid(fung.mod))) + geom_point() + 
  labs(x="D1", y="Residual fungal richness") + 
  theme(text = element_text(size=15), axis.text = element_text(colour="black"))
ggsave("Resid fung by D1.png", device="png", width=12, height=10, units="cm")

par(mfrow=c(2,2));plot(lm(resid(fung.mod) ~ qD1, Res_fil_mer_wideD));par(mfrow=c(1,1))


## Spearman rank correlations for paper ####
cor(Res_fil_mer_wideD[,c(7:8,11:15)], method="spearman")
#                   FUNG_R_RICH_BLAST BACT_R_RICH          qD0         qD1         qD2        D1_D0       D2_D1
# FUNG_R_RICH_BLAST       1.000000000 0.469187285 -0.068391654 -0.04602738 -0.02373118 -0.008993242  0.05209921
# BACT_R_RICH             0.469187285 1.000000000  0.007333572  0.06392251  0.08389211  0.050789497  0.13842444
# qD0                    -0.068391654 0.007333572  1.000000000  0.44259456  0.26641465  0.059072096 -0.09052815
# qD1                    -0.046027383 0.063922510  0.442594556  1.00000000  0.94564959  0.849948537  0.65775135
# qD2                    -0.023731180 0.083892111  0.266414647  0.94564959  1.00000000  0.886304086  0.84016979
# D1_D0                  -0.008993242 0.050789497  0.059072096  0.84994854  0.88630409  1.000000000  0.80350586
# D2_D1                   0.052099208 0.138424437 -0.090528149  0.65775135  0.84016979  0.803505860  1.00000000


# network ####
library(SpiecEasi)
library(igraph)
library(Matrix)

BACT_NET <- BACT_ORD
dim(BACT_NET) #[1]   436 58164

# only samples with psd
BACT_NET <- BACT_NET[substring(rownames(BACT_NET),2) %in% GMEP$REP_ID,]
# only OTUs that appear in 50% of sites
BACT_NET <- BACT_NET[,colSums(BACT_NET>0)>0.8*nrow(BACT_NET)]
dim(BACT_NET) #[1]  359 1265

FUNG_NET <- FUNGI_ORD
dim(FUNG_NET) #[1]  437 8407

FUNG_NET <- FUNG_NET[substring(rownames(FUNG_NET),2) %in% GMEP$REP_ID,]
FUNG_NET <- FUNG_NET[,colSums(FUNG_NET>0)>0.25*nrow(FUNG_NET)]
dim(FUNG_NET) #[1] 361  223
FUNG_NET <- FUNG_NET[substring(rownames(FUNG_NET),2) %in% substring(rownames(BACT_NET),2),]
dim(FUNG_NET) #[1] 359  223

PSD_NET <- GMEP[GMEP$REPEATED < 2,c(1,8:123)]
rownames(PSD_NET) <- paste0("X",PSD_NET$REP_ID)
PSD_NET <- PSD_NET[,-1]
PSD_NET <- PSD_NET[rownames(PSD_NET) %in% rownames(FUNG_NET),]

FUNG_NET <- FUNG_NET[rownames(FUNG_NET) %in% rownames(PSD_NET),]
FUNG_NET <- FUNG_NET[rownames(PSD_NET),]

BACT_NET <- BACT_NET[rownames(BACT_NET) %in% rownames(PSD_NET),]
BACT_NET <- BACT_NET[rownames(PSD_NET),]


all.equal(rownames(PSD_NET),rownames(FUNG_NET))
PSD_NET <- as.matrix(PSD_NET)
PSD_NET <- 1000*round(PSD_NET,3)
PSD_NET <- PSD_NET[,colSums(PSD_NET>0)>0.5*nrow(PSD_NET)]

# network construction
se1 <- spiec.easi(list(PSD_NET, BACT_NET, FUNG_NET), 
                  method = "mb", sel.criterion = "bstars")

refit1 <- getRefit(se1)
elist <- summary(symBeta(getOptBeta(se1), mode = "maxabs"))
diag(elist) <- 0
weights <- elist[,3]

se1.ig <- adj2igraph(refit1,
                     edge.attr = list(weight = weights))

co <- layout_with_fr(se1.ig)
col <- c(rep("dodgerblue",114),rep("gray",1265),rep("firebrick3",223))
plot(se1.ig, layout = co, vertex.label = NA, vertex.color = col,
     vertex.size = 2, alpha = 0.1)

test <- as_adjacency_matrix(se1.ig, attr = "weight", sparse = FALSE)

ncol(PSD_NET) #114
ncol(BACT_NET) #1265
ncol(FUNG_NET) #223
t2 <- test[1:114,115:ncol(test)]
sum(t2 !=0) #307

summary(c(t2))
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.1316688  0.0000000  0.0000000  0.0001723  0.0000000  0.8600919 
hist(rowSums(t2))
summary(rowSums(t2))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.14649  0.00000  0.05153  0.25642  0.33352  1.63847 
sort(rowSums(t2))
# only cols 
t2_no0 <- c(t2)[c(t2) != 0]
summary(t2_no0)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.131669  0.004112  0.024892  0.095219  0.131068  0.860092

test2 <- test
test2[1:114,1:114] <- 0
test2[115:1380,115:1380] <- 0
test2[1381:1602,1381:1602] <- 0

test2.ig <- adj2igraph(test2)
co <- layout_with_fr(test2.ig)
plot(test2.ig, layout = co, vertex.color = col, vertex.label = NA,
     vertex.size = 5)

test3 <- test2
test3[115:1101,1102:1273] <- 0
test3[1102:1273,117:1101] <- 0

rownames(test3) <- c(paste0("P",substring(colnames(PSD_NET),2,5)),
                     paste0("B",substring(colnames(BACT_NET),4)),
                     paste0("F",substring(colnames(FUNG_NET), 4)))
colnames(test3) <- c(paste0("P",substring(colnames(PSD_NET),2,5)),
                     paste0("B",substring(colnames(BACT_NET),4)),
                     paste0("F",substring(colnames(FUNG_NET), 4)))
test4 <- test3[rowSums(test3) != 0, colSums(test3) != 0]

plot(adj2igraph(test4))

library(qgraph)
qgraph(test4)

col_test <- ifelse(grepl("P",rownames(test4)), "red",
                   ifelse(grepl("B", rownames(test4)), "blue", "yellow"))

qgraph(test4, color = col_test)



# spearman rank correlationscolnames(PSD_NET)
testcor <- cor(cbind(PSD_NET, FUNG_NET), method = "spearman")
dim(testcor)
hist(testcor)
length(c(testcor)[abs(c(testcor))>0.5])
length(c(testcor)[abs(c(testcor))>0.8])

# set all texture to texture and fungi to fungi correlations to 0
testcor[1:114,1:114] <- 0
testcor[115:337,115:337] <- 0
length(c(testcor)[abs(c(testcor))>0.8]) #0
length(c(testcor)[abs(c(testcor))>0.5]) #0
hist(testcor)


testcor <- cor(cbind(PSD_NET, BACT_NET), method = "spearman")
dim(testcor)
hist(testcor)
length(c(testcor)[abs(c(testcor))>0.5])
length(c(testcor)[abs(c(testcor))>0.8])


# set all texture to texture and bact to bact correlations to 0
testcor[1:116,1:116] <- 0
testcor[117:1381,117:1381] <- 0
length(c(testcor)[abs(c(testcor))>0.8]) #0
length(c(testcor)[abs(c(testcor))>0.5]) #2
hist(testcor)
