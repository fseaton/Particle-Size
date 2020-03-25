## Structural equation model using output from Multifractals v3.R
## going to use lavaan.survey for this
library(lavaan.survey)
library(tidyverse)

SEM_data <- select(Ancil, REP_ID, SQUARE, PLOTYEAR, ELEVATION_M, BROAD_HABITAT_SURVEYOR, 
                   SAAR_PRECIP, SOIL_GROUP)

summary(Res_fil_mer_wideD)

SEM_data <- merge(Res_fil_mer_wideD, SEM_data, by="REP_ID", all=FALSE)

STOCH <- select(Soil, REP_ID, C_FE_NTOTAL, C_FE_PTOTAL, C_FE_POLSEN)

SEM_data <- merge(SEM_data, STOCH, by="REP_ID")

# calculate median grain size
PSD_cum <- PSD
for (i in 2:ncol(PSD_cum)){
  PSD_cum[,i] <- rowSums(PSD[,seq(1,i)])
}
Sizes <- as.numeric(substring(colnames(PSD_cum),2))
med_grain <- apply(PSD_cum, 1, function(x) min(Sizes[x>=0.5]))
hist(log(med_grain))

mg <- data.frame(REP_ID = names(med_grain),
                 med_grain = med_grain)
SEM_data <- merge(SEM_data, mg, by= "REP_ID", all = TRUE)

summary(SEM_data)
SEM_data <- filter(SEM_data, C_FE_CTOTAL > 0.02)
summary(SEM_data)

SEM_data$NC <- 10*SEM_data$C_FE_NTOTAL/SEM_data$C_FE_CTOTAL

SEM_data$CARB <- log(SEM_data$C_FE_CTOTAL)
SEM_data$MDGR <- log(SEM_data$med_grain)
SEM_data$BACT <- SEM_data$BACT_R_RICH/1000
SEM_data$FUNG <- SEM_data$FUNG_R_RICH_BLAST/100
SEM_data$ELEV <- SEM_data$ELEVATION_M/100
SEM_data$PREC <- SEM_data$SAAR_PRECIP/1000

summary(SEM_data)

# quick plot
theme_set(theme_bw() + theme(panel.grid = element_blank()))
ggplot(SEM_data, aes(x = med_grain, y = qD1)) + geom_point() + scale_x_log10() + 
  labs(x = expression("Median grain size ("*mu*"m)"), y = expression(italic(D)[1]))
ggsave("D1 vs median grain size.png", path = "Graphs/", width = 10, height =8, units = "cm")

p <- egg::ggarrange(ggplot(SEM_data, aes(x = med_grain, y = qD0)) + geom_point() + scale_x_log10() + 
                 labs(x = expression("Median grain size ("*mu*"m)"), y = expression(italic(D)[0])),
               ggplot(SEM_data, aes(x = med_grain, y = qD1)) + geom_point() + scale_x_log10() + 
                 labs(x = expression("Median grain size ("*mu*"m)"), y = expression(italic(D)[1])),
               ggplot(SEM_data, aes(x = med_grain, y = qD2)) + geom_point() + scale_x_log10() + 
                 labs(x = expression("Median grain size ("*mu*"m)"), y = expression(italic(D)[2])),
               nrow = 1)
ggsave("D012 vs median grain size.png", path = "Graphs/", 
       plot = p, width = 25, height = 10, units = "cm")

## Base model
Basemod <- '
CARB ~ PREC + ELEV + qD1
NC ~ CARB + ELEV + PREC + qD1
C_B_PH_CACL2 ~ NC + CARB + PREC + ELEV + qD1
BACT ~ NC + CARB + C_B_PH_CACL2 + qD1
FUNG ~ BACT + NC + CARB + C_B_PH_CACL2 + qD1'

BaseSEM <- sem(Basemod, data=SEM_data)
summary(BaseSEM, rsquare = TRUE)

## correction for survey design
SEM_data$SQUARE <- as.factor(SEM_data$SQUARE)
Srvy <- svydesign(ids=~SQUARE, prob=~1, data=SEM_data)
BaseSEM <- lavaan.survey(BaseSEM, Srvy)

summary(BaseSEM)
sink("SEM output.txt", split=TRUE)
print("Base SEM")
summary(BaseSEM, rsquare=TRUE, standardized = TRUE)
sink()

## remove fungi ~ pH link
Mod1 <- '
CARB ~ PREC + ELEV + qD1
NC ~ CARB + ELEV + PREC + qD1
C_B_PH_CACL2 ~ NC + CARB + PREC + ELEV + qD1
BACT ~ NC + CARB + C_B_PH_CACL2 + qD1
FUNG ~ BACT + NC + CARB + qD1'
Fit1 <- sem(Mod1, data=SEM_data)
Fit1 <- lavaan.survey(Fit1, Srvy)
summary(Fit1, standardized=TRUE, rsquare = TRUE)

## normality check ####
source(***REMOVED***)
source(***REMOVED***)
Fit1_b <- sem(Mod1, data=SEM_data, meanstructure=T)
Resid <- residuals_lavaan(Fit1_b)
head(Resid)
psych::multi.hist(Resid)

par(mfrow=c(2,3))
apply(Resid, 2, function(x){car::qqp(x)})
par(mfrow=c(1,1))
## not predicting high enough values for carbon and pH
## bacteria isn't getting low enough values

#Mardia's test
mt <- MVN::mvn(Resid)
mt

## residuals not normal

## residuals vs fitted?
Fitt <- fitted_lavaan(Fit1_b)
head(Fitt)

par(mfrow=c(2,2))
for (i in 1:ncol(Resid)){
  plot(Fitt[,i], Resid[,i], main=colnames(Fitt)[i], xlab="Fitted values", ylab="Residuals",
       col=SEM_data$Habitat)
}
par(mfrow = c(1,1))


## No nitrogen ####
Mod2 <- '
CARB ~ PREC + ELEV + qD1
C_B_PH_CACL2 ~ CARB + PREC + ELEV + qD1
BACT ~ CARB + C_B_PH_CACL2 + qD1
FUNG ~ BACT + CARB + C_B_PH_CACL2 + qD1'

Fit2 <- sem(Mod2, data=SEM_data)
summary(Fit2, rsquare = TRUE)

## correction for survey design
Fit2 <- lavaan.survey(Fit2, Srvy)

summary(Fit2)
sink("SEM output.txt", append = TRUE)
print("No nitrogen")
summary(Fit2, rsquare=TRUE, standardized=TRUE)
sink()

## remove fungi ~ pH link
Mod3 <- '
CARB ~ PREC + ELEV + qD1
C_B_PH_CACL2 ~ CARB + PREC + ELEV + qD1
BACT ~ CARB + C_B_PH_CACL2 + qD1
FUNG ~ BACT + CARB + qD1'
Fit3 <- sem(Mod3, data=SEM_data)
Fit3 <- lavaan.survey(Fit3, Srvy)

sink("SEM output.txt", append=TRUE, split=TRUE)
print("No nitrogen no fungi ~ pH")
summary(Fit3, standardized=TRUE, rsquare = TRUE)
sink()

## normality check ####
Fit1_b <- sem(Mod1, data=SEM_data, meanstructure=T)
Resid <- residuals_lavaan(Fit1_b)
head(Resid)
psych::multi.hist(Resid)

par(mfrow=c(2,3))
apply(Resid, 2, function(x){car::qqp(x)})
par(mfrow=c(1,1))
## not predicting high enough values for carbon and pH
## bacteria isn't getting low enough values

#Mardia's test
mt <- MVN::mvn(Resid)
mt

## residuals not normal

## residuals vs fitted?
Fitt <- fitted_lavaan(Fit1_b)
head(Fitt)

par(mfrow=c(2,2))
for (i in 1:ncol(Resid)){
  plot(Fitt[,i], Resid[,i], main=colnames(Fitt)[i], xlab="Fitted values", ylab="Residuals",
       col=SEM_data$Habitat)
}
par(mfrow = c(1,1))




## now fit SEM for only grassland and arable ####
SEM_grass <- filter(SEM_data, Habitat == "Arable"|Habitat == "Impr Grass"|
                      Habitat == "Neutr Grass"|Habitat == "Acid Grass")
summary(SEM_grass)


GrassSEM <- sem(Basemod, data=SEM_grass)
summary(GrassSEM, rsquare = TRUE)

## correction for survey design
Srvy_grass <- svydesign(ids=~SQUARE, prob=~1, data=SEM_grass)
GrassSEM <- lavaan.survey(GrassSEM, Srvy_grass)

summary(GrassSEM)
sink("SEM output.txt", append = TRUE)
print("Grass SEM")
summary(GrassSEM, standardized=TRUE, rsquare = TRUE)
sink()

Grass_b <- sem(Basemod, data=SEM_grass, meanstructure=T)
Resid <- residuals_lavaan(Grass_b)
head(Resid)
psych::multi.hist(Resid)

par(mfrow=c(2,2))
apply(Resid, 2, function(x){car::qqp(x)})
par(mfrow=c(1,1))
## not predicting high enough values for carbon and pH
## bacteria isn't getting low enough values

#Mardia's test
MVN::mvn(Resid)

## fotted vc residuals
Fitt <- fitted_lavaan(Grass_b)
head(Fitt)

par(mfrow=c(2,2))
for (i in 1:ncol(Resid)){
  plot(Fitt[,i], Resid[,i], main=colnames(Fitt)[i], xlab="Fitted values", ylab="Residuals",
       col=SEM_grass$Habitat)
}
par(mfrow = c(1,1))




### Binary variable for habitat intensity ####
SEM_data$HAB_INT <- 0

SEM_data$HAB_INT <- ifelse(SEM_data$Habitat == "Arable"|SEM_data$Habitat == "Improved Grass"|
                             SEM_data$Habitat == "Neutral Grass",1, 0)
table(SEM_data$HAB_INT)
plot(HAB_INT ~ Habitat, SEM_data)


Mod4 <- '
CARB ~ PREC + ELEV + qD1 + HAB_INT
C_B_PH_CACL2 ~ CARB + ELEV + qD1 + HAB_INT
BACT ~ CARB + C_B_PH_CACL2 + qD1 + HAB_INT
FUNG ~ BACT + qD1 + HAB_INT'
Fit4 <- sem(Mod4, data=SEM_data)
Srvy <- svydesign(ids=~SQUARE, prob=~1, data=SEM_data)
Fit4 <- lavaan.survey(Fit4, Srvy)

sink("SEM output.txt", append=TRUE, split=TRUE)
print("Habitat intensity No nitrogen no fungi ~ pH")
summary(Fit4, standardized=TRUE, rsquare = TRUE)
sink()

Fit4_b <- sem(Mod4, data=SEM_data, meanstructure=T)
Resid <- residuals_lavaan(Fit4_b)
head(Resid)
psych::multi.hist(Resid)

par(mfrow=c(2,2))
apply(Resid, 2, function(x){car::qqp(x)})
par(mfrow=c(1,1))
## not predicting high enough values for carbon and pH

#Mardia's test
MVN::mvn(Resid)

## fitted vs residuals
Fitt <- fitted_lavaan(Fit4_b)
head(Fitt)

par(mfrow=c(2,2))
for (i in 1:ncol(Resid)){
  plot(Fitt[,i], Resid[,i], main=colnames(Fitt)[i], xlab="Fitted values", ylab="Residuals",
       col=SEM_data$Habitat)
}
par(mfrow = c(1,1))


# with median grain size
Mod5 <- '
CARB ~ PREC + ELEV + qD1 + HAB_INT + MDGR
C_B_PH_CACL2 ~ CARB + ELEV + qD1 + HAB_INT + MDGR
BACT ~ CARB + C_B_PH_CACL2 + qD1 + HAB_INT + MDGR
FUNG ~ BACT + qD1 + HAB_INT + MDGR'
Fit5 <- sem(Mod5, data=SEM_data)
Srvy <- svydesign(ids=~SQUARE, prob=~1, data=SEM_data)
Fit5 <- lavaan.survey(Fit5, Srvy)

summary(Fit5, standardized=TRUE, rsquare = TRUE)

Fit5_b <- sem(Mod5, data=SEM_data, meanstructure=T)
Resid <- residuals_lavaan(Fit5_b)
head(Resid)
psych::multi.hist(Resid)

par(mfrow=c(2,2))
apply(Resid, 2, function(x){car::qqp(x)})
par(mfrow=c(1,1))
## not predicting high enough values for carbon and pH

#Mardia's test
MVN::mvn(Resid)

## fitted vs residuals
Fitt <- fitted_lavaan(Fit5_b)
head(Fitt)

par(mfrow=c(2,2))
for (i in 1:ncol(Resid)){
  plot(Fitt[,i], Resid[,i], main=colnames(Fitt)[i], xlab="Fitted values", ylab="Residuals",
       col=SEM_data$Habitat)
}
par(mfrow = c(1,1))

# with median grain size on carbon and fungi only
Mod6 <- '
CARB ~ PREC + ELEV + qD1 + HAB_INT + MDGR
C_B_PH_CACL2 ~ CARB + ELEV + qD1 + HAB_INT
BACT ~ CARB + C_B_PH_CACL2 + qD1 + HAB_INT
FUNG ~ BACT + qD1 + HAB_INT + MDGR'
Fit6 <- sem(Mod6, data=SEM_data)
Fit6 <- lavaan.survey(Fit6, Srvy)

summary(Fit6, standardized=TRUE, rsquare = TRUE)

Fit6_b <- sem(Mod6, data=SEM_data, meanstructure=T)
Resid <- residuals_lavaan(Fit6_b)
head(Resid)
psych::multi.hist(Resid)

par(mfrow=c(2,2))
apply(Resid, 2, function(x){car::qqp(x)})
par(mfrow=c(1,1))
## not predicting high enough values for carbon and pH

#Mardia's test
MVN::mvn(Resid)

## fitted vs residuals
Fitt <- fitted_lavaan(Fit6_b)
head(Fitt)

par(mfrow=c(2,2))
for (i in 1:ncol(Resid)){
  plot(Fitt[,i], Resid[,i], main=colnames(Fitt)[i], xlab="Fitted values", ylab="Residuals",
       col=SEM_data$Habitat)
}
par(mfrow = c(1,1))

# with median grain size not textural heterogeneity
Mod7 <- '
CARB ~ PREC + ELEV + HAB_INT + MDGR
C_B_PH_CACL2 ~ CARB + ELEV + HAB_INT + MDGR
BACT ~ CARB + C_B_PH_CACL2 + HAB_INT + MDGR
FUNG ~ BACT + HAB_INT + MDGR'
Fit7 <- sem(Mod7, data=SEM_data)
Fit7 <- lavaan.survey(Fit7, Srvy)

summary(Fit7, standardized=TRUE, rsquare = TRUE)

Fit7_b <- sem(Mod7, data=SEM_data, meanstructure=T)
Resid <- residuals_lavaan(Fit7_b)
head(Resid)
psych::multi.hist(Resid)

par(mfrow=c(2,2))
apply(Resid, 2, function(x){car::qqp(x)})
par(mfrow=c(1,1))
## not predicting high enough values for carbon and pH

#Mardia's test
MVN::mvn(Resid)

## fitted vs residuals
Fitt <- fitted_lavaan(Fit7_b)
head(Fitt)

par(mfrow=c(2,2))
for (i in 1:ncol(Resid)){
  plot(Fitt[,i], Resid[,i], main=colnames(Fitt)[i], xlab="Fitted values", ylab="Residuals",
       col=SEM_data$Habitat)
}
par(mfrow = c(1,1))


# with median grain size not textural heterogeneity
Mod8 <- '
CARB ~ PREC + ELEV + HAB_INT
C_B_PH_CACL2 ~ CARB + ELEV + HAB_INT
BACT ~ CARB + C_B_PH_CACL2 + HAB_INT
FUNG ~ BACT + HAB_INT'
Fit8 <- sem(Mod8, data=SEM_data)
Fit8 <- lavaan.survey(Fit8, Srvy)

summary(Fit8, standardized=TRUE, rsquare = TRUE)

Fit8_b <- sem(Mod8, data=SEM_data, meanstructure=T)
Resid <- residuals_lavaan(Fit8_b)
head(Resid)
psych::multi.hist(Resid)

par(mfrow=c(2,2))
apply(Resid, 2, function(x){car::qqp(x)})
par(mfrow=c(1,1))
## not predicting high enough values for carbon and pH

#Mardia's test
MVN::mvn(Resid)

## fitted vs residuals
Fitt <- fitted_lavaan(Fit8_b)
head(Fitt)

par(mfrow=c(2,2))
for (i in 1:ncol(Resid)){
  plot(Fitt[,i], Resid[,i], main=colnames(Fitt)[i], xlab="Fitted values", ylab="Residuals",
       col=SEM_data$Habitat)
}
par(mfrow = c(1,1))


aictab.lavaan(list(Fit4,Fit5,Fit6,Fit7, Fit8), c("Fit4","Fit5","Fit6","Fit7","Fit8"))

Mod9 <- '
CARB ~ PREC + ELEV + qD1 + HAB_INT
C_B_PH_CACL2 ~ CARB + ELEV + qD1 + HAB_INT
BACT ~ CARB + C_B_PH_CACL2 + qD1 + HAB_INT + PREC
FUNG ~ BACT + qD1 + HAB_INT + PREC'
Fit9 <- sem(Mod9, data=SEM_data)
Fit9 <- lavaan.survey(Fit9, Srvy)

summary(Fit9)
