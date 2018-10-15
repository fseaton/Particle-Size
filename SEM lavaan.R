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

summary(SEM_data)
SEM_data <- filter(SEM_data, C_FE_CTOTAL > 0.02)
summary(SEM_data)

SEM_data$NC <- 10*SEM_data$C_FE_NTOTAL/SEM_data$C_FE_CTOTAL

SEM_data$CARB <- log(SEM_data$C_FE_CTOTAL)
SEM_data$BACT <- SEM_data$BACT_R_RICH/1000
SEM_data$FUNG <- SEM_data$FUNG_R_RICH_BLAST/100
SEM_data$ELEV <- SEM_data$ELEVATION_M/100
SEM_data$PREC <- SEM_data$SAAR_PRECIP/1000

summary(SEM_data)

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

SEM_data$HAB_INT <- ifelse(SEM_data$Habitat == "Arable"|SEM_data$Habitat == "Impr Grass"|
                             SEM_data$Habitat == "Neutr Grass",1, 0)
summary(SEM_data$HAB_INT)
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

