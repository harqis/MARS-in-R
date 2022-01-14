#
# Tommi Kivinen
# Tilastotieteen kandidaatintyossa kaytetty R-koodi
# Tyo nahtavissa: https://urn.fi/URN:NBN:fi:tuni-202104263685
# 2021
#

#
# data
#

library(readxl)
data1 <- read_excel("Tilastotiede/Kandi/dataeturauhastutkimus.xlsx", 
                   sheet = "file04", skip = 1)
View(data1)
attach(data1)

rekisteri <- read_excel("Tilastotiede/Kandi/dataeturauhastutkimus.xlsx", 
                    sheet = "file01 regist", skip = 1)

rauhaskoot <- read_excel("Tilastotiede/Kandi/dataeturauhastutkimus.xlsx", 
                          sheet = "file05 tumour ass", skip = 1)
# Muutama joilla ei ole kokotietoja

# 1.12. esittely
#bmi
bmi <- rekisteri[,6]/((rekisteri[,5]*0.01)^2)
rekisteri <- cbind(rekisteri, bmi)

allData <- merge(rekisteri, data1, by="Patient_no")
allData <- allData[,-8]
allData <- merge(allData, rauhaskoot, by="Patient_no")
allData <- subset(allData, select=-c(Visit))
allData <- allData[-512,]

allData_no_na <- allData[complete.cases(allData),]
#allData_no_na <- as.matrix(allData_no_na)

# PSA histogrammi log-muunnoksen jalkeen
hist(log(allData_no_na[,8]))

# mars sisaltava pakkaus
library(mda)

# ei toimi, tarkista!!!
# Error in mars(allData[, -8], allData[, 8]) : 
# NA/NaN/Inf in foreign function call (arg 6)
# In addition: Warning messages:
# 1: In storage.mode(tagx) <- "integer" : NAs introduced by coercion
# 2: In storage.mode(x) <- "double" : NAs introduced by coercion

#potilas nro 21015 on ongelma, yksi 'NA'

mars(allData_no_na[,-8], allData_no_na[,8])

# selittaja matriisi

selittajat <- allData[,7:17]
selittajat <- cbind(allData[,7], selittajat)
selittajat <- cbind(selittajat, allData[,20:21])
selittajat <- selittajat[complete.cases(selittajat),]
selittajat <- selittajat[, -1]

# selitettava eli log(PSA)
selitettava <- selittajat[, 2]

selitettava1 <- log(selitettava1) #kaikki ei na

selittajat <- selittajat[, -2]
#numeerinen matriisi, kayta tata marsiin
selittajat1 <- as.matrix(selittajat)

#
# testeja
#

sapply(allData_no_na, is.numeric)

hist(PSA) # hyvin vino

library(MASS)
testlm1 <- lm(PSA ~ 1, data=data1)

boxcox(testlm1) # suosittelee lambda=0 eli logaritmimuunnosta

hist(log(PSA))

logPSA <- log(PSA)

#
# mars
#

library(mda)

a1 <- mars(selittajat1, selitettava1)

a1$selected.terms

a1$factor[a1$s,]

showcuts <- function(obj)
{
  tmp <- obj$cuts[obj$sel, ]
  dimnames(tmp) <- list(NULL, names(selittajat1))
  tmp
}

showcuts(a1)


## examine the fitted functions
par(mfrow=c(3,4), pty="s")
Xp <- matrix(sapply(selittajat1[1:12], mean), nrow(selittajat1), ncol(selittajat1), byrow=TRUE)
for(i in 1:12) {
  xr <- sapply(selittajat1, range)
  Xp[,i] <- seq(xr[1,i], xr[2,i], len=nrow(selittajat1))
  Xf <- predict(a1, Xp)
  plot(Xp[ ,i], Xf, xlab=names(selittajat1)[i], ylab="", type="l")
}

#
# Uusi yritys marsilla, talla kertaa muuttujat samassa matriisissa
#

aineisto <- cbind(selittajat1, selitettava1)

fit1 <- mars(aineisto[,1:12], aineisto[,13])

## tutkitaan fitattuja funktioita
par(mfrow=c(3,4),mar=c(2, 2, 2, 2), pty="s")
for(i in 1:12) 
{
  Xp <- matrix(sapply(aineisto[,1:12], mean), nrow(aineisto), ncol(aineisto) - 1, byrow=TRUE)
  xr <- sapply(aineisto, range)
  Xp[,i] <- seq(xr[1,i], xr[2,i], len=nrow(aineisto))
  Xf <- predict(fit2, Xp)
  plot(Xp[ ,i], Xf, xlab=names(aineisto)[,i], ylab="", type="l")
}

ages <- allData[,3]-allData[,4]

#
# Uusi yritys earth-pakkauksella
# Lopulliseen tyohon kaytetty koodi 
#

#earth-pakkaus
library(earth)

colnames(aineisto) <- c("Ika", "BMI", "Testosteroni", "AFOS", "Kreatiniini", 
                        "Hematokriitti", "Hemoglobiini", "Valkosolut", "Punasolut",
                        "Korpuskulaarinen_tilavuus", "Korpuskulaarinen_hemoglobiini",
                        "Leveys", "Pituus", "logPSA")

fit3 <- earth(aineisto[,1:13], aineisto[,14], trace=2)
fit3.1 <- earth(aineisto[,1:13], aineisto[,14], nfold = 2, ncross = 30, varmod.method = "lm")
fit4 <- earth(aineisto[,1:13], aineisto[,14], degree=2)
fit4.1 <- earth(aineisto[,1:13], aineisto[,14], degree=2, nfold = 2, ncross = 30, varmod.method = "lm")

fit5 <- earth(aineisto[,1:13], aineisto[,14], degree=3)
fit5.1 <- earth(aineisto[,1:13], aineisto[,14], degree=3, nfold = 2, ncross = 30, varmod.method = "lm")
fit6 <- earth(aineisto[,1:13], aineisto[,14], degree=4) # paras selitysaste, 0.33
fit6.1 <- earth(aineisto[,1:13], aineisto[,14], degree=4, nfold = 2, ncross = 30, varmod.method = "lm")
fit7 <- earth(aineisto[,1:13], aineisto[,14], degree=5) #ei parane, samat kuin fit6

fit.lm1 <- lm(aineisto[,13]~aineisto[,3])

plot(fit3)
malli1_summary <- summary(fit3)
format(fit3)
evimp(fit3)
plotmo(fit3)
plotmo(fit3, degree1 = c(2), caption = "", xlab = "AFOS (U/L)", ylab = "log(PSA)")
plotmo(fit3, degree1 = c(5), caption = "", main = "Eturauhasen pituussuuntainen koko",
       xlab = "Pituus (cm)", ylab = "log(PSA)")

plot(fit4)
summary(fit4)
format(fit4)
plotmo(fit4)
plotmo(fit4, degree1 = c(), degree2 = c(1), caption = "", ylab = "PSA")
evimp(fit4)

plot(fit5)
summary(fit5)
plotmo(fit5)
plotmo(fit5, degree1 = c(), degree2 = c(1), caption = "", ylab = "PSA")
evimp(fit5)

summary(fit6)
plotmo(fit6)

apply(aineisto, 2, mean)
apply(aineisto, 2, sd)

#
#tavallinen lineaarinen regressio vertailuksi
#

vertailu1 <- lm(logPSA ~ bmi+Testosterone+Alkaline+Creatinine+Haematocrit+age+Haemoglobin+White_cells+Erythrocytes+Corpuscular_volume+Corpuscular_haemaglobin+Width+Length, data=aineisto)
vertailu2 <- lm(logPSA ~ bmi+Alkaline+age+Length, data=aineisto)
vertailu3 <- lm(logPSA ~ bmi+Alkaline+Length, data=aineisto)

summary(vertailu2)

write.table(malli1_summary, file = "malli1_summary.txt", sep = ",", quote = FALSE, row.names = F)
