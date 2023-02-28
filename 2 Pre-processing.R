#### Pre-processing ####

library(caret)
library(readxl)
library(tidyr)

bp1 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/BloodPressure_First_Second_Readings.xlsx", sheet = 1)
bp2 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/BloodPressure_First_Second_Readings.xlsx", sheet = 2)
cau <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/CardiacUltrasound_CIMT.xlsx", sheet = 1)
ch1 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/Cholesterol_2007-2018_NO_PCODES.xlsx", sheet = 1)
ch2 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/Cholesterol_2007-2018_NO_PCODES.xlsx", sheet = 2)
dx1 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_Batch2_23062021/DEXA.xlsx", sheet = 1)
dx2 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_Batch2_23062021/DEXA.xlsx", sheet = 2)
dx3 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_Batch2_23062021/DEXA.xlsx", sheet = 3)
dx4 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_Batch2_23062021/DEXA.xlsx", sheet = 4)
dx5 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_Batch2_23062021/DEXA.xlsx", sheet = 5)
dx6 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_Batch2_23062021/DEXA.xlsx", sheet = 6)
dx7 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_Batch2_23062021/DEXA.xlsx", sheet = 7)
dx8 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_Batch2_23062021/DEXA.xlsx", sheet = 8)
dx9 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_Batch2_23062021/DEXA.xlsx", sheet = 9)
dxx <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_Batch2_23062021/DEXA.xlsx", sheet = 10)
fih <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Other/DS00027_Frailty_Index_HATS/Frailty_Index.xlsx", sheet = 1)
gi1 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/Glucose_Insulin_2007-2018_NO_PCODES.xlsx", sheet = 1)
gi2 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/Glucose_Insulin_2007-2018_NO_PCODES.xlsx", sheet = 2)
lf1 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/LungFunction_Spirometry_Vitalograph.xlsx", sheet = 1)
lf2 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/LungFunction_Spirometry_Vitalograph.xlsx", sheet = 2)
mh1 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/Measured_HeightWeight.xlsx", sheet = 1)
mh2 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/Measured_HeightWeight.xlsx", sheet = 2)
p11 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/Phenotypes_Batch1.xlsx", sheet = 1)
p12 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/Phenotypes_Batch1.xlsx", sheet = 2)
p21 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/Phenotypes_Batch2.xlsx", sheet = 1)
p22 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_Batch1_28052021/Phenotypes_Batch2.xlsx", sheet = 2)

colnames(x = gi1)[2:3] <- c("VisitDate", "Glucose plasma")
colnames(x = ch1)[1:6] <- c("ParticipantID", "VisitDate", "Total cholesterol", "High density lipoprotein", "Low density lipoprotein", "Direct LDL cholesterol")
colnames(x = bp1)[3:6] <- c("Systolic blood pressure second measurement", "Diastolic blood pressure second measurement", "Systolic blood pressure third measurement", "Diastolic blood pressure third measurement")
colnames(x = lf1)[4:8] <- c("Number of blows during the lung function assessment", "Lung function predicted forced expiratory volume at 1 second", "Lung function measured forced expiratory volume at 1 second", "Lung function measured forced vital capacity", "Room temperature during the lung function assessment")
colnames(x = dx1)[3:ncol(x = dx1)] <- dx1[1, 3:ncol(x = dx1)]
colnames(x = dx2)[3:ncol(x = dx2)] <- dx2[1, 3:ncol(x = dx2)]
colnames(x = dx3)[3:ncol(x = dx3)] <- dx3[1, 3:ncol(x = dx3)]
colnames(x = dx4)[3:ncol(x = dx4)] <- dx4[1, 3:ncol(x = dx4)]
colnames(x = dx5)[3:ncol(x = dx5)] <- dx5[1, 3:ncol(x = dx5)]
colnames(x = dx6)[3:ncol(x = dx6)] <- dx6[1, 3:ncol(x = dx6)]
colnames(x = dx7)[3:ncol(x = dx7)] <- dx7[1, 3:ncol(x = dx7)]
colnames(x = dx8)[3:ncol(x = dx8)] <- dx8[1, 3:ncol(x = dx8)]
colnames(x = dx9)[3:ncol(x = dx9)] <- dx9[1, 3:ncol(x = dx9)]

dx1 <- dx1[-1, ]
dx2 <- dx2[-1, ]
dx3 <- dx3[-1, ]
dx4 <- dx4[-1, ]
dx5 <- dx5[-1, ]
dx6 <- dx6[-1, ]
dx7 <- dx7[-1, ]
dx8 <- dx8[-1, ]
dx9 <- dx9[-1, ]

dx1[3:ncol(x = dx1)] <- sapply(X = dx1[3:ncol(x = dx1)], FUN = as.numeric)
dx2[3:ncol(x = dx2)] <- sapply(X = dx2[3:ncol(x = dx2)], FUN = as.numeric)
dx3[3:ncol(x = dx3)] <- sapply(X = dx3[3:ncol(x = dx3)], FUN = as.numeric)
#dx4[3:ncol(x = dx4)] <- sapply(X = dx4[3:ncol(x = dx4)], FUN = as.numeric)
dx5[3:ncol(x = dx5)] <- sapply(X = dx5[3:ncol(x = dx5)], FUN = as.numeric)
dx6[3:ncol(x = dx6)] <- sapply(X = dx6[3:ncol(x = dx6)], FUN = as.numeric)
dx7[3:ncol(x = dx7)] <- sapply(X = dx7[3:ncol(x = dx7)], FUN = as.numeric)
dx8[3:ncol(x = dx8)] <- sapply(X = dx8[3:ncol(x = dx8)], FUN = as.numeric)
dx9[3:ncol(x = dx9)] <- sapply(X = dx9[3:ncol(x = dx9)], FUN = as.numeric)

dx1$`Percentage of fat in largest visceral fat region / Percentage of fat in android region` <- dx1$`Percentage of fat in largest visceral fat region` / dx1$`Percentage of fat in android region` 
dx1$`Percentage of fat in largest visceral fat region / Percentage of fat in gynoid region` <- dx1$`Percentage of fat in largest visceral fat region` / dx1$`Percentage of fat in gynoid region`

dx1 <- predict(object = preProcess(x = dx1, method = c("nzv", "corr"), na.remove = TRUE), newdata = dx1)
dx2 <- predict(object = preProcess(x = dx2, method = c("nzv", "corr"), na.remove = TRUE), newdata = dx2)
dx3 <- predict(object = preProcess(x = dx3, method = c("nzv", "corr"), na.remove = TRUE), newdata = dx3)
#dx4 <- predict(object = preProcess(x = dx4, method = c("nzv", "corr"), na.remove = TRUE), newdata = dx4)
dx5 <- predict(object = preProcess(x = dx5, method = c("nzv", "corr"), na.remove = TRUE), newdata = dx5)
dx6 <- predict(object = preProcess(x = dx6, method = c("nzv", "corr"), na.remove = TRUE), newdata = dx6)
dx7 <- predict(object = preProcess(x = dx7, method = c("nzv", "corr"), na.remove = TRUE), newdata = dx7)
dx8 <- predict(object = preProcess(x = dx8, method = c("nzv", "corr"), na.remove = TRUE), newdata = dx8)
dx9 <- predict(object = preProcess(x = dx9, method = c("nzv", "corr"), na.remove = TRUE), newdata = dx9)

gi1$Insulin <- as.numeric(x = gi1$Insulin)

mhw <- merge(x = mh1, y = mh2, by = c(1:2), sort = FALSE)
mhw$`Body Mass Index` <- round(x = mhw$Weight / (mhw$Height / 100)^2, digits = 1)

lf1 <- pivot_longer(data = lf1[-1, -3], cols = 3:7, names_to = "Description", values_to = "Result")
ch1 <- pivot_longer(data = ch1, cols = 3:ncol(x = ch1), names_to = "Description", values_to = "Result")
gi1 <- pivot_longer(data = gi1, cols = 3:ncol(x = gi1), names_to = "Description", values_to = "Result")
mhw <- pivot_longer(data = mhw, cols = 3:ncol(x = mhw), names_to = "Description", values_to = "Result")
dx1 <- pivot_longer(data = dx1, cols = 3:ncol(x = dx1), names_to = "Description", values_to = "Result")
dx2 <- pivot_longer(data = dx2, cols = 3:ncol(x = dx2), names_to = "Description", values_to = "Result")
dx3 <- pivot_longer(data = dx3, cols = 3:ncol(x = dx3), names_to = "Description", values_to = "Result")
#dx4 <- pivot_longer(data = dx4, cols = 3:ncol(x = dx4), names_to = "Description", values_to = "Result")
dx5 <- pivot_longer(data = dx5, cols = 3:ncol(x = dx5), names_to = "Description", values_to = "Result")
dx6 <- pivot_longer(data = dx6, cols = 3:ncol(x = dx6), names_to = "Description", values_to = "Result")
dx7 <- pivot_longer(data = dx7, cols = 3:ncol(x = dx7), names_to = "Description", values_to = "Result")
dx8 <- pivot_longer(data = dx8, cols = 3:ncol(x = dx8), names_to = "Description", values_to = "Result")
dx9 <- pivot_longer(data = dx9, cols = 3:ncol(x = dx9), names_to = "Description", values_to = "Result")
bp1 <- pivot_longer(data = bp1[-1, ], cols = 3:ncol(x = bp1), names_to = "Description", values_to = "Result")

bp2 <- bp2[, c("ParticipantID", "VisitDate", "Description", "Result")]
cau <- cau[, c("ParticipantID", "VisitDate", "Description", "Result")]
ch2 <- ch2[, c("ParticipantID", "VisitDate", "Description", "Result")]
dxx <- dxx[, c("ParticipantID", "VisitDate", "Description", "Result")]
gi2 <- gi2[, c("ParticipantID", "VisitDate", "Description", "Result")]
lf2 <- lf2[, c("ParticipantID", "VisitDate", "Description", "Result")]
p11 <- p11[, c("ParticipantID", "VisitDate", "Description", "Result")]
p12 <- p12[, c("ParticipantID", "VisitDate", "Description", "Result")]
p21 <- p21[, c("ParticipantID", "VisitDate", "Description", "Result")]
p22 <- p22[, c("ParticipantID", "VisitDate", "Description", "Result")]

bp1$Result <- as.numeric(x = bp1$Result)
bpa <- rbind.data.frame(bp1, bp2, stringsAsFactors = FALSE)
bpa <- split(x = bpa, f = bpa$Description)

sbp2 <- bpa$`Systolic blood pressure second measurement`
sbp3 <- bpa$`Systolic blood pressure third measurement`
dbp2 <- bpa$`Diastolic blood pressure second measurement`
dbp3 <- bpa$`Diastolic blood pressure third measurement`

sbp <- merge(x = sbp2, y = sbp3, by = c(1:2), all = TRUE, sort = FALSE)
sbp$Result <- rowMeans(x = sbp[, c(4, 6)], na.rm = TRUE)
sbp$Description <- "Systolic blood pressure"
sbp <- sbp[, c("ParticipantID", "VisitDate", "Description", "Result")]

dbp <- merge(x = dbp2, y = dbp3, by = c(1:2), all = TRUE, sort = FALSE)
dbp$Result <- rowMeans(x = dbp[, c(4, 6)], na.rm = TRUE)
dbp$Description <- "Diastolic blood pressure"
dbp <- dbp[, c("ParticipantID", "VisitDate", "Description", "Result")]

pp <- merge(x = sbp, y = dbp, by = c(1:2), all = TRUE, sort = FALSE)
pp$Description <- "Pulse pressure"
pp$Result <- pp$Result.x - pp$Result.y
pp <- pp[, c("ParticipantID", "VisitDate", "Description", "Result")]

dxa <- rbind.data.frame(dx1, dx5, dx9, stringsAsFactors = FALSE)
dxa$VisitDate <- as.Date.POSIXct(x = dxa$VisitDate, format = "%Y/%m/%d")

fih$Description <- "Frailty Index"

colnames(x = fih)[4] <- "Result"

fih <- fih[, c("ParticipantID", "VisitDate", "Description", "Result")]

g_i <- rbind.data.frame(gi1, gi2, stringsAsFactors = FALSE)

fasting <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_04022021/FastingStatus_LastEaten.xlsx", sheet = 3)
fasting$ResponseDate <- as.Date(x = fasting$ResponseDate, format = "%d/%m/%Y")

g_i <- merge(x = g_i, y = fasting[, c(3:4, 6)], by.x = c("ParticipantID", "VisitDate"), by.y = c("ParticipantID", "ResponseDate"), all = TRUE, sort = FALSE)
g_i <- g_i[g_i$ResponseDescription == "Fasted", -5]

lf1$Result <- as.numeric(x = lf1$Result)

wc <- p11[p11$Description == "Measured hip circumference", ]
hc <- p11[p11$Description == "Measured waist circumference", ]

whr <- merge(x = wc, y = hc, by = c(1:2), all = TRUE, sort = FALSE)
whr$Description <- "Waist-to-hip ratio"
whr$Result <- whr$Result.y / whr$Result.x
whr <- whr[, c("ParticipantID", "VisitDate", "Description", "Result")]

mta <- rbind.data.frame(sbp, dbp, pp, cau, ch1, ch2, dxa, fih, g_i, lf1, lf2, mhw, p11, whr, p12, p21, p22, stringsAsFactors = FALSE)

#### ASCVD ####

library(dplyr)
library(ggplot2)
library(ggsci)
library(psych)
library(reshape2)

subject <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_22032021.xlsx", sheet = 1)

phenos <- merge(x = mta, y = subject[, c("STUDY_NO", "YEAR_BIRTH", "SEX", "ETHNICITY")], by = 1, all = TRUE)
phenos$AGE <- as.numeric(x = format(x = as.Date(phenos$VisitDate, format = "%Y-%m-%d"),"%Y")) - phenos$YEAR_BIRTH
phenos <- phenos[!is.na(x = phenos$AGE), ]

age <- phenos[, c("ParticipantID", "VisitDate", "SEX", "ETHNICITY", "AGE")][!duplicated(x = phenos[, c("ParticipantID", "VisitDate", "SEX", "ETHNICITY", "AGE")]), ]

tch <- phenos[phenos$Description == "Total cholesterol", c("ParticipantID", "VisitDate", "Result")]
tch[tch == 999] <- NA

hdl <- phenos[phenos$Description == "High density lipoprotein", c("ParticipantID", "VisitDate", "Result")]
hdl[hdl == 999 | hdl == 777] <- NA

sbp[sbp == 999] <- NA

htn <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Data/Medication/Medication_raw.csv")
htn$`Hypertension.Drugs` <- ifelse(test = rowSums(x = htn[, c("Alpha.Adrenoceptor.Blocking.Drugs", "Angiotensin.Converting.Enzyme.Inhibitors", "Angiotensin.II.Receptor.Antagonists", "Beta.Adrenoceptor.Blocking.Drugs", "Calcium.Channel.Blockers", "Centrally.Acting.Antihypertensive.Drugs", "Loop.Diuretics", "Peripheral.Vasodilators...Related.Drugs", "Pot.Sparing.Diuretics.Aldosterone.Antag", "Potassium.Sparing.Diuretics...Compounds", "Thiazides.And.Related.Diuretics")]) > 0, yes = 1, no = 0)

smokestatus <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Other/Smoking/alltwins_smokestatus_2020.csv", header = TRUE, check.names = FALSE)
smokestatus <- melt(data = smokestatus, id.vars = "StudyNo")
smokestatus[which(x = smokestatus$value > 2), ]$value <- 0
smokestatus$value <- ifelse(test = smokestatus$value == 1, yes = 1, no = 0)

diab1 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_28022022.xlsx", sheet = 1)
diab2 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_28022022.xlsx", sheet = 2)
diab3 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/Max_28022022.xlsx", sheet = 3)
diab <- as.data.frame(x = mapply(FUN = c, diab1[-1, c(1, 3)], diab2[-1, c(2, 5)], diab2[-1, c(2, 6)], diab3[-1, c(2, 5)], diab3[-1, c(2, 6)]))
diab[diab >= 999905] <- NA
diab <- na.omit(object = diab)

phenos <- Reduce(function(x, y, ...)  merge(x = x, y = y, by = c(1, 2), all = TRUE, ...), list(age, tch, hdl, sbp[, c("ParticipantID", "VisitDate", "Result")], htn[, c("ParticipantID", "Date", "Hypertension.Drugs")]))
phenos <- phenos[!is.na(x = phenos$ParticipantID), ]

for (i in unique(x = phenos$ParticipantID)) {
  
  for (j in 1:nrow(x = phenos[which(x = phenos$ParticipantID == i), ])) {
    
    if (is.na(x = phenos[which(x = phenos$ParticipantID == i), ]$Hypertension.Drugs[1])) {
      
      phenos[which(x = phenos$ParticipantID == i), ]$Hypertension.Drugs[1] <- first(x = na.omit(object = phenos[which(x = phenos$ParticipantID == i), ]$Hypertension.Drugs))

    } else if (is.na(x = phenos[which(x = phenos$ParticipantID == i), ]$Hypertension.Drugs[j])) {
      
      diff <- NULL
      
      for (k in 1:nrow(x = phenos[which(x = phenos$ParticipantID == i), ])) {
        
        diff <- c(diff, abs(x = difftime(time1 = phenos[which(x = phenos$ParticipantID == i), ][j, 2], time2 = phenos[which(x = phenos$ParticipantID == i), ][k, 2]))) }
      
      k = which.min(x = diff[which(x = phenos[which(x = phenos$ParticipantID == i), ]$Hypertension.Drugs >= 0)])
      
      phenos[which(x = phenos$ParticipantID == i), ]$Hypertension.Drugs[j] <- phenos[which(x = phenos$ParticipantID == i), ][which(x = phenos[which(x = phenos$ParticipantID == i), ]$Hypertension.Drugs >= 0), ]$Hypertension.Drugs[k] } } }

phenos$variable <- as.numeric(x = format(x = phenos$VisitDate, format = "%Y"))
phenos <- merge(x = phenos, y = smokestatus, by.x = c(1, 10), by.y = c(1:2), all.x = TRUE)
phenos <- merge(x = phenos, y = diab, by = 1, all.x = TRUE)
phenos$Diabetic <- ifelse(test = phenos$Q16_67 < phenos$AGE, yes = 1, no = 0)
phenos[is.na(phenos$Diabetic), ]$Diabetic <- 0
phenos <- na.omit(object = phenos[, -12])

colnames(x = phenos)[c(2, 7:9, 11)] <- c("Year", "TotalCholesterol", "HDLcholesterol", "SBP", "Smoker10")

women <- phenos[phenos$SEX == "F", ]
women_white <- women[women$ETHNICITY == "White", ]

MeanXB = -29.18

Age = women_white$AGE
TotChol = women_white$TotalCholesterol * 38.67
HDL = women_white$HDLcholesterol * 38.67
Treated = women_white$Hypertension.Drugs
Untreated = abs(x = women_white$Hypertension.Drugs - 1)
SBP = women_white$SBP
CurrSmoker = women_white$Smoker10
Diabetes = women_white$Diabetic

IndXB = (-29.799 * log(x = Age)) + (4.884 * log(x = Age) * log(x = Age)) + (13.540 * log(x = TotChol)) - (3.114 * log(x = Age) * log(x = TotChol)) - (13.578 * log(x = HDL)) + (3.149 * log(x = Age) * log(x = HDL)) + (Treated * 2.019 * log(x = SBP)) + (Untreated * 1.957 * log(x = SBP)) + (7.574 * CurrSmoker) - (1.665 * log(x = Age) * CurrSmoker) + (0.661 * Diabetes)

S10 = 0.9665

women_white$ASCVD <- 1-S10^exp(x = IndXB - MeanXB)
women_black <- women[women$ETHNICITY == "Black", ]

MeanXB = 86.61

Age = women_black$AGE
TotChol = women_black$TotalCholesterol * 38.67
HDL = women_black$HDLcholesterol * 38.67
Treated = women_black$Hypertension.Drugs
Untreated = abs(x = women_black$Hypertension.Drugs - 1)
SBP = women_black$SBP
CurrSmoker = women_black$Smoker10
Diabetes = women_black$Diabetic

IndXB = (17.114 * log(x = Age)) + (0.940 * log(x = TotChol)) - (18.920 * log(x = HDL)) + (4.475 * log(x = Age) * log(x = HDL)) + (Treated * 29.291 * log(x = SBP)) - (6.432 * log(x = Age) * (Treated * log(x = SBP))) + (Untreated * 27.820 * log(x = SBP)) - (6.087 * log(x = Age) * (Untreated * log(x = SBP))) + (0.691 * CurrSmoker) + (0.874 * Diabetes)

S10 = 0.9533

women_black$ASCVD <- 1-S10^exp(x = IndXB - MeanXB)

men <- phenos[phenos$SEX == "M", ]
men_white <- men[men$ETHNICITY == "White", ]

MeanXB = 61.18

Age = men_white$AGE
TotChol = men_white$TotalCholesterol * 38.67
HDL = men_white$HDLcholesterol * 38.67
Treated = men_white$Hypertension.Drugs
Untreated = abs(x = men_white$Hypertension.Drugs - 1)
SBP = men_white$SBP
CurrSmoker = men_white$Smoker10
Diabetes = men_white$Diabetic

IndXB = (12.344 * log(x = Age)) + (11.853 * log(x = TotChol)) - (2.664 * log(x = Age) * log(x = TotChol)) - (7.990 * log(x = HDL)) + (1.769 * log(x = Age) * log(x = HDL)) + (Treated * 1.797 * log(x = SBP)) + (Untreated * 1.764 * log(x = SBP)) + (7.837 * CurrSmoker) - (1.795 * log(x = Age) * CurrSmoker) + (0.658 * Diabetes)

S10 = 0.9144

men_white$ASCVD <- 1-S10^exp(x = IndXB - MeanXB)
men_black <- men[men$ETHNICITY == "Black", ]

MeanXB = 19.54

Age = men_black$AGE
TotChol = men_black$TotalCholesterol * 38.67
HDL = men_black$HDLcholesterol * 38.67
Treated = men_black$Hypertension.Drugs
Untreated = abs(x = men_black$Hypertension.Drugs - 1)
SBP = men_black$SBP
CurrSmoker = men_black$Smoker10
Diabetes = men_black$Diabetic

IndXB = (2.469 * log(x = Age)) + (0.302 * log(x = TotChol)) - (0.307 * log(x = HDL)) + (Treated * 1.916 * log(x = SBP)) + (Untreated * 1.809 * log(x = SBP)) + (0.549 * CurrSmoker) + (0.645 * Diabetes)

S10 = 0.8954

men_black$ASCVD <- 1-S10^exp(x = IndXB - MeanXB)

ascvd <- rbind(women_white, women_black, men_white, men_black)
ascvd_2 <- pivot_longer(data = ascvd, cols = 13, names_to = "Description", values_to = "Result")[, c("ParticipantID", "VisitDate", "Description", "Result")]
ascvd_f <- ascvd[ascvd$SEX == "F", ]
ascvd_m <- ascvd[ascvd$SEX == "M", ]

describe_f <- paste0(round(x = describe(x = ascvd_f[, -c(1:5, 10:12)])[, 3],  digits = 2), "-", round(x = describe(x = ascvd_f[, -c(1:5, 10:12)])[, 4],  digits = 2))
describe_f <- c(nrow(x = ascvd_f), describe_f[1:4], colSums(x = ascvd_f[, 10:12]), describe_f[5])
describe_m <- paste0(round(x = describe(x = ascvd_m[, -c(1:5, 10:12)])[, 3],  digits = 2), "-", round(x = describe(x = ascvd_m[, -c(1:5, 10:12)])[, 4],  digits = 2))
describe_m <- c(nrow(x = ascvd_m), describe_m[1:4], colSums(x = ascvd_m[, 10:12]), describe_m[5])
describe <- rbind(describe_f, describe_m)

rownames(x = describe) <- c("Females", "Males")
colnames(x = describe) <- c("N", "Age", "Total cholesterol", "High density lipoprotein", "Systolic blood pressure", "Hypertension drugs", "Smoke status", "Diabetes", "ASCVD risk score")

write.table(x = t(x = describe), file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/2_ASCVD_risk_score.csv", quote = FALSE, sep = ",", row.names = TRUE, col.names = NA)

ascvd <- ascvd %>% filter(Year >= 2001 & Year <= 2010) %>% group_by(ParticipantID) %>% arrange(ParticipantID, VisitDate) %>% mutate(VisitOrder = row_number()) %>% filter(VisitOrder == 1)

for (i in 1:10) {
  
  outliers <- ascvd$ASCVD %in% boxplot.stats(x = ascvd$ASCVD)$out
  
  ascvd[which(x = outliers), "ASCVD"] <- NA
  ascvd <- na.omit(object = ascvd) }

f <- ggplot(data = ascvd, mapping = aes(x = ASCVD*100, fill = SEX)) + geom_density(alpha = 0.5, colour = "black") + labs(x = "ASCVD risk score (%)", y = "Density") + theme_classic(base_size = 18) + scale_fill_npg()

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/2_ASCVD_risk_score.pdf", plot = f, width = 11, height = 8.5)

mta <- rbind(mta, ascvd_2)
mta <- mta[!duplicated(x = mta), ]

write.table(x = mta[with(data = mta, expr = order(Description, ParticipantID)), ], file = "~/Documents/KCL/PhD/TwinsUK/Output/Master_table.txt", quote = FALSE, sep = "\t", row.names = FALSE)
