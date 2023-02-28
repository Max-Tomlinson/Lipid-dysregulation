#### Co-variate analysis ####

setwd(dir = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Covariates")

library(lme4)
library(lubridate)
library(pbnm)
library(readxl)

#### Seasonality ####

days_in_year <- function(year) { 365 + (year %% 4 == 0) - (year %% 100 == 0) + (year %% 400 == 0) } 

covars <- read.table(file = "Adipose_TechCovars_19052017.txt", header = TRUE, quote = "", stringsAsFactors = FALSE)

date <- dmy(format(x = as.Date(x = covars$date, format = "d%y%m%d"), "%d-%m-%Y"))
  
days <- yday(x = date)
  
year <- as.numeric(x = gsub(pattern = "-.*", replacement = "", x = date))
  
prop <- days/days_in_year(year = year)
  
covars$CosSeasonality <- cos(x = 2 * pi * prop)
covars$SinSeasonality <- sin(x = 2 * pi * prop)
covars$DayofYear <- days
  
write.table(x = covars, file = "Adipose_TechCovars_19052017.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

#### Subject-level ####

classes <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Clinical Outcomes/Lipid_regulating_drugs.csv", na = "NA")

subject <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_22032021.xlsx", sheet = 1)

summary_table <- NULL

for (i in 7:8) {
  
  class_i <- classes[!is.na(x = classes[, i]), ]
  
  covars_i <- merge(x = class_i, y = subject, by.x = "ParticipantID", by.y = "STUDY_NO")
  covars_i[, c("YEAR_BIRTH")] <- scale(x = covars_i[, c("YEAR_BIRTH")])
  covars_i[, i] <- as.factor(x = covars_i[, i])
  
  colnames(x = covars_i)[i] <- "Class"
  
  covars_i_ethnicity <- covars_i[!is.na(x = covars_i[, "ETHNICITY"]), ]
  
  fix <- rbind("Year Birth" = summary(object = glm(formula = Class ~ YEAR_BIRTH, family = "binomial", data = covars_i))$coefficients[8],
                      "Sex" = summary(object = glm(formula = Class ~ factor(x = SEX), family = "binomial", data = covars_i))$coefficients[8],
                 "Zygosity" = summary(object = glm(formula = Class ~ factor(x = ACTUAL_ZYGOSITY), family = "binomial", data = covars_i))$coefficients[11])
  fix <- signif(x = fix, digits = 3)
  
  fmX01 <- glm(formula = Class ~ YEAR_BIRTH, data = covars_i, family = "binomial")
  fmX02 <- glm(formula = Class ~ YEAR_BIRTH, data = covars_i_ethnicity, family = "binomial")
  
  fmX1 <- glmer(formula = Class ~ YEAR_BIRTH + (1|Family_No), data = covars_i, family = "binomial")
  fmX2 <- glmer(formula = Class ~ YEAR_BIRTH + (1|Type_Desc), data = covars_i, family = "binomial")
  fmX3 <- glmer(formula = Class ~ YEAR_BIRTH + (1|ETHNICITY), data = covars_i_ethnicity, family = "binomial")
  
  pbnm_X1 <- pbnm(m1 = fmX1, m0 = fmX01, nsim = 5000, tasks = 10, cores = 2, seed = 3400835)
  pbnm_X2 <- pbnm(m1 = fmX2, m0 = fmX01, nsim = 5000, tasks = 10, cores = 2, seed = 3400835)
  pbnm_X3 <- pbnm(m1 = fmX3, m0 = fmX02, nsim = 5000, tasks = 10, cores = 2, seed = 3400835)

  p1 <- ifelse(test = is.na(x = summary(object = pbnm_X1)$`P(>=obs)`), yes = mean(x = c(summary(object = pbnm_X1)$`#(>=obs)/runs`, summary(object = pbnm_X1)$`#(>=obs)+noEst/runs`)), no = summary(object = pbnm_X1)$`P(>=obs)`)
  p2 <- ifelse(test = is.na(x = summary(object = pbnm_X2)$`P(>=obs)`), yes = mean(x = c(summary(object = pbnm_X2)$`#(>=obs)/runs`, summary(object = pbnm_X2)$`#(>=obs)+noEst/runs`)), no = summary(object = pbnm_X2)$`P(>=obs)`)
  p3 <- ifelse(test = is.na(x = summary(object = pbnm_X3)$`P(>=obs)`), yes = mean(x = c(summary(object = pbnm_X3)$`#(>=obs)/runs`, summary(object = pbnm_X3)$`#(>=obs)+noEst/runs`)), no = summary(object = pbnm_X3)$`P(>=obs)`)
  
  ran <- rbind("Family" = p1, "Type" = p2, "Ethnicity" = p3)
  
  column <- rbind(fix, ran)
  
  colnames(x = column) <- paste0(gsub(pattern = "_", replacement = " ", x = colnames(x = class_i)[i]), " Pr(>|z|)")
  
  summary_table <- cbind(summary_table, column) }

write.csv(x = summary_table, file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/3_Co-variate_analysis_subject.csv")

#### Gene expression ####

adip <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Covariates/Adipose tissue cell composition data/TwinsUK.adipocytes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
macr <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Covariates/Adipose tissue cell composition data/TwinsUK.macrophages.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mvec <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Covariates/Adipose tissue cell composition data/TwinsUK.MVEC.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

covars <- read.table(file = "Adipose_TechCovars_19052017.txt", header = TRUE, quote = "", stringsAsFactors = FALSE)

summary_table <- NULL
  
for (i in 7:8) {
  
  class_i <- classes[!is.na(x = classes[, i]), ]
  
  ids <- Reduce(f = intersect, x = list(class_i$ParticipantID, covars$TwinID, adip$Column, macr$Column, mvec$Column))
  
  class_i <- subset(x = class_i, ParticipantID %in% ids)
  
  covars_i <- subset(x = covars, TwinID %in% ids)
  
  adip_i <- subset(x = adip, Column %in% ids)
  macr_i <- subset(x = macr, Column %in% ids)
  mvec_i <- subset(x = mvec, Column %in% ids)
  
  covars_i <- Reduce(function(x, y, ...) merge(x = x, y = y, by = 1, all = TRUE), x = list(class_i, covars_i, adip_i[, 1:2], macr_i[, 1:2], mvec_i[, 1:2]))
  covars_i[, i] <- as.factor(x = covars_i[, i])
  
  colnames(x = covars_i)[i] <- "Class"
  
  fix <- rbind(    "BMI" = summary(object = glm(formula = Class ~ BMI, family = "binomial", data = covars_i))$coefficients[8],
                   "Age" = summary(object = glm(formula = Class ~ AGE, family = "binomial", data = covars_i))$coefficients[8],
               "GC mean" = summary(object = glm(formula = Class ~ GC_mean, family = "binomial", data = covars_i))$coefficients[8],
      "Duplication rate" = summary(object = glm(formula = Class ~ DupRate, family = "binomial", data = covars_i))$coefficients[8],
         "Fragment mean" = summary(object = glm(formula = Class ~ FragMean, family = "binomial", data = covars_i))$coefficients[8],
    "Insert size median" = summary(object = glm(formula = Class ~ INSERT_SIZE_MEDIAN, family = "binomial", data = covars_i))$coefficients[8],
       "Cos seasonality" = summary(object = glm(formula = Class ~ CosSeasonality, family = "binomial", data = covars_i))$coefficients[8],
       "Sin seasonality" = summary(object = glm(formula = Class ~ SinSeasonality, family = "binomial", data = covars_i))$coefficients[8],
           "Day of year" = summary(object = glm(formula = Class ~ DayofYear, family = "binomial", data = covars_i))$coefficients[8],
            "Adipocytes" = summary(object = glm(formula = Class ~ adipocytes, family = "binomial", data = covars_i))$coefficients[8],
           "Macrophages" = summary(object = glm(formula = Class ~ macrophages, family = "binomial", data = covars_i))$coefficients[8],
                  "MVEC" = summary(object = glm(formula = Class ~ MVEC, family = "binomial", data = covars_i))$coefficients[8])
  fix <- signif(x = fix, digits = 3)
  
  fmX0 <- glm(formula = Class ~ AGE, data = covars_i, family = "binomial")
  
  fmX1 <- glmer(formula = Class ~ AGE + (1|PrimerIndex), data = covars_i, family = "binomial")
  fmX2 <- glmer(formula = Class ~ AGE + (1|date), data = covars_i, family = "binomial")
  fmX3 <- glmer(formula = Class ~ AGE + (1|Zygosity), data = covars_i, family = "binomial")
  fmX4 <- glmer(formula = Class ~ AGE + (1|FatBatch), data = covars_i, family = "binomial")
  
  pbnm_X1 <- pbnm(m1 = fmX1, m0 = fmX0, nsim = 5000, tasks = 10, cores = 2, seed = 3400835)
  pbnm_X2 <- pbnm(m1 = fmX2, m0 = fmX0, nsim = 5000, tasks = 10, cores = 2, seed = 3400835)
  pbnm_X3 <- pbnm(m1 = fmX3, m0 = fmX0, nsim = 5000, tasks = 10, cores = 2, seed = 3400835)
  pbnm_X4 <- pbnm(m1 = fmX4, m0 = fmX0, nsim = 5000, tasks = 10, cores = 2, seed = 3400835)
  
  p1 <- ifelse(test = is.na(x = summary(object = pbnm_X1)$`P(>=obs)`), yes = mean(x = c(summary(object = pbnm_X1)$`#(>=obs)/runs`, summary(object = pbnm_X1)$`#(>=obs)+noEst/runs`)), no = summary(object = pbnm_X1)$`P(>=obs)`)
  p2 <- ifelse(test = is.na(x = summary(object = pbnm_X2)$`P(>=obs)`), yes = mean(x = c(summary(object = pbnm_X2)$`#(>=obs)/runs`, summary(object = pbnm_X2)$`#(>=obs)+noEst/runs`)), no = summary(object = pbnm_X2)$`P(>=obs)`)
  p3 <- ifelse(test = is.na(x = summary(object = pbnm_X3)$`P(>=obs)`), yes = mean(x = c(summary(object = pbnm_X3)$`#(>=obs)/runs`, summary(object = pbnm_X3)$`#(>=obs)+noEst/runs`)), no = summary(object = pbnm_X3)$`P(>=obs)`)
  p4 <- ifelse(test = is.na(x = summary(object = pbnm_X4)$`P(>=obs)`), yes = mean(x = c(summary(object = pbnm_X4)$`#(>=obs)/runs`, summary(object = pbnm_X4)$`#(>=obs)+noEst/runs`)), no = summary(object = pbnm_X4)$`P(>=obs)`)
  
  ran <- rbind("Primer index" = p1, "Date" = p2, "Zygosity" = p3, "Batch" = p4)

  column <- rbind(fix, ran)
  
  colnames(x = column) <- paste0(gsub(pattern = "_", replacement = " ", x = colnames(x = class_i)[i]), " Pr(>|z|)")
  
  summary_table <- cbind(summary_table, column) }
  
write.csv(x = summary_table, file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/3_Co-variate_analysis_expression.csv")

