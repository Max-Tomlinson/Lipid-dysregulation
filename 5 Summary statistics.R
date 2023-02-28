#### Summary statistics ####

library(dplyr)
library(parallel)
library(readxl)

classes <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Clinical Outcomes/Lipid_regulating_drugs.csv", na = "NA")

subject1 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_22032021.xlsx", sheet = 1)

subject2 <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Other/Paxgene/TwinsOmicSummary_20190304_ForR.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

master1 <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Output/Master_table.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

phenotypes <- Reduce(function(x, y, ...) merge(x, y, by = 1, all = TRUE, sort = FALSE, ...), list(classes[, c("ParticipantID", "Baseline", "End_line")], subject1[, c("STUDY_NO", "Family_No", "YEAR_BIRTH", "SEX", "ACTUAL_ZYGOSITY")], subject2[, c("Study_Number", "DOB")], master1))
phenotypes$VisitDate <- as.Date(x = phenotypes$VisitDate, format = "%Y-%m-%d")
phenotypes$VisitYear <- as.numeric(x = format(x = phenotypes$VisitDate, "%Y"))
phenotypes$AGE <- round(x = as.numeric(x = ((difftime(time1 = phenotypes$VisitDate, time2 = as.Date(x = phenotypes$DOB, format = "%d/%m/%Y"), units = "days")) / 365.2425)), digits = 2)

bmi <- phenotypes[which(x = phenotypes$Description == "Body Mass Index"), ]
bmi$BMI <- bmi$Result

phenotypes <- merge(x = phenotypes, y = bmi[, c("ParticipantID", "VisitDate", "BMI")], by = c("ParticipantID", "VisitDate"), sort = FALSE)
phenotypes <- phenotypes[!is.na(x = phenotypes$Result), -9]
phenotypes <- split(x = phenotypes, f = phenotypes$Description)

master1 <- NULL
master2 <- NULL

for (i in 1:length(x = phenotypes)) {
  
  phenotype <- phenotypes[[i]]
  phenotype[which(x = phenotype$Result %in% boxplot.stats(x = phenotype$Result)$out), "Result"] <- NA
  phenotype <- phenotype[-which(x = is.na(x = phenotype$Result)), ]
  
  base <- phenotype %>% filter(VisitYear >= 2001 & VisitYear <= 2010) %>% group_by(ParticipantID) %>% arrange(ParticipantID, desc(x = VisitDate)) %>% mutate(VisitOrder = row_number()) %>% filter(VisitOrder == 1)
  
  if (nrow(x = base) >= 30 & var(x = base$Result, na.rm = TRUE) > 0) {
    
    master1 <- rbind(master1, base) } }

master0 <- master1[!duplicated(x = master1[, c(1:8, 11:14)]), ]
master1 <- split(x = master1, f = master1$Description)

results  <- mclapply(X = 1:length(x = master1), FUN = function(i) {
  
  tmp <- data.frame(master1[[i]])
  tmp <- tmp[!is.na(x = tmp$End_line), ]
  
  pheno = master1[[i]]$Description[1]
  
  results <- tryCatch( {
    
    case = sum(tmp$End_line == 1, na.rm = TRUE)
    cont = sum(tmp$End_line == 0, na.rm = TRUE)

    male = sum(tmp$SEX == "M", na.rm = TRUE)
    fmle = sum(tmp$SEX == "F", na.rm = TRUE)
    
    mz = sum(tmp$ACTUAL_ZYGOSITY == "MZ", na.rm = TRUE)
    dz = sum(tmp$ACTUAL_ZYGOSITY == "DZ", na.rm = TRUE)
    
    visit_case = as.numeric(x = tapply(X = tmp$VisitYear, INDEX = tmp$End_line, FUN = mean)[2])
    visit_cont = as.numeric(x = tapply(X = tmp$VisitYear, INDEX = tmp$End_line, FUN = mean)[1])
    visit_p <- function(...) {
      obj <- try(expr = t.test(x = tmp[tmp$End_line == 0, ]$VisitYear, y = tmp[tmp$End_line == 1, ]$VisitYear)[3], silent = TRUE)
      if (is(object = obj, class2 = "try-error")) return("Error") else return(obj$p.value)
    }
    visit_p <- visit_p()
    
    year_case = as.numeric(x = tapply(X = tmp$YEAR_BIRTH, INDEX = tmp$End_line, FUN = mean)[2])
    year_cont = as.numeric(x = tapply(X = tmp$YEAR_BIRTH, INDEX = tmp$End_line, FUN = mean)[1])
    year_p    = as.numeric(x = t.test(x = tmp[tmp$End_line == 0, ]$YEAR_BIRTH, y = tmp[tmp$End_line == 1, ]$YEAR_BIRTH)[3])

    res_case = as.numeric(x = tapply(X = tmp$Result, INDEX = tmp$End_line, FUN = mean)[2])
    res_cont = as.numeric(x = tapply(X = tmp$Result, INDEX = tmp$End_line, FUN = mean)[1])
    res_p    = as.numeric(x = t.test(x = tmp[tmp$End_line == 0, ]$Result, y = tmp[tmp$End_line == 1, ]$Result)[3])
    
    age_case = as.numeric(x = tapply(X = tmp$AGE, INDEX = tmp$End_line, FUN = mean)[2])
    age_cont = as.numeric(x = tapply(X = tmp$AGE, INDEX = tmp$End_line, FUN = mean)[1])
    age_p    = as.numeric(x = t.test(x = tmp[tmp$End_line == 0, ]$AGE, y = tmp[tmp$End_line == 1, ]$AGE)[3])
    
    bmi_case = as.numeric(x = tapply(X = tmp$BMI, INDEX = tmp$End_line, FUN = mean)[2])
    bmi_cont = as.numeric(x = tapply(X = tmp$BMI, INDEX = tmp$End_line, FUN = mean)[1])
    bmi_p    = as.numeric(x = t.test(x = tmp[tmp$End_line == 0, ]$BMI, y = tmp[tmp$End_line == 1, ]$BMI)[3])
    
    list(pheno, case, cont, male, fmle, mz, dz,
         visit_case, visit_cont, visit_p,
         year_case, year_cont, year_p,
         res_case, res_cont, res_p,
         age_case, age_cont, age_p,
         bmi_case, bmi_cont, bmi_p) } )
  
  results
  }, mc.cores = 12)

results <- Filter(f = Negate(f = anyNA), x = results)

pheno       <- sapply(X = results, FUN = function(r){ r[1] })
case        <- sapply(X = results, FUN = function(r){ r[2] })
cont        <- sapply(X = results, FUN = function(r){ r[3] })
male        <- sapply(X = results, FUN = function(r){ r[4] })
fmle        <- sapply(X = results, FUN = function(r){ r[5] })
mz          <- sapply(X = results, FUN = function(r){ r[6] })
dz          <- sapply(X = results, FUN = function(r){ r[7] })
visit_case  <- sapply(X = results, FUN = function(r){ r[8] })
visit_cont  <- sapply(X = results, FUN = function(r){ r[9] })
visit_p     <- sapply(X = results, FUN = function(r){ r[10] })
year_case   <- sapply(X = results, FUN = function(r){ r[11] })
year_cont   <- sapply(X = results, FUN = function(r){ r[12] })
year_p      <- sapply(X = results, FUN = function(r){ r[13] })
res_case    <- sapply(X = results, FUN = function(r){ r[14] })
res_cont    <- sapply(X = results, FUN = function(r){ r[15] })
res_p       <- sapply(X = results, FUN = function(r){ r[16] })
age_case    <- sapply(X = results, FUN = function(r){ r[17] })
age_cont    <- sapply(X = results, FUN = function(r){ r[18] })
age_p       <- sapply(X = results, FUN = function(r){ r[19] })
bmi_case    <- sapply(X = results, FUN = function(r){ r[20] })
bmi_cont    <- sapply(X = results, FUN = function(r){ r[21] })
bmi_p       <- sapply(X = results, FUN = function(r){ r[22] })

results <- data.frame(cbind(pheno, case, cont, male, fmle, mz, dz,
                            visit_case, visit_cont, visit_p,
                            year_case, year_cont, year_p,
                            res_case, res_cont, res_p,
                            age_case, age_cont, age_p,
                            bmi_case, bmi_cont, bmi_p))

master0 <- master0[!duplicated(x = master0[, 1]), ]
master0 <- master0[!is.na(x = master0$End_line), ]

age <- c("Age", sum(master0$End_line == 1, na.rm = TRUE), sum(master0$End_line == 0, na.rm = TRUE), sum(master0$SEX == "M", na.rm = TRUE), sum(master0$SEX == "F", na.rm = TRUE), sum(master0$ACTUAL_ZYGOSITY == "MZ", na.rm = TRUE), sum(master0$ACTUAL_ZYGOSITY == "DZ", na.rm = TRUE),
         as.numeric(x = tapply(X = master0$VisitYear, INDEX = master0$End_line, FUN = mean)[2]), as.numeric(x = tapply(X = master0$VisitYear, INDEX = master0$End_line, FUN = mean)[1]), as.numeric(x = t.test(x = master0[master0$End_line == 1, ]$VisitYear, y = master0[master0$End_line == 0, ]$VisitYear)[3]),
         as.numeric(x = tapply(X = master0$YEAR_BIRTH, INDEX = master0$End_line, FUN = mean)[2]), as.numeric(x = tapply(X = master0$YEAR_BIRTH, INDEX = master0$End_line, FUN = mean)[1]), as.numeric(x = t.test(x = master0[master0$End_line == 1, ]$YEAR_BIRTH, y = master0[master0$End_line == 0, ]$YEAR_BIRTH)[3]),
         as.numeric(x = tapply(X = master0$AGE, INDEX = master0$End_line, FUN = mean)[2]), as.numeric(x = tapply(X = master0$AGE, INDEX = master0$End_line, FUN = mean)[1]), as.numeric(x = t.test(x = master0[master0$End_line == 1, ]$AGE, y = master0[master0$End_line == 0, ]$AGE)[3]), "NA", "NA", "NA", "NA", "NA", "NA")

bmi <- c("BMI", sum(master0$End_line == 1, na.rm = TRUE), sum(master0$End_line == 0, na.rm = TRUE), sum(master0$SEX == "M", na.rm = TRUE), sum(master0$SEX == "F", na.rm = TRUE), sum(master0$ACTUAL_ZYGOSITY == "MZ", na.rm = TRUE), sum(master0$ACTUAL_ZYGOSITY == "DZ", na.rm = TRUE),
         as.numeric(x = tapply(X = master0$VisitYear, INDEX = master0$End_line, FUN = mean)[2]), as.numeric(x = tapply(X = master0$VisitYear, INDEX = master0$End_line, FUN = mean)[1]), as.numeric(x = t.test(x = master0[master0$End_line == 1, ]$VisitYear, y = master0[master0$End_line == 0, ]$VisitYear)[3]),
         as.numeric(x = tapply(X = master0$YEAR_BIRTH, INDEX = master0$End_line, FUN = mean)[2]), as.numeric(x = tapply(X = master0$YEAR_BIRTH, INDEX = master0$End_line, FUN = mean)[1]), as.numeric(x = t.test(x = master0[master0$End_line == 1, ]$YEAR_BIRTH, y = master0[master0$End_line == 0, ]$YEAR_BIRTH)[3]),
         as.numeric(x = tapply(X = master0$BMI, INDEX = master0$End_line, FUN = mean)[2]), as.numeric(x = tapply(X = master0$BMI, INDEX = master0$End_line, FUN = mean)[1]), as.numeric(x = t.test(x = master0[master0$End_line == 1, ]$BMI, y = master0[master0$End_line == 0, ]$BMI)[3]), "NA", "NA", "NA", "NA", "NA", "NA")

results <- rbind(age, bmi, results)
results[, 2:22] <- sapply(X = results[, 2:22], FUN = as.numeric)
results[, 2:22] <- signif(x = results[, 2:22], digits = 4)

colnames(x = results) <- c("Phenotype", "Cases", "Controls", "Male", "Female", "MZs", "DZs",
                           "Mean year of visit case", "Mean year of visit control", "Year of visit p-value",
                           "Mean year of birth case", "Mean year of birth control", "Year of birth p-value",
                           "Mean result case", "Mean result control", "Result p-value", 
                           "Mean age case", "Mean age control", "Age p-value", 
                           "Mean BMI case", "Mean BMI control", "BMI p-value")

write.table(x = as.matrix(x = results), file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/5_Summary_phenotypes.txt", quote = FALSE, sep = "\t", row.names = FALSE)

baseline <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Output/Baseline_q5.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

results <- results[c(1:2, which(x = match(x = gsub(pattern = '[[:punct:] ]', replacement = '.', x = results$Phenotype), table = colnames(x = baseline)) > 0)), ]

write.table(x = as.matrix(x = results[, c(1:3, 14:16)]), file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/5_Summary_phenotypes_short.txt", quote = FALSE, sep = "\t", row.names = FALSE)
