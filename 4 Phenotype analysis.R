#### Phenotype analysis ####

library(dplyr)
library(lme4)
library(parallel)
library(qvalue)
library(readxl)
library(reshape2)

classes <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Clinical Outcomes/Lipid_regulating_drugs.csv", na = "NA")
master1 <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Output/Master_table.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
subject1 <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Other/Paxgene/TwinsOmicSummary_20190304_ForR.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
subject2 <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Phenobase/MaxTomlinson_22032021.xlsx", sheet = 1)

phenotypes <- Reduce(function(x, y, ...) merge(x, y, by = 1, all = TRUE, sort = FALSE, ...), list(classes[, c("ParticipantID", "Baseline", "End_line")], subject1[, c("Study_Number", "DOB")], subject2[, c("STUDY_NO", "Family_No", "YEAR_BIRTH", "SEX", "ACTUAL_ZYGOSITY")], master1))
phenotypes$VisitDate <- as.Date(x = phenotypes$VisitDate, format = "%Y-%m-%d")
phenotypes$VisitYear <- as.numeric(x = format(x = phenotypes$VisitDate, "%Y"))
phenotypes$AGE <- round(x = as.numeric(x = ((difftime(time1 = phenotypes$VisitDate, time2 = as.Date(x = phenotypes$DOB, format = "%d/%m/%Y"), units = "days")) / 365.2425)), digits = 2)

bmi <- phenotypes[which(x = phenotypes$Description == "Body Mass Index"), ]
bmi$BMI <- bmi$Result

phenotypes <- merge(x = phenotypes, y = bmi[, c("ParticipantID", "VisitDate", "BMI")], by = c("ParticipantID", "VisitDate"), sort = FALSE)
phenotypes <- phenotypes[!is.na(x = phenotypes$Result), -5]
phenotypes <- split(x = phenotypes, f = phenotypes$Description)

master1 <- NULL
master2 <- NULL

for (i in 1:length(x = phenotypes)) {
  
  phenotype <- phenotypes[[i]]
  phenotype[which(x = phenotype$Result %in% boxplot.stats(x = phenotype$Result)$out), "Result"] <- NA
  phenotype <- phenotype[-which(x = is.na(x = phenotype$Result)), ]

  base <- phenotype %>% filter(VisitYear >= 2001 & VisitYear <= 2010) %>% group_by(ParticipantID) %>% arrange(ParticipantID, desc(x = VisitDate)) %>% mutate(VisitOrder = row_number()) %>% filter(VisitOrder == 1)
  
  if (nrow(x = base) >= 30 & var(x = base$Result, na.rm = TRUE) > 0) {
    
    master1 <- rbind(master1, base) }
  
  post <- phenotype %>% filter(VisitYear >= 2001) %>% group_by(ParticipantID) %>% arrange(ParticipantID, VisitDate) %>% mutate(VisitOrder = row_number()) %>% arrange(ParticipantID, desc(x = VisitDate)) %>% mutate(VisitOrder2 = row_number())
  post <- post %>% filter(VisitOrder == 1 | VisitOrder2 == 1) %>% filter(n() > 1) %>% arrange(ParticipantID, VisitDate) %>% mutate(VisitOrder = row_number())
  
  if (nrow(x = post) >= 30 & var(x = post$Result, na.rm = TRUE) > 0) {
    
    post_wide1 <- reshape2::dcast(data = post, formula = ParticipantID ~ VisitOrder, value.var = "VisitDate")
    post_wide2 <- reshape2::dcast(data = post, formula = ParticipantID ~ VisitOrder, value.var = "Result")
    post_wide <- Reduce(function(x, y, ...) merge(x, y, by = 1, all = TRUE, ...), list(post, post_wide1, post_wide2))
    post_wide <- post_wide[which(x = post_wide$VisitOrder == 1), ]
    
    colnames(x = post_wide)[16:19] <- c("T1", "T2", "V1", "V2")
    
    post_wide$timdiff <- post_wide$T2 - post_wide$T1
    post_wide$rawdiff <- post_wide$V2 - post_wide$V1
    post_wide$rate <- post_wide$rawdiff / post_wide$timdiff
    
    change <- lm(formula = rawdiff ~ timdiff + V1, data = post_wide)

    resids <- residuals(object = change)

    post_wide <- cbind(post_wide, as.data.frame(x = resids))
    
    master2 <- rbind(master2, post_wide) } }

baseline <- reshape2::dcast(data = master1[, c(1, 9:10)], formula = ParticipantID ~ Description, value.var = "Result", fun.aggregate = function(x) if(length(x = x) == 0) NA_real_ else sum(x, na.rm = TRUE))

residual <- reshape2::dcast(data = master2[, c(1, 9, 23)], formula = ParticipantID ~ Description, value.var = "resids", fun.aggregate = function(x) if(length(x = x) == 0) NA_real_ else sum(x, na.rm = TRUE))

write.table(x = baseline, file = "~/Documents/KCL/PhD/TwinsUK/Output/Baseline.txt", quote = FALSE, sep = "\t", col.names = make.names(names = colnames(x = baseline), unique = TRUE))
write.table(x = residual, file = "~/Documents/KCL/PhD/TwinsUK/Output/Residual.txt", quote = FALSE, sep = "\t", col.names = make.names(names = colnames(x = residual), unique = TRUE))

master1 <- split(x = master1, f = master1$Description)

results  <- mclapply(X = 1:length(x = master1), FUN = function(i) {
  
  tmp <- data.frame(master1[[i]])
  tmp[, c("Result", "AGE", "BMI", "VisitYear", "YEAR_BIRTH")] <- lapply(X = tmp[, c("Result", "AGE", "BMI", "VisitYear", "YEAR_BIRTH")], FUN = scale)

  pheno = master1[[i]]$Description[1]
  
  formula1 <- paste0("Baseline ~ Result + AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + factor(x = SEX) + (1|Family_No)")
  formula2 <- paste0("Baseline ~ AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + factor(x = SEX) + (1|Family_No)")
  formula3 <- paste0("End_line ~ Result + AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + factor(x = SEX) + (1|Family_No)")
  formula4 <- paste0("End_line ~ AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + factor(x = SEX) + (1|Family_No)")
  
  if (nrow(x = tmp[tmp$SEX == "M", ]) == 0) {
    
    formula1 <- paste0("Baseline ~ Result + AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + (1|Family_No)")
    formula2 <- paste0("Baseline ~ AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + (1|Family_No)")
    formula3 <- paste0("End_line ~ Result + AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + (1|Family_No)")
    formula4 <- paste0("End_line ~ AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + (1|Family_No)") }
  
  results <- tryCatch( {
    
    fm1 = fm2 = fm3 = fm4 = NULL
    fm1 <- glmer(formula = formula1, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm2 <- glmer(formula = formula2, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm3 <- glmer(formula = formula3, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm4 <- glmer(formula = formula4, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    
    p1 <- anova(fm1, fm2)$Pr[2]
    p2 <- anova(fm3, fm4)$Pr[2]
    
    list(pheno = pheno, fm1 = fm1, fm2 = fm2, fm3 = fm3, fm4 = fm4, p1 = p1, p2 = p2, status = "OK", message = as.matrix(x = c("None", "None"), nrow = 2), case = sum(tmp$End_line == 1, na.rm = TRUE), control = sum(tmp$End_line == 0, na.rm = TRUE))
    
  }, warning = function(w) {
    
    fm1 <- glmer(formula = formula1, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm2 <- glmer(formula = formula2, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm3 <- glmer(formula = formula3, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm4 <- glmer(formula = formula4, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    
    p1 <- anova(fm1, fm2)$Pr[2]
    p2 <- anova(fm3, fm4)$Pr[2]
    
    list(pheno = pheno, fm1 = fm1, fm2 = fm2, fm3 = fm3, fm4 = fm4, p1 = p1, p2 = p2, status = "Warning", message = w, case = sum(tmp$End_line == 1, na.rm = TRUE), control = sum(tmp$End_line == 0, na.rm = TRUE))
    
  }, error = function(e) {
    
    list(pheno = pheno, fm1 = NA, fm2 = NA, fm3 = NA, fm4 = NA, p = NA, status = "Error", message = e, case = sum(tmp$End_line == 1, na.rm = TRUE), control = sum(tmp$End_line == 0, na.rm = TRUE)) } )
  
  results
  
  }, mc.cores = 12)

results <- Filter(f = Negate(f = anyNA), x = results)

pheno   <- sapply(X = results, FUN = function(r){ r$pheno })
case    <- sapply(X = results, FUN = function(r){ r$case })
control <- sapply(X = results, FUN = function(r){ r$control })
p1      <- sapply(X = results, FUN = function(r){ r$p1 })
p2      <- sapply(X = results, FUN = function(r){ r$p2 })
status  <- sapply(X = results, FUN = function(r){ r$status })
message <- sapply(X = results, FUN = function(r){ r$message })
beta1   <- sapply(X = results, FUN = function(r){ coef(object = summary(object = r$fm1))[2] })
beta2   <- sapply(X = results, FUN = function(r){ coef(object = summary(object = r$fm3))[2] })
se1     <- sapply(X = results, FUN = function(r){ coef(object = summary(object = r$fm1))[9] })
se2     <- sapply(X = results, FUN = function(r){ coef(object = summary(object = r$fm3))[9] })
resids1 <- sapply(X = results, FUN = function(r){ residuals(object = r$fm3) })

results <- data.frame(cbind(pheno, case, control, beta1, se1, p1, beta2, se2, p2))
results[, 2:9] <- sapply(X = results[, 2:9], FUN = as.numeric)
results[, 2:9] <- signif(x = results[, 2:9], digits = 4)
results$q1 <- qvalue(p = results$p1)$qvalues
results$q2 <- qvalue(p = results$p2)$qvalues

baseline_results <- results

master2 <- split(x = master2, f = master2$Description)

results  <- mclapply(X = 1:length(x = master2), FUN = function(i) {
  
  tmp <- data.frame(master2[[i]])
  tmp[, c("resids", "AGE", "BMI", "VisitYear", "YEAR_BIRTH")] <- lapply(X = tmp[, c("resids", "AGE", "BMI", "VisitYear", "YEAR_BIRTH")], FUN = scale)
  
  pheno = master2[[i]]$Description[1]
  
  formula1 <- paste0("Baseline ~ resids + AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + factor(x = SEX) + (1|Family_No)")
  formula2 <- paste0("Baseline ~ AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + factor(x = SEX) + (1|Family_No)")
  formula3 <- paste0("End_line ~ resids + AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + factor(x = SEX) + (1|Family_No)")
  formula4 <- paste0("End_line ~ AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + factor(x = SEX) + (1|Family_No)")
  
  if (nrow(x = tmp[tmp$SEX == "M", ]) == 0) {
    
    formula1 <- paste0("Baseline ~ resids + AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + (1|Family_No)")
    formula2 <- paste0("Baseline ~ AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + (1|Family_No)")
    formula3 <- paste0("End_line ~ resids + AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + (1|Family_No)")
    formula4 <- paste0("End_line ~ AGE + BMI + VisitYear + YEAR_BIRTH + factor(x = ACTUAL_ZYGOSITY) + (1|Family_No)") }
  
  results <- tryCatch( {
    
    fm1 = fm2 = fm3 = fm4 = NULL
    fm1 <- glmer(formula = formula1, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm2 <- glmer(formula = formula2, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm3 <- glmer(formula = formula3, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm4 <- glmer(formula = formula4, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    
    p1 <- anova(fm1, fm2)$Pr[2]
    p2 <- anova(fm3, fm4)$Pr[2]
    
    list(pheno = pheno, fm1 = fm1, fm2 = fm2, fm3 = fm3, fm4 = fm4, p1 = p1, p2 = p2, status = "OK", message = as.matrix(x = c("None", "None"), nrow = 2), case = sum(tmp$End_line == 1, na.rm = TRUE), control = sum(tmp$End_line == 0, na.rm = TRUE))
    
  }, warning = function(w) {
    
    fm1 <- glmer(formula = formula1, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm2 <- glmer(formula = formula2, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm3 <- glmer(formula = formula3, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    fm4 <- glmer(formula = formula4, data = tmp, family = binomial, control = glmerControl(optimizer = "optimx", optCtrl = list(method = 'nlminb')))
    
    p1 <- anova(fm1, fm2)$Pr[2]
    p2 <- anova(fm3, fm4)$Pr[2]
    
    list(pheno = pheno, fm1 = fm1, fm2 = fm2, fm3 = fm3, fm4 = fm4, p1 = p1, p2 = p2, status = "Warning", message = w, case = sum(tmp$End_line == 1, na.rm = TRUE), control = sum(tmp$End_line == 0, na.rm = TRUE))
    
  }, error = function(e) {
    
    list(pheno = pheno, fm1 = NA, fm2 = NA, fm3 = NA, fm4 = NA, p = NA, status = "Error", message = e, case = sum(tmp$End_line == 1, na.rm = TRUE), control = sum(tmp$End_line == 0, na.rm = TRUE)) } )
  
  results
  
  }, mc.cores = 12)

results <- Filter(f = Negate(f = anyNA), x = results)

pheno   <- sapply(X = results, FUN = function(r){ r$pheno })
case    <- sapply(X = results, FUN = function(r){ r$case })
control <- sapply(X = results, FUN = function(r){ r$control })
p1      <- sapply(X = results, FUN = function(r){ r$p1 })
p2      <- sapply(X = results, FUN = function(r){ r$p2 })
status  <- sapply(X = results, FUN = function(r){ r$status })
message <- sapply(X = results, FUN = function(r){ r$message })
beta1   <- sapply(X = results, FUN = function(r){ coef(object = summary(object = r$fm1))[2] })
beta2   <- sapply(X = results, FUN = function(r){ coef(object = summary(object = r$fm3))[2] })
se1     <- sapply(X = results, FUN = function(r){ coef(object = summary(object = r$fm1))[9] })
se2     <- sapply(X = results, FUN = function(r){ coef(object = summary(object = r$fm3))[9] })
resids2 <- sapply(X = results, FUN = function(r){ residuals(object = r$fm3) })

results <- data.frame(cbind(pheno, case, control, beta1, se1, p1, beta2, se2, p2))
results[, 2:9] <- sapply(X = results[, 2:9], FUN = as.numeric)
results[, 2:9] <- signif(x = results[, 2:9], digits = 4)
results$q1 <- qvalue(p = results$p1)$qvalues
results$q2 <- qvalue(p = results$p2)$qvalues

residual_results <- results

baseline_results <- baseline_results[baseline_results$q2 < 0.05, ]

residual_results <- residual_results[residual_results$q2 < 0.05, ]

write.table(x = as.matrix(x = baseline_results), file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/4_Phenotype_analysis_baseline.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(x = as.matrix(x = residual_results), file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/4_Phenotype_analysis_residual.txt", row.names = FALSE, sep = "\t", quote = FALSE)

baseline <- baseline[, c(1, which(x = match(x = colnames(x = baseline), table = baseline_results$pheno) > 0))]

residual <- residual[, c(1, which(x = match(x = colnames(x = residual), table = residual_results$pheno) > 0))]

write.table(x = baseline, file = "~/Documents/KCL/PhD/TwinsUK/Output/Baseline_q5.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = make.names(names = colnames(x = baseline), unique = TRUE))
write.table(x = residual, file = "~/Documents/KCL/PhD/TwinsUK/Output/Residual_q5.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = make.names(names = colnames(x = residual), unique = TRUE))
