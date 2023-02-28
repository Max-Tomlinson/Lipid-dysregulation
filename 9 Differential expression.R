#### Differential expression ####

library(DESeq2)
library(IHW)

diab <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Other/Diabetes/AllTUK_DefinitiveDiabetesStatus_UpdatedJuly2018[1].txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

classes <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Clinical Outcomes/Lipid_regulating_drugs_censored.csv")

gene.count <- read.delim(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Raw Counts/swap.merged.fat.gene.count.txt", check.names = FALSE, stringsAsFactors = FALSE)
gene.count <- as.data.frame(x = gene.count[, -c(1:6, which(x = colnames(x = gene.count) == "26882"), which(x = colnames(x = gene.count) == "51512"), which(x = colnames(x = gene.count) == "99282"))], row.names = gene.count[, 4])

covariates <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Covariates/Adipose_TechCovars_19052017.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

geneinfo <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Gencode_V19_GeneInfo_23052017.txt", header = TRUE, stringsAsFactors = FALSE)
geneinfo <- merge(x = geneinfo[, 1:2], y = gene.count[, 0], by.x = "ensemblID", by.y = 0)

summary_table1 <- NULL

for (i in 2:ncol(x = classes)) {
  
  class_i <- classes[!is.na(x = classes[, i]), c(1, i)]
  
  ids <- Reduce(f = intersect, x = list(class_i$ParticipantID, colnames(x = gene.count), covariates$TwinID, diab$TwinID))
  
  class_i <- subset(x = class_i, ParticipantID %in% ids)
  
  countData <- t(x = subset(x = t(x = gene.count), colnames(x = gene.count) %in% ids))
  countData <- countData[order(rownames(x = countData)), order(as.integer(x = colnames(x = countData)))]
  
  covars_i <- subset(x = covariates, TwinID %in% ids)
  
  colData <- merge(x = covars_i[, c("TwinID", "AGE", "BMI")], y = diab[, c("TwinID", "Diab_01")], by = 1)
  colData <- merge(x = colData, y = class_i, by = 1, all = TRUE)
  colData[is.na(x = colData$Diab_01), ]$Diab_01 <- 0
  colData[, c("AGE", "BMI")] <- scale(x = colData[, c("AGE", "BMI")])
  colData[4:5] <- lapply(X = colData[4:5], FUN = factor)
  
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = as.formula(object = paste0(" ~ AGE + BMI + Diab_01 + ", colnames(x = colData)[5])))
  dds <- dds[rowSums(x = counts(object = dds) >= 10) >= (min(table(class_i[, 2]))), ]
  dds <- DESeq(object = dds)
  
  res <- results(object = dds, contrast = c(colnames(x = colData)[5], "1", "0"), filterFun = ihw)
  res <- merge(x = geneinfo[, 1:2], y = as.data.frame(x = res), by.x = "ensemblID", by.y = 0, sort = FALSE)
  
  write.csv(x = res[order(res$padj), ], file = paste0("~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/9_", gsub(pattern = "\\.", replacement = "_", x = colnames(x = class_i)[2], ".csv"), "_deseq.csv"), row.names = FALSE)
  
  summary_deseq <- cbind(length(x = ids), table(class_i[, 2])[1], table(class_i[, 2])[2], table(res$padj < 0.1)[2])
  
  colnames(x = summary_deseq) <- c("N", "Controls", "Cases", "FDR 10%")
  rownames(x = summary_deseq) <- colnames(x = classes)[i]

  summary_table1 <- rbind(summary_table1, summary_deseq) }

write.csv(x = summary_table1, file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/9_Summary_deseq.csv")