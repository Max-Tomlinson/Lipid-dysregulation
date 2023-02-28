#### Gene expression association study ####

library(gprofiler2)
library(qvalue)

gencode <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EUROBATS/Gencode_V19_GeneInfo_23052017.txt", header = TRUE)

res <- NULL

for (i in list.files(path = "~/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/Baseline", pattern = ".txt", full.names = TRUE, recursive = TRUE)) {
  
  lmem <- read.table(file = i, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  lmem <- merge(x = gencode[, 1:2], y = lmem, by.x = "ensemblID", by.y = "gene")

  val <- nrow(x = lmem[qvalue(p = lmem$p)$qvalues < 0.05, ])
  
  res <- rbind(res, cbind(lmem$n[1], val))
  
  sig <- lmem[qvalue(p = lmem$p)$qvalues < 0.05, ]
  sig <- sig[order(sig$p), ]
  
  if (nrow(x = sig) > 4 &
      i != "/Users/MaxTomlinson/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/Baseline/All_tissue_in_gynoid_region_20230208.txt" &
      i != "/Users/MaxTomlinson/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/Baseline/Cardiac_ultrasound___Intima_media_thickness_of_anterior_wall_of_left_carotid_artery__20230208.txt" &
      i != "/Users/MaxTomlinson/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/Baseline/ECG_RR_interval_based_on_automatic_measurements__The_time_elapsing_between_two_consecutive_R_waves_in_the_electrocardiogram_20230208.txt") {
    
    gostres <- gost(query = sig$GeneName, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = lmem$GeneName, sources = c("KEGG"))
    
    publish_gosttable(gostres = gostres, filename = paste0("~/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/", paste0(gsub(pattern = "\\..*", replacement = "", x = gsub(pattern = ".*/", replacement = "", x = i)), ".pdf"))) }
  
  if (
    i == "/Users/MaxTomlinson/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/Baseline/All_tissue_in_gynoid_region_20230208.txt" |
    i == "/Users/MaxTomlinson/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/Baseline/Cardiac_ultrasound___Intima_media_thickness_of_anterior_wall_of_left_carotid_artery__20230208.txt" |
    i == "/Users/MaxTomlinson/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/Baseline/ECG_RR_interval_based_on_automatic_measurements__The_time_elapsing_between_two_consecutive_R_waves_in_the_electrocardiogram_20230208.txt") {
    
    gostres <- gost(query = sig$GeneName, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = lmem$GeneName)
    
    publish_gosttable(gostres = gostres, filename = paste0("~/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/", paste0(gsub(pattern = "\\..*", replacement = "", x = gsub(pattern = ".*/", replacement = "", x = i)), ".pdf"))) } }

res <- as.data.frame(x = res)

rownames(x = res) <- list.files(path = "~/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/Baseline", pattern = ".txt", full.names = TRUE, recursive = TRUE)

names <- gsub(pattern = ".*/", replacement = "", x = rownames(x = res))
names <- gsub(pattern = "\\..*", replacement = "", x = names)
names <- gsub(pattern = "_", replacement = " ", x = names)
names <- gsub(pattern = "\\  .*", replacement = "", x = names)
names <- gsub(pattern = "20230208", replacement = "", x = names)
names <- gsub(pattern = "20230209", replacement = "", x = names)
names <- trimws(x = names)

rownames(x = res) <- make.unique(names = names)
colnames(x = res) <- c("Sample", "FDR 5%")

res$`FDR 5%` <- ifelse(test = res$`FDR 5%` == 0, yes = NA, no = res$`FDR 5%`)
res1 <- na.omit(object = res)

res <- NULL

for (i in list.files(path = "~/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/Residual/", pattern = ".txt", full.names = TRUE, recursive = TRUE)) {
  
  lmem <- read.table(file = i, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  lmem <- merge(x = gencode[, 1:2], y = lmem, by.x = "ensemblID", by.y = "gene")
  
  val <- nrow(x = lmem[qvalue(p = lmem$p)$qvalues < 0.05, ])
  
  res <- rbind(res, cbind(lmem$n[1], val)) 
  
  sig <- lmem[qvalue(p = lmem$p)$qvalues < 0.05, ]
  sig <- sig[order(sig$p), ]
  
  if (nrow(x = sig) > 1) {

    gostres <- gost(query = sig$GeneName, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = lmem$GeneName, sources = c("KEGG"))
    
    publish_gosttable(gostres = gostres, filename = paste0("~/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/", paste0(gsub(pattern = "\\..*", replacement = "", x = gsub(pattern = ".*/", replacement = "", x = i)), ".pdf"))) } }

res <- as.data.frame(x = res)

rownames(x = res) <- list.files(path = "~/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose/Residual", pattern = ".txt", full.names = TRUE, recursive = TRUE)

names <- gsub(pattern = ".*/", replacement = "", x = rownames(x = res))
names <- gsub(pattern = "\\..*", replacement = "", x = names)
names <- gsub(pattern = "_", replacement = " ", x = names)
names <- gsub(pattern = "20230208", replacement = "", x = names)
names <- gsub(pattern = "20230209", replacement = "", x = names)
names <- trimws(x = names)

rownames(x = res) <- make.unique(names = names)
colnames(x = res) <- c("Sample", "FDR 5%")

res$`FDR 5%` <- ifelse(test = res$`FDR 5%` == 0, yes = NA, no = res$`FDR 5%`)

res2 <- na.omit(object = res)

write.csv(x = res1, file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/6_GEAS_Baseline.csv", quote = FALSE, row.names = TRUE)
write.csv(x = res2, file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/6_GEAS_Residual.csv", quote = FALSE, row.names = TRUE)

write.csv(x = rbind(res1, res2), file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/6_GEAS.csv", quote = FALSE, row.names = TRUE)
