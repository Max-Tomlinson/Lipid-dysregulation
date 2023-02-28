#### Molecular markers ####

library(data.table)
library(ggplot2)
library(ggsci)
library(gprofiler2)
library(rstatix)
library(scales)
library(qvalue)

gencode <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EUROBATS/Gencode_V19_GeneInfo_23052017.txt", header = TRUE)

tests <- 0

diff.exp <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/9_End_line_deseq.csv")
diff.exp$log2FoldChange.abs <- abs(x = diff.exp$log2FoldChange)

df1 <- NULL

for (k in list.files(path = "~/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose", pattern = ".txt", full.names = TRUE, recursive = TRUE)) { 
  
  lmem <- read.table(file = k, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  lmem <- merge(x = gencode[, 1:2], y = lmem, by.x = "ensemblID", by.y = "gene")
  lmem$q <- qvalue(p = lmem$p)$qvalue
  lmem$significant <- lmem$q < 0.05
  
  diff.lmem <- as.data.frame(x = matrix(data = NA, nrow = 1, ncol = 15))
  
  colnames(x = diff.lmem) <- c("estimate", "estimate1", "estimate2", ".y.", "group1", "group2", "n1", "n2", "statistic", "p", "df", "conf.low", "conf.high", "method", "alternative")
  
  if (nrow(x = lmem[lmem$significant == TRUE, ]) > 1) {
    
    diff.lmem <- merge(x = diff.exp, y = lmem, by = c("GeneName", "ensemblID"))
    diff.lmem <- t_test(data = diff.lmem, formula = log2FoldChange.abs ~ significant, p.adjust.method = "bonferroni", alternative = "less", detailed = TRUE)
    
    tests = tests + 1 }
  
  df1 <- rbind(df1, diff.lmem) }

df1 <- as.data.frame(x = df1)

rownames(x = df1) <- list.files(path = "~/Documents/KCL/PhD/TwinsUK/Output/GEAS/Adipose", pattern = ".txt", full.names = TRUE, recursive = TRUE)

df1 <- na.omit(object = df1)
df1$significant <- df1$p < 0.05/tests

df2 <- df1[df1$significant == TRUE, ]

names <- gsub(pattern = ".*/", replacement = "", x = rownames(x = df1))
names <- gsub(pattern = "\\..*", replacement = "", x = names)
names <- gsub(pattern = "_", replacement = " ", x = names)
names <- gsub(pattern = "\\  .*", replacement = "", x = names)
names <- gsub(pattern = "20230130", replacement = "", x = names)
names <- gsub(pattern = "20230131", replacement = "", x = names)
names <- trimws(x = names)

rownames(x = df1) <- make.unique(names = names)

write.csv(x = df1, file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/10_DE_results.csv", quote = FALSE, row.names = TRUE)

terms <- NULL

ids   <- NULL

for (i in rownames(x = df2)) {
  
  results <- read.table(file = i, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  results <- merge(x = gencode[, 1:2], y = results, by.x = "ensemblID", by.y = "gene")
  
  sig <- results[qvalue(p = results$p)$qvalues < 0.05, ]
  sig <- unique(x = sig[order(sig$p), ]$GeneName)
  
  ids <- c(ids, sig)
  
  if (length(x = sig) > 2) {
    
    gostres <- gost(query = sig, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = results$GeneName, sources = c("CORUM", "HP", "GO", "KEGG", "REAC")) }
  
  terms <- c(terms, unique(x = gostres$result$term_name)) }

gostres <- gost(query = unique(x = ids), organism = "hsapiens", ordered_query = FALSE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = results$GeneName, sources = "KEGG")

publish_gosttable(gostres = gostres, filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/gosttables/10_gosttable_all.pdf")

df <- as.data.frame(x = table(ids))

for (i in max(df$Freq):1) {
  
  query <- unique(x = df[df$Freq >= i, ])
  
  if (nrow(x = query) > 10 & nrow(x = query) < 542) {
    
    gostres <- gost(query = query$ids, organism = "hsapiens", ordered_query = FALSE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = results$GeneName)
    
    publish_gosttable(gostres = gostres, filename = paste0("~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/gosttables/10_gosttable_", i, ".pdf")) } 
  
  if (nrow(x = query) >= 542) {
    
    gostres <- gost(query = query$ids, organism = "hsapiens", ordered_query = FALSE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = results$GeneName, sources = "KEGG")
    
    publish_gosttable(gostres = gostres, filename = paste0("~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/gosttables/10_gosttable_", i, ".pdf")) } }

term_table <- as.data.frame(x = table(terms))

f10 <- ggplot(data = term_table[term_table$Freq >= 7, ], aes(x = reorder(x = terms, -Freq), y = Freq, fill = reorder(x = terms, -Freq))) + geom_bar(stat = "identity", position = position_dodge()) + 
  labs(title = " ", x = " ", y = " ") + guides(fill = "none") +
  theme_classic(base_size = 18) + theme(axis.text.x = element_text(angle = 270, hjust = 0), axis.title.x = element_blank()) + scale_y_sqrt()

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/10_Term_barplot.pdf", plot = f10, width = 11, height = 8.5)
