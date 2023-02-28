#### Robust PCA ####

library(DESeq2)
library(ggplot2)
library(ggsci)
library(rrcov)
library(scales)

classes <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Clinical Outcomes/Lipid_regulating_drugs.csv")

gene.count <- read.delim(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Raw Counts/swap.merged.fat.gene.count.txt", check.names = FALSE, stringsAsFactors = FALSE)
gene.count <- as.data.frame(x = gene.count[, -c(1:6)], row.names = gene.count[, 4])

covariates <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/EuroBATS/Covariates/Adipose_TechCovars_19052017.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ids <- Reduce(f = intersect, x = list(colnames(x = gene.count), covariates$TwinID))

countData <- t(x = subset(x = t(x = gene.count), colnames(x = gene.count) %in% ids))

colData <- subset(x = covariates, TwinID %in% ids)
colData[, c("AGE", "BMI")] <- scale(x = colData[, c("AGE", "BMI")])

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = as.formula(object = paste0(" ~ BMI + AGE")))

keep <- rowSums(x = counts(object = dds) >= 10) >= 92

dds <- dds[keep,]
dds <- estimateSizeFactors(object = dds)

vsd <- vst(object = dds, blind = FALSE)

pca <- PcaGrid(x = t(x = assay(x = vsd)))

data <- cbind.data.frame(`Score distance` = pca$sd, Index = names(x = pca$sd))
data <- merge(x = data, y = classes[, c(1, 8)], by.x = "Index", by.y = "ParticipantID")
data$Group <- ifelse(test = data$End_line == 1, yes = "Drug", no = "Control")

p1 <- ggplot(data = data, aes(x = factor(x = Index), y = `Score distance`)) + 
  geom_point() + geom_hline(yintercept = pca$cutoff.sd, colour = "red") + geom_text(mapping = aes(label = ifelse(test = `Score distance` > pca$cutoff.sd, yes = as.character(x = Index), no = '')), hjust = 1.25, vjust = 0) + 
  scale_colour_jama() + scale_fill_jama() + scale_y_continuous(breaks = pretty_breaks()) + 
  labs(title = "Robust PCA plot of variance stabilising transformed counts", x = "Sample index", y = "Score distance") + 
  theme_classic(base_size = 18) + theme(axis.text.x = element_blank(), axis.ticks = element_blank())

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/8_Robust_PCA_vst.pdf", plot = p1, width = 11, height = 8.5)

data <- data[-which(x = is.na(x = data$End_line)), ]

p1 <- ggplot(data = data, aes(x = factor(x = Index), y = `Score distance`, colour = Group)) + 
  geom_point() + geom_hline(yintercept = pca$cutoff.sd, colour = "red") + geom_text(mapping = aes(label = ifelse(test = `Score distance` > pca$cutoff.sd, yes = as.character(x = Index), no = '')), hjust = 1.25, vjust = 0) + 
  scale_colour_jama() + scale_fill_jama() + scale_y_continuous(breaks = pretty_breaks()) + 
  labs(title = "Robust PCA plot of variance stabilising transformed counts", x = "Sample index", y = "Score distance") + 
  theme_classic(base_size = 18) + theme(axis.text.x = element_blank(), axis.ticks = element_blank())

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/8_Robust_PCA_vst_group.pdf", plot = p1, width = 11, height = 8.5)
