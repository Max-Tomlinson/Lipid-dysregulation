#### The incidence of lipid regulation ####

library(ggplot2)
library(ggsci)
library(lubridate)
library(readxl)
library(scales)
library(tidyr)

medication <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Medication/Prescriptions/TwinsUK_Prescribed_Medication_V1.0_2021.xlsx", sheet = 1, na = "NA")
medication <- medication[!duplicated(x = medication[, c(1, 8)], fromLast = FALSE), ]
medication$Lipid_regulating_drugs <- ifelse(test = grepl(pattern = "Lipid-Regulating Drugs", x = medication$Class, ignore.case = TRUE), yes = 1, no = 0)
medication$Lipid_regulating_drugs <- ifelse(test = medication$Source == "Q11B" & medication$Date > "2010-12-20", yes = NA, no = medication$Lipid_regulating_drugs)
medication <- na.omit(object = medication)

lrd <- spread(data = medication[, c(1, 8:9)], key = "Source", value = "Lipid_regulating_drugs")
lrd$Baseline <- ifelse(test = rowMeans(x = lrd[, 3:4], na.rm = TRUE) > 0, yes = 1, no = 0)
lrd$End_line <- ifelse(test = rowMeans(x = lrd[, c(2, 5:6)], na.rm = TRUE) > 0, yes = 1, no = 0)
lrd$Final_call <- ifelse(test = rowMeans(x = lrd[, -1], na.rm = TRUE) > 0, yes = 1, no = 0)

for (i in which(x = lrd$Baseline == 1)) { lrd$End_line[i] <- NA }
for (i in which(x = lrd$End_line == 0)) { lrd$Baseline[i] <- 0 }

dates <- as.data.frame(x = spread(data = medication[, c(1:2, 8)], key = "Source", value = "Date"))

range <- NULL

for (i in 2:6) {

  years <- year(x = range(as.Date(x = dates[, i]), na.rm = TRUE))
  years <- paste0(years[1], "-", years[2])
  
  range <- rbind(range, years) }

range <- rbind(range, "2001-2010", "2014-2019", "2001-2019")

data <- NULL

for (i in 2:9) {
  
  name <- gsub(pattern = "_", replacement = " ", colnames(x = lrd)[i])
  
  case <- sum(lrd[, i] > 0, na.rm = TRUE)
  
  cont <- sum(lrd[, i] == 0, na.rm = TRUE)
  
  prop <- round(x = case/cont * 100, digits = 2)
  
  data <- rbind(data, c(name, case, cont, prop)) }

data <- cbind(range, data)
data <- as.data.frame(x = data[c(2:5, 1, 6:8), c(2, 1, 3:5)])

write.table(x = data, file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/1_The_proportion_of_lipid_regulating_drug_users.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = c("Source", "Dates", "Cases", "Controls", "Proportion"))

data <- data.frame(Condition = data[6:7, ]$V1, Date_range = data[6:7, ]$V2, Value = as.numeric(x = data[6:7, ]$V5))

f1 <- ggplot(data = data, aes(x = Date_range, y = Value, colour = Condition, fill = Condition)) + geom_col(width = 0.5) + 
  labs(y = "Proportion of lipid-regulating drug users (%)") + 
  theme_classic(base_size = 18) + theme(axis.title.y = element_blank()) + 
  scale_colour_manual(values = c("#3399CC", "#99FFFF")) + scale_fill_manual(values = c("#3399CC", "#99FFFF")) + scale_y_continuous(breaks = pretty_breaks()) + coord_flip() + guides(colour = "none", fill = "none")

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/1_The_proportion_of_lipid_regulating_drug_users.pdf", plot = f1, width = 11, height = 8.5)

#### Co-variate analysis ####

cov1 <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Output/Lipid-regulating drugs/Lipid-regulating_drugs_subject.csv", row.names = 1)
cov2 <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Output/Lipid-regulating drugs/Lipid-regulating_drugs_expression.csv", row.names = 1)
cov <- rbind(cov1, cov2)

colnames(x = cov) = c("Baseline Pr(>|z|)", "End line Pr(>|z|)")

write.csv(x = cov[order(cov$`End line Pr(>|z|)`), ], file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/3_The_covariate_effects_on_lipid_dysregulation.csv", quote = FALSE, row.names = TRUE)

data <- data.frame(Condition = c(rep(x = "Baseline", nrow(x = cov)), rep(x = "Follow-up", nrow(x = cov))), Covariate = c(rownames(x = cov), rownames(x = cov)), Value = c(-log10(x = cov$`Baseline Pr(>|z|)`), -log10(x = cov$`End line Pr(>|z|)`)), Group = c(rep(x = "Subject", nrow(x = cov1)), rep(x = "Expression", nrow(x = cov2)), rep(x = "Subject", nrow(x = cov1)), rep(x = "Expression", nrow(x = cov2))))
data <- data[!data$Covariate=="Zygosity1", ]
data <- data[data$Value>0, ]

f3 <- ggplot(data = data, aes(x = reorder(x = Covariate, -Value), y = Value, fill = Condition)) + geom_bar(stat = "identity", position = position_dodge()) + 
  geom_hline(yintercept = -log10(x = 0.05), linetype = "dashed", colour = "red") + labs(x = "Date range", y = "-log10(p-value)") + theme_classic(base_size = 18) + 
  theme(axis.text.x = element_text(angle = 315, hjust = 0), axis.title.x = element_blank()) + scale_fill_manual(values = c("#3399CC", "#99FFFF")) + scale_y_sqrt() + facet_wrap(Group ~ ., scales = "free")

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/3_The_covariate_effects_on_lipid-regulating_drug_use.pdf", plot = f3, width = 11, height = 8.5)

#### The proportion of censored participants with lipid regulation ####

lrd1 <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Clinical Outcomes/Lipid_regulating_drugs.csv")
lrd2 <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Clinical Outcomes/Lipid_regulating_drugs_censored.csv")

lrd1 <- colSums(x = lrd1[, 7:8], na.rm = TRUE)
lrd2 <- colSums(x = lrd2[, 2:3], na.rm = TRUE)
lrd <- rbind(Censored = lrd1 - lrd2, Uncensored = lrd2, Total = lrd1, Proportion = round(x = (lrd1 - lrd2)/lrd1, digits = 2))

colnames(x = lrd)[2] <- "Follow-up"

write.table(x = lrd, file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/7_The_proportion_of_censored_participants.csv", quote = FALSE, sep = ",", row.names = TRUE, col.names = NA)

lrd <- as.data.frame(x = t(lrd))

data <- data.frame(Class = rownames(x = lrd), Value = lrd$Proportion)

f7 <- ggplot(data = data, aes(x = Class, y = Value, colour = Class, fill = Class)) + 
  geom_bar(stat = "identity", width = 0.5) +
  labs(y = "Proportion of censored participants (%)") +
  theme_classic(base_size = 18) + theme(axis.title.y = element_blank()) + 
  scale_colour_manual(values = c("#3399CC", "#99FFFF")) + scale_fill_manual(values = c("#3399CC", "#99FFFF")) + scale_y_continuous(breaks = pretty_breaks()) + coord_flip() + guides(colour = "none", fill = "none")

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/7_The_proportion_of_censored_participants.pdf", plot = f7, width = 11, height = 8.5)

#### The most significant phenotype ####

library(ggpubr)
library(ggsignif)
library(rstatix)

phenotypes <- master1

phenotype1 <- phenotypes$`Apolipoprotein B in Serum`
phenotype2 <- phenotypes$`Trunk fat / limb fat`
phenotype3 <- phenotypes$`Systolic blood pressure`
phenotype4 <- phenotypes$`Frailty Index`

phenotypes_x4 <- rbind(phenotype1, phenotype2, phenotype3, phenotype4)

phenotype1$End_line <- ifelse(test = phenotype1$End_line == 1, yes = "Drug", no = "Control")
phenotypes_x4$End_line <- ifelse(test = phenotypes_x4$End_line == 1, yes = "Drug", no = "Control")

phenotype1 <- as.data.frame(x = phenotype1)
phenotypes_x4 <- as.data.frame(x = phenotypes_x4)

phenotype1$End_line <- as.factor(x = phenotype1$End_line)
phenotypes_x4$End_line <- as.factor(x = phenotypes_x4$End_line)

phenotype1 <- phenotype1[-which(x = is.na(x = phenotype1$End_line)), ]
phenotypes_x4 <- phenotypes_x4[-which(x = is.na(x = phenotypes_x4$End_line)), ]

anova <- phenotype1 %>% anova_test(Result ~ End_line)

ttest <- t_test(data = phenotype1, formula = Result ~ End_line, p.adjust.method = "bonferroni", detailed = TRUE)
ttest <- ttest %>% add_xy_position(x = "End_line")

f4.1 <- ggplot(data = phenotype1, aes(x = End_line, y = Result)) +
  geom_violin(aes(fill = End_line)) + geom_boxplot(width = 0.25) + 
  geom_signif(comparisons = list(c("Control", "Drug")), annotations = "****") +
  labs(title = "ApoB in serum at baseline by future lipid-regulating drug use", x = "Future lipid-regulating drug use", y = "Apolipoprotein B in serum at baseline (g/L)", subtitle = get_test_label(stat.test = anova[1, ], detailed = TRUE), caption = get_pwc_label(stat.test = ttest), colour = "Class") +
  scale_fill_jama() + scale_y_continuous(breaks = pretty_breaks()) +
  guides(fill = "none") +
  theme_classic(base_size = 18)

f4.2 <- ggplot(data = phenotypes_x4, aes(x = End_line, y = Result)) +
  geom_violin(aes(fill = End_line)) + geom_boxplot(width = 0.25) +
  scale_fill_jama() + scale_y_continuous(breaks = pretty_breaks()) +
  guides(fill = "none") +
  theme_classic(base_size = 18) + theme(text = element_text(size = 20), axis.title = element_blank()) +
  facet_wrap(facets = ~ Description, scales = "free")

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/4_Phenotypes_1_vs_end_line_lipid_dysregulation.pdf", plot = f4.1, width = 11, height = 8.5)
ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/4_Phenotypes_4_vs_end_line_lipid_dysregulation.pdf", plot = f4.2, width = 11, height = 8.5)

#### The most significant change in phenotype ####

phenotypes <- master2
phenotype1 <- phenotypes$`Low density lipoprotein`
phenotype2 <- phenotypes$ASCVD
phenotype3 <- phenotypes$`Body Mass Index`
phenotype4 <- phenotypes$`Fat tissue in trunk region`
phenotypes_x4 <- rbind(phenotype1, phenotype2, phenotype3, phenotype4)

phenotype1$End_line <- ifelse(test = phenotype1$End_line == 1, yes = "Drug", no = "Control")
phenotypes_x4$End_line <- ifelse(test = phenotypes_x4$End_line == 1, yes = "Drug", no = "Control")

phenotype1 <- as.data.frame(x = phenotype1)
phenotypes_x4 <- as.data.frame(x = phenotypes_x4)

phenotype1$End_line <- as.factor(x = phenotype1$End_line)
phenotypes_x4$End_line <- as.factor(x = phenotypes_x4$End_line)

anova <- phenotype1 %>% anova_test(resids ~ End_line)

ttest <- t_test(data = phenotype1, formula = resids ~ End_line, p.adjust.method = "bonferroni", detailed = TRUE)
ttest <- ttest %>% add_xy_position(x = "End_line")

f4.3 <- ggplot(data = phenotype1, aes(x = timdiff, y = rawdiff, colour = End_line, fill = End_line)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  labs(title = "Change in low-density lipoprotein by future lipid-regulating drug use", x = "Future lipid-regulating drug use", y = "Change in low-density lipoprotein (mmol/L)", subtitle = get_test_label(stat.test = anova[1, ], detailed = TRUE), caption = get_pwc_label(stat.test = ttest), colour = "Class") +
  scale_color_jama() + scale_fill_jama() + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks()) +
  guides(color = "none", fill = "none") +
  theme_classic(base_size = 18)

f4.4 <- ggplot(data = phenotypes_x4, aes(x = timdiff, y = rawdiff, colour = End_line, fill = End_line)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  labs(x = "Time difference between visits (days)", y = "Raw difference in result between visits") +
  scale_color_jama() + scale_fill_jama() + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks()) +
  guides(color = "none", fill = "none") +
  theme_classic(base_size = 18) + facet_wrap(facets = ~ Description, scales = "free")

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/4_Phenotypes_smooth_1_vs_end_line_lipid_dysregulation.pdf", plot = f4.4, width = 11, height = 8.5)
ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/4_Phenotypes_smooth_4_vs_end_line_lipid_dysregulation.pdf", plot = f4.5, width = 11, height = 8.5)

#### The VatSat Ratio ####

vatsat <- master2$`Percentage of fat in largest visceral fat region / Percentage of fat in android region`

change <- lm(formula = rawdiff ~ timdiff + V1, data = vatsat)

predicted <- predict(object = change)

residuals <- residuals(object = change)

vatsat <- cbind(vatsat, as.data.frame(x = predicted), as.data.frame(x = residuals))
vatsat <- vatsat[-which(x = is.na(x = vatsat$End_line)), ]
vatsat$End_line <- ifelse(test = vatsat$End_line == 1, yes = "Drug", no = "Control")
vatsat$End_line <- factor(x = vatsat$End_line)

f4.5 <- ggplot(data = vatsat, aes(x = timdiff, y = rawdiff, group = End_line, colour = End_line, fill = End_line)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.6) +
  labs(x = "Time difference between visits (days)", y = "Raw difference in visceral fat: android fat") +
  scale_color_jama() + scale_fill_jama() + scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks()) +
  theme_classic(base_size = 18) + facet_wrap(~ End_line) + guides(color = "none", fill = "none")

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/4_Time_difference_vs_raw_difference_in_vat_sat_ratio.pdf", plot = f4.5, width = 11, height = 8.5)

f4.6 <- ggplot(data = vatsat, mapping = aes(x = timdiff, y = rawdiff)) +
  labs(x = "Time difference between visits (days)", y = "Raw difference in visceral fat: android fat") +
  scale_color_gsea() + scale_fill_jama() + 
  scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks()) +
  guides(color = "none") +
  theme_classic(base_size = 18) + 
  geom_segment(mapping = aes(xend = timdiff, yend = predicted), alpha = 0.5) +
  geom_point(aes(color = residuals)) +
  geom_point(mapping = aes(y = predicted), shape = 1)

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/4_Time_difference_vs_raw_difference_in_vat_sat_ratio_resids.pdf", plot = f4.6, width = 11, height = 8.5)

#### Differential expression ####

library(gprofiler2)

results_base <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/9_Baseline_deseq.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
results_post <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/9_End_line_deseq.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

sig.base <- results_base[which(x = results_base$pvalue < 0.05), ]
sig.post <- results_post[which(x = results_post$pvalue < 0.05), ]

sig.res <- subset(x = sig.post, !(GeneName %in% sig.base$GeneName))

gostres <- gost(query = sig.res$GeneName, organism = "hsapiens", ordered_query = TRUE, correction_method = "gSCS", custom_bg = results_post$GeneName, sources = c("HP", "GO", "KEGG", "REAC", "WP"))

publish_gosttable(gostres = gostres, filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/9_DESeq2_gosttable.pdf")

p <- gostplot(gostres = gostres, interactive = FALSE)

publish_gostplot(p = p, filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/9_DESeq2_gostplot.pdf", highlight_terms = c("GO:0005759", "HP:0001943", "HP:0030956", "HP:0011015", "HP:0011675", "HP:0002395", "KEGG:00020", "KEGG:00280", "KEGG:04614", "KEGG:00062", "KEGG:04934", "KEGG:00910", "WP:WP4236", "WP:WP5220", "WP:WP236"), width = 11, height = 8.5)

#### Differential expression volcano ####

library(ggrepel)

de <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Tables/9_End_line_deseq.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
de$DE <- "No"
de$DE[de$log2FoldChange > 0.25 & de$padj < 0.1] <- "Up"
de$DE[de$log2FoldChange < -0.25 & de$padj < 0.1] <- "Down"
de$delabel <- NA
de$delabel[de$DE != "No"] <- de$GeneName[de$DE != "No"]

f9 <- ggplot(data = de, aes(x = log2FoldChange, y = -log10(x = padj), col = DE, label = delabel)) +
  geom_point() + geom_text_repel() + geom_vline(xintercept = c(-0.25, 0.25), col = "black") + geom_hline(yintercept = -log10(x = 0.1), col = "black") +
  labs(x = "Log2 fold change", y = "-log10 adjusted p-value") +
  scale_color_manual(values = c("#CC0000", "#999FFF", "#0066CC")) + 
  scale_x_continuous(limits = c(-max(de$log2FoldChange), max(de$log2FoldChange))) + scale_y_continuous(trans = 'log10', limits = c(1e-02, max(-log10(x = de$padj)))) + 
  theme_classic(base_size = 18) + 
  guides(color = "none")

ggsave(filename = "~/Documents/KCL/PhD/TwinsUK/Results/Chapter 1/Figures/9_Volcano_plot.pdf", plot = f9, width = 11, height = 8.5)
