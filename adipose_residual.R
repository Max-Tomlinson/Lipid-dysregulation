#### Gene Expression Association Study  ####

# Set Working Directory
setwd(dir = "/scratch/prj/dtr/Groups_WorkSpace/KerrinSmall/Max/EuroBATS/data")

# Extract Command Line Arguments
args <- commandArgs(trailingOnly = TRUE)

# Loading Packages
library(lme4)
library(parallel)

# Data Input
genes <- read.table(file = args[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
covar <- read.table(file = args[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pheno <- read.table(file = args[3], header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set Row Names and Matrix Transpose
row.names(x = genes) <- genes$gene
genes <- t(x = genes[, 7:ncol(x = genes)])

# Normality Test, Column Names and Coerce to Data Frame
shapiro_results <- NULL
for (i in 1:ncol(x = genes)){
  shapiro_gene <- shapiro.test(x = genes[, i])
  shapiro_results <- rbind(shapiro_results, shapiro_gene)
}

cols <- colnames(x = genes)
shapiro_results <- as.data.frame(x = shapiro_results, row.names = cols)[, 1:2]

# Abnormal and Set Row Names
abnormal <- shapiro_results$p.value < 0.05 / ncol(x = genes)
abnormal_genes <- as.vector(x = which(x = abnormal))
abnormal_genes <- row.names(x = shapiro_results[abnormal_genes, ])

# Value Matching and Ordering Permutation
genes <- genes[, -match(x = abnormal_genes, table = colnames(x = genes))]
genes <- genes[order(as.numeric(x = row.names(x = genes))), ]

# Index
covar <- covar[, c("TwinID", "AGE", "BMI", "INSERT_SIZE_MEDIAN", "GC_mean",
                   "PrimerIndex", "date", "Zygosity", "Family", "FatBatch")]

# Merge Data Frames
pheno <- merge(x = pheno, y = covar, by = 1)

# Scaling and Centering of Matrix-like Objects
pheno[, 2:(ncol(x = pheno)-5)] <- lapply(X = pheno[, 2:(ncol(x = pheno)-5)], FUN = scale)

# Which are TRUE?
pheno <- pheno[-which(x = is.na(x = pheno[, args[4]])), ]

# Value Matching
ind <- match(x = row.names(x = genes), table = pheno$ParticipantID)
genes <- genes[which(x = ind > 0), ]
ind <- match(x = pheno$ParticipantID, table = row.names(x = genes))
pheno <- pheno[which(x = ind > 0), ]

#### Run GEAS ####

# Parallel linear mixed-effect models
results.GEAS  <- mclapply(X = 1:ncol(x = genes), FUN = function(i){
  # Select data
  ctmp <- colnames(x = genes)[i]
  tmp <- data.frame(genes[, i], pheno[, args[4]], pheno[, c((ncol(x = pheno)-8):ncol(x = pheno))])
  names(x = tmp)[1:2] <- c(ctmp, "Phenotype")[1:2]
  formula1 <- paste0("Phenotype ~ ", ctmp, " + AGE + BMI + INSERT_SIZE_MEDIAN + GC_mean + (1|PrimerIndex) + (1|date) + (1|Zygosity) + (1|Family) + (1|FatBatch)")
  formula2 <- paste0("Phenotype ~ AGE + BMI + INSERT_SIZE_MEDIAN + GC_mean + (1|PrimerIndex) + (1|date) + (1|Zygosity) + (1|Family) + (1|FatBatch)")
  # Run linear mixed-effect model
  results.GEAS <- tryCatch({
    fm1 = fm2 = NULL
    fm1 <- lmer(formula = formula1, data = tmp, REML = F)
    fm2 <- lmer(formula = formula2, data = tmp, REML = F)
    p <- anova(fm1, fm2)$Pr[2]
    list(gene = ctmp, fm1 = fm1, fm2 = fm2, p = p, status = "OK", message = as.matrix(x = c("None", "None"), nrow = 2), n = dim(x = tmp)[1])
  # Warnings
  }, warning = function(w){
    fm1 <- lmer(formula = formula1, data = tmp, REML = F)
    fm2 <- lmer(formula = formula2, data = tmp, REML = F)
    p <- anova(fm1, fm2)$Pr[2]
    list(gene = ctmp, fm1 = fm1, fm2 = fm2, p = p, status = "Warning", message = w, n = dim(x = tmp)[1])
  # Errors
  }, error = function(e){
    list(gene = ctmp, fm1 = NA, fm2 = NA, p = NA, status = "Error", message = e, n = dim(x = tmp)[1])
  })
  results.GEAS
}, mc.cores = 12)

#### Extract results ####

# For loop
results.error <- NULL
for (i in 1:length(x = results.GEAS)){
  if (results.GEAS[[i]]$status == "Error"){
    col <- (x = which(x = colnames(x = genes) == as.character(x = i$gene)))
    results.error <- c(results.error, col)
    print(x = i)
  } 
}

# Condition Handling and Recovery
tryCatch({results.GEAS <- results.GEAS[-(results.error)]}, error = function(e) print(x = "No linear mixed-effects models with error status"))

# Apply a Function over a List or Vector
gene    <- sapply(X = results.GEAS, FUN = function(r){ r$gene })
n       <- sapply(X = results.GEAS, FUN = function(r){ r$n })
p       <- sapply(X = results.GEAS, FUN = function(r){ r$p })
status  <- sapply(X = results.GEAS, FUN = function(r){ r$status })
message <- sapply(X = results.GEAS, FUN = function(r){ r$message })
beta    <- sapply(X = results.GEAS, FUN = function(r){ coef(object = summary(object = r$fm1))[2] })
se      <- sapply(X = results.GEAS, FUN = function(r){ coef(object = summary(object = r$fm1))[8] })

# Data Frames and Matrices
results.GEAS <- data.frame(cbind(gene, n, beta, se, p, status, message[1, ], message[2, ]))
results.GEAS <- as.matrix(x = results.GEAS)

# Data Output
write.table(x = results.GEAS, file = paste0("../output/adipose/residual/", as.character(x = gsub(pattern = "\\.", replacement = "_", x = args[4])), format(x = Sys.Date(), "_%Y%m%d"), ".txt"), row.names = FALSE, sep = "\t", quote = FALSE)

