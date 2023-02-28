#### Censoring #### 

class <- read.csv(file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Clinical Outcomes/Lipid_regulating_drugs.csv")
class$Family <- substr(x = class$ParticipantID, start = 1, stop = nchar(x = class$ParticipantID) - 1)

df <- 1:nrow(x = class)

for (j in 7:8) {
  
  df2 <- NULL
  
  for (i in unique(x = class$Family)) {
    
    twinpair <- class[which(x = class$Family == i), ]
    
    if (nrow(x = twinpair) == 2) {
      
      if (is.na(x = twinpair[1, j]) | is.na(x = twinpair[2, j])) {
        
        twinpair <- twinpair }
      
      else if (twinpair[1, j] == 1){twinpair[2, j] <- NA }
      else if (twinpair[2, j] == 1){twinpair[1, j] <- NA }
      else if (twinpair[1, j] == 0 & twinpair[2, j] == 0){ twinpair[2, j] <- NA } }
  
    df2 <- rbind(df2, twinpair) }
  
  df2 <- as.data.frame(x = df2[, j])
  df <- cbind(df, df2) }

df[, 1] <- class$ParticipantID

colnames(x = df) <- colnames(x = class)[c(1, 7:8)]

write.csv(x = df, file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Clinical Outcomes/Lipid_regulating_drugs_censored.csv", row.names = FALSE)