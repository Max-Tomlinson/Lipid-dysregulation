#### Classification #### 

library(readxl)
library(tidyr)

medication <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Medication/Prescriptions/TwinsUK_Prescribed_Medication_V1.0_2021.xlsx", sheet = 1, na = "NA")
medication <- medication[!duplicated(x = medication[, c(1, 8)], fromLast = FALSE), ]
medication$Lipid_regulating_drugs <- ifelse(test = grepl(pattern = "Lipid-Regulating Drugs", x = medication$Class, ignore.case = TRUE), yes = 1, no = 0)
medication$Lipid_regulating_drugs <- ifelse(test = medication$Source == "Q11B" & medication$Date > "2010-12-20", yes = NA, no = medication$Lipid_regulating_drugs)
medication <- na.omit(object = medication)

lrd <- spread(data = medication[, c(1, 8:9)], key = "Source", value = "Lipid_regulating_drugs")
lrd$Baseline <- ifelse(test = rowMeans(x = lrd[, 3:4], na.rm = TRUE) > 0, yes = 1, no = 0)
lrd$End_line <- ifelse(test = rowMeans(x = lrd[, c(2, 5:6)], na.rm = TRUE) > 0, yes = 1, no = 0)

for (i in which(x = lrd$Baseline == 1)) { lrd$End_line[i] <- NA }
for (i in which(x = lrd$End_line == 0)) { lrd$Baseline[i] <- 0 }

write.table(x = lrd, file = "~/Documents/KCL/PhD/TwinsUK/Data/Phenotypes/Clinical Outcomes/Lipid_regulating_drugs.csv", quote = FALSE, sep = ",", row.names = FALSE)
