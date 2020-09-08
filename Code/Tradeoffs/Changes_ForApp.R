ss_long_v2 <- read.csv("Coding Projects/NN-Analysis/NN-Dimensionality/Data/tradeoffs_specscores_long.csv",
                       row.names = 1)
ss_long_old <- read.csv("Coding Projects/NN-Analysis/NN-Dimensionality/Data/tradeoffs_specscores_long_old.csv", 
                        row.names = 1)

library(tidyverse)


compdf <- full_join(ss_longv2, ss_long_old, by = c("site", "Taxon" = "Species"))

write.csv(compdf, "Coding Projects/NN-Analysis/Data/tradeoffs_coef_changes.csv", row.names = FALSE)
