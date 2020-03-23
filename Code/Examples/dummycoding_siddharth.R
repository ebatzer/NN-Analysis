library(tidyverse) # Data manipulation
library(data.table) # Fast reading of data tables
library(fastDummies) # Fast generation of dummy variables

# Community data wide format
cover.wide <- read.csv("../../Data/cover_wide_tradeoffs.csv")

# Checking that all cover values are numeric
cols <- c(colnames(cover.wide)[-c(1:6)])
cover.wide[-c(1:6)] <- apply(cover.wide[-c(1:6)], MARGIN = 2, FUN = as.numeric)

# Separating dataframe into cover and attribute matrices
attr.mat = cover.wide[,1:7] %>% column_to_rownames(var="X")
cover.vals = data.frame(cover.wide[,-c(2:7)]) %>% column_to_rownames(var="X")

# Creating new dummy variables that refer to treatment effects (3 new columns for the 4 treatments of interest)
attr.mat = dummy_columns(attr.mat, 
                         select_columns = "trt",
                         remove_first_dummy = TRUE) %>%
  # Remove NA values
  replace_na(list(trt_K = 0, trt_N = 0, trt_P = 0)) %>%
  
  # Create numeric treatment values, which are the interaction between treatment and year_trt + 1, then log transform
  mutate(trt_K_num =log(as.numeric(trt_K * year_trt) + 1),
         trt_P_num =log(as.numeric(trt_P * year_trt) + 1),
         trt_N_num =log(as.numeric(trt_N * year_trt) + 1))

# Fitting model
# Loop over sites
# Filter to different subsets of the community matrix and the attribute matrix

output_list <- list()
counter <- 1

for(chose_site in sitenames){
  
  # Subset to a single site
  com.subset <- cover.mat[attr.mat$site_code == chosen_site,]
  
  # Remove all zero columns
  com.subset <- com.subset[,colSums(com.subset) > 0]
  
  # Filter attribute matrix to a single site
  attr.subset <- data.frame(attr.mat) %>%
    filter(site_code == sites) %>%
    mutate(year_trt = as.factor(year_trt),
           plot = as.factor(plot),
           block = as.factor(block))
  
  # Fit the following model
  mod_rrpp <- lm.rrpp(com.subset ~ year_trt + block + trt_K_num + trt_P_num + trt_N_num,
                      data = attr.subset,# Attribute matrix consisting of explanatory variables
                      iter = num_iter,# Specify number of iterations
                      print.progress = FALSE)
  
  # Extract coefficients, model significance
  output = list(siteinfo = attr.subset,
                aovtable = anova(mod_rrpp, effect.type = "F",
                                 error = c("Residuals", "Residuals", "Residuals",
                                           "Residuals", "Residuals")),
                specscores = coef(mod_rrpp)) # Response coefficients

  output_list[[counter]] = output
  counter = counter + 1
}