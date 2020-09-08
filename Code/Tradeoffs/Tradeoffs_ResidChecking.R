# Checking on the log-linear pattern present

# Loading necessary packages
library(tidyverse) # Data manipulation
library(vegan) # General ecology functions
library(data.table) # Fast reading of data tables
library(testthat) # Unit testing
library(fastDummies) # Fast generation of dummy variables
library(permute)
library(mvabund)

# For a given site

# To get an estimate of the nutrient treatment effect -
  # Generate residuals of model that includes both year and plot effects
  # Plot relative to treatment * year
cover.mat <- read.csv("../../Data/tradeoffs_comm.csv", stringsAsFactors = FALSE)

attr.mat <- read.csv("../../Data/tradeoffs_attr.csv", stringsAsFactors = FALSE)

sites <- "trel.us"

# Subset to a single site
com.subset <- cover.mat[attr.mat$site_code == sites,]

# Remove all zero columns
com.subset <- com.subset[,colSums(com.subset) > 0]

# Set values less than 1, but greater than 0 to 1
com.subset[com.subset > 0 & com.subset < 1] = 1
com.subset = round(com.subset)

# Standardize to log scale
com.subset <- decostand(com.subset, method = "log")

# Pull out relevant plot attributes for each site
attr.subset <- data.frame(attr.mat) %>%
  filter(site_code == sites) %>%
  mutate(year_trt = as.factor(year_trt),
         plot = as.factor(plot),
         block = as.factor(block))

# Check that these two matrices are the same size
expect_true(nrow(attr.subset) == nrow(com.subset))
expect_true(nrow(attr.subset) > 0)

mod_null <- lm.rrpp(as.matrix(com.subset) ~  year_trt + plot,
                   data = attr.subset)

resid(mod_null) %>% cbind(attr.subset) %>%
  pivot_longer(cols = -colnames(attr.subset)) %>%
  #  filter(trt != "Control") %>%
  ggplot(aes(x = as.numeric(year_trt),
             y = value,
             color = trt)) +
  geom_point() +
  facet_wrap(~name, scales = "free_y") +
  geom_smooth(method ="lm", se = FALSE)

plot(rda(com.subset ~ trt_N_num + trt_P_num + trt_K_num + plot + year_trt, attr.subset), scaling = "sites")

myrda <- rda(com.subset ~ trt_N_num + trt_P_num + trt_K_num + plot + year_trt, attr.subset)
coef(myrda)

