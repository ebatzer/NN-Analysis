# Set working directory:
setwd("where my data is")

# Call libraries:
library(tidyr)
library(dplyr)
library(vegan)

# Read in dataframe
mydata <- read.csv("mydataframe.csv", stringsAsFactors = F, header = T)
head(mydata)

# Convert dataframe from long to wide form
speciesdata <- spread(mydata, # Dataframe to spread
                      key = species, # What column name contains species IDs
                      value = cover, # What column name contains species abundances
                      fill = 0) # Fills any cell that doesn't contain an observation with a 0 instead of NA

# Select sites to run PermANOVA PermDisp on:
site <- c("site_name") # Choose any number of site names in a vector 

# Dissimilarity index method:
dispmethod <- "bray" # Can also use "jaccard" or any other of choosing

# Number of ID columns:
idcols <- 7 # How many columns at the start of your dataset contain ID data (not species cover)

# Subsetting species data to sites of choice:
comdat <- speciesdata[speciesdata$site_code %in% sites,]

# Generating vegetation distance matrix
distmat <- vegdist(comdat[,c((idcols + 1):ncol(comdat))], # Runs on species data (not columns 1 to "idcols")
                   method = dispmethod, # Dissimilarity method (don't change here)
                   binary = FALSE) # Change to true for jaccard



### PermANOVA ###
adonis(distmat ~ # Distance matrix
         trt * year, # Model specification - what factors to compare? 
       comdat) # Data where model specification columns are contained



### Beta dispersion ###
com_betadispersion <- betadisper(distmat, # Veg dissimilarity matrix
                 comdat$year) #  Factors used to group community types


# Plot dispersion graph
plot(com_betadispersion)

# Boxplot of group dispersions
boxplot(com_betadispersion)

# Permutation test for beta dispersion
permutest(com_betadispersion, permutations = 999)

# Tukey post-hoc test of dispersion means
TukeyHSD(com_betadispersion)