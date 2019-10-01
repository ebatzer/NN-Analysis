##### COMPOSITIONAL ANALYSIS
## Based on analyses from 2017 NutNet workshop and previous efforts
## JD Bakker et al
## 190802

## Partition each site separately and then compare across sites to correlate with site variables
## Focus on unfenced plots
## Compare Yr0 (pretreatment) through Yr3
## Focus on 'treatment' rather than NPK fertilizer combinations

# 0.0 Set working directory ----
# setwd(choose.dir()) # choose root folder that contains 'Taxonomic.Adjustments.function.190125.R' and subfolders for 'data' and 'graphs'

# 1.0 LOAD ITEMS -----------------------------------------------------------
# 1.1 Load packages ----
library(MASS)
library(plyr)
library(labdsv)
library(vegan)
#library(grid)
library(tidyverse)
library(geometry)
library(stringr)
library(ggfortify)
library(emmeans)
library(lme4)
library(lmerTest)
library("ggdendro")
library("betapart")
library(maps)

# 1.2 Load functions ---------------------------------------------------------------
CV <-function(x) 100 * (sd(x,na.rm=TRUE) / mean(x,na.rm=TRUE))

# Function to tally species in a site (taxon x year matrix)
spplist <- function(datafile = datafile, site_code = site_code) {
  temp <- matrify(ddply(datafile[datafile$site_code == site_code,], .(Taxon, year), summarize, N = length(max_cover > 0)))
  temp <- merge(x = temp, by.x = "row.names", all.x = TRUE,
                y = unique(datafile[,c("Taxon", "Family", "live")]), by.y = "Taxon", all.y = FALSE)
  temp <- temp[order(temp$Family, temp$Row.names), ]
  temp
}

# Function to tally species in a site (taxon x block matrix)
spplist2 <- function(datafile = datafile, site_code = site_code) {
  temp <- matrify(ddply(datafile[datafile$site_code == site_code,], .(Taxon, block), summarize, N = length(max_cover > 0)))
  temp <- merge(x = temp, by.x = "row.names", all.x = TRUE,
                y = unique(datafile[,c("Taxon", "Family", "live")]), by.y = "Taxon", all.y = FALSE)
  temp <- temp[order(temp$Family, temp$Row.names), ]
  temp
}

# Function to tally plots in a site (block x year matrix)
plotlist <- function(datafile = datafile, site_code = site_code) {
  temp <- matrify(ddply(datafile[datafile$site_code == site_code,], .(plot, year), summarize, 
                        S = length(max_cover > 0)))
  temp <- merge(x = unique(datafile[datafile$site_code == site_code, c("block", "plot", "subplot", "trt")]),
                by.x = "plot", all.x = FALSE, y = temp, by.y = "row.names", all.y = TRUE)
  temp <- temp[order(temp$block, temp$plot), ]
  temp[temp == 0] <- "."
  temp
}
getwd()
# Function to conduct site-specific taxonomic adjustments
source("Taxonomic.Adjustments.function.190802.R")

# Function to conduct stepAIC for a response variable
stepAIC.function <- function(data = data, y = y) {
  stepAIC(object = lm(data[[y]] ~ 1, data = data),
          scope = list(upper = ~ total_mass + total_mass_CV + log.S + ANN_TEMP_RANGE + MAP,
                       lower = ~ 1),
          direction = "both", k = log(nrow(data)))
}
#currently removed MAT, TEMP_VAR, TEMP_WET_Q, MAP_VAR, management, N_Dep

# Function to create biplot from PCA
pca.plot <- function(pca.object, x = x, y = y, expansion = 1, cutoff = 0.3, lab.x.title = NULL, lab.y.title = NULL) {
  PCAvalues <- pca.object$scores
  PCAloadings <- pca.object$loadings[ , c(x, y)]
  PCAloadings.max <- apply(abs(PCAloadings), 1, FUN = max)
  PCAloadings <- PCAloadings[PCAloadings.max > cutoff, ]
  PCAloadings
  ggplot(data = PCAvalues, aes_string(x = x, y = y)) +
    geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = PCAloadings[ , x] * expansion, yend = PCAloadings[ , y] * expansion),
                 arrow = arrow(length = unit(1/2, "picas")), color = "black") +
    geom_point(size = 1.5, alpha = 0.5, shape = 21, colour = "black", fill = "dark grey") +
    annotate("text", x = PCAloadings[ , x] * (1.1 * expansion), y = PCAloadings[ , y] * (1.1 * expansion), label = rownames(PCAloadings), size = 3) +
    theme_bw(base_size = 10) +
    labs(x = lab.x.title, y = lab.y.title) +
    coord_fixed()
}

# Function to create violin plots for subsets (All, Control, Pretreatment)
subset.plot <- function(data = data, x = x, y = y, lab.x.title = "Subset", lab.y.title = NULL) {
  ggplot(data = data, aes_string(x = x, y = y)) + 
    geom_violin(colour = "grey50", draw_quantiles = c(0.5)) + 
    geom_jitter(width = 0.1, height = 0, size = 1.25, shape = 21, colour = "black", fill = "dark grey") +
    annotate("text", x = subset.letters$Subset, y = subset.letters$y, label = subset.letters$letters, size = 4) +
    theme_bw(base_size = 10) + labs(x = lab.x.title, y = lab.y.title)
}

# Function to create scatter plot relating dissimilarity metric (y) to an explanatory variable (x) and species richness (lines)
scatter.plot <- function(data.actual, data.pred, x = x, y = y, S = S, lab.x.title = NULL, lab.y.title = NULL, y.lims = c(0,1)) {
  ggplot(data = data.actual, aes_string(x = x, y = y)) +
    geom_line(data = data.pred, aes(linetype = factor(S)), size = 1) +
    geom_point(aes(fill = log(S)), size = 1.25, shape = 21, colour = "black") +
    scale_fill_gradient(low = "white", high = "black", guide = FALSE) +
    guides(linetype = FALSE) +
    lims(y = y.lims) +
    theme_bw(base_size = 10) + labs(x = lab.x.title, y = lab.y.title) 
}

# Function to save graphic to standard 2x3" size
standard.ggsave <- function(filename = filename) {
  ggsave(filename, height = 1.85, width = 3, units = "in", dpi = 600)
}

# 1.3 Load data ----
original.triplet <- read.csv("../Data/full-cover-02-August-2019.csv", header = T) # 222626 x 18
triplet <- original.triplet # backup

original.site.covars <- read.csv("../Data/comb-by-plot-clim-soil-diversity-02-Aug-2019.csv") # 18520 x 90
## 180831 dataset had to be adjusted - check if fixes have been made here ##
site.covars.orig <- original.site.covars

map.dat <- read.csv("../Data/sites-02-August-2019.csv", header = TRUE)  #140 x 12

#treatment codes
trts <- data.frame(trt = c("Control", "N", "P", "K", "PK", "NK", "NP", "NPK", "Fence", "NPK+Fence"),
                   N = c("No", "Yes", "No", "No", "No", "Yes", "Yes", "Yes", "No", "Yes"),
                   P = c("No", "No", "Yes", "No", "Yes", "No", "Yes", "Yes", "No", "Yes"),
                   K = c("No", "No", "No", "Yes", "Yes", "Yes", "No", "Yes", "No", "Yes"),
                   Fence = c("No", "No", "No", "No", "No", "No", "No", "No", "Yes", "Yes"),
                   Num.Fert = c(0, 1, 1, 1, 2, 2, 2, 3, 0, 3))

# 2.0 PRE-PROCESSING OF COVER DATA ------------------------------------------------------------
# 2.1 Choose sites ----

#sevi.us: plots completely randomly assigned to treatments. Verified with S. Collins on 181024.
#Assigning plots to 5 'pseudo-blocks' (as contiguous as feasible) for analysis.
site.covars.orig$block[site.covars.orig$site_code == "sevi.us" & site.covars.orig$plot %in% c(2, 3, 4, 12, 13, 17, 21, 26)] <- 2
site.covars.orig$block[site.covars.orig$site_code == "sevi.us" & site.covars.orig$plot %in% c(5, 8, 9, 10, 15, 19, 20, 30)] <- 3
site.covars.orig$block[site.covars.orig$site_code == "sevi.us" & site.covars.orig$plot %in% c(14, 18, 23, 24, 25, 31, 34, 35)] <- 4
site.covars.orig$block[site.covars.orig$site_code == "sevi.us" & site.covars.orig$plot %in% c(28, 29, 33, 36, 37, 38, 39, 40)] <- 5
triplet$block[triplet$site_code == "sevi.us" & triplet$plot %in% c(2, 3, 4, 12, 13, 17, 21, 26)] <- 2
triplet$block[triplet$site_code == "sevi.us" & triplet$plot %in% c(5, 8, 9, 10, 15, 19, 20, 30)] <- 3
triplet$block[triplet$site_code == "sevi.us" & triplet$plot %in% c(14, 18, 23, 24, 25, 31, 34, 35)] <- 4
triplet$block[triplet$site_code == "sevi.us" & triplet$plot %in% c(28, 29, 33, 36, 37, 38, 39, 40)] <- 5

## drop observational sites and those with < 3 years of post-treatment data
siteyear <- ddply(triplet, . (site_code), summarize,
                  years = length(unique(year_trt)),
                  min_year_trt = min(unique(year_trt)),
                  max_year_trt = max(unique(year_trt)),
                  blocks = length(unique(block)))
exptsites.all <- siteyear[siteyear$years > 1 & siteyear$min_year_trt == 0 &
                            siteyear$max_year_trt >= 3, ] # select sites with 3 years of post-treatment data (n = 62)

#drop sites with too few blocks or with missing data
exptsites.sub <- exptsites.all
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "sval.no",] #2016: year_trt = 999999?
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "azi.cn",] #no Yr3 data
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "badlau.de",] #plot 5 missing Yr3 data
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "barta.us",] #no Yr3 data
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "bldr.us",] #only two blocks
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "bnch.us",] #missing two plots in Yr2
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "doane.us",] #only one block has Yr 1-3 data
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "hopl.us",] #plots 6, 16, 26 missing Yr0 data
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "lake.us",] #plot 5 missing Yr3 data
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "mcla.us",] #plot 14 missing Yr2 data
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "pape.de",] #only one block
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "saana.fi",] #half of plots not measured in Yr3
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "sava.us",] #only two blocks
exptsites.sub <- exptsites.sub[exptsites.sub$site_code != "sgs.us",] #multiple control in Yr0 but other trts missing that yr (verified with D. Blumenthal on 140718)

# 2.2 Taxonomic QC ----
data1.all <- Taxonomic.Adjustments(datafile = triplet[triplet$site_code %in% exptsites.all$site_code , ]) # 170176 x 9
# includes sites with too few blocks or missing data (most of which have not had taxonomic QC)

# 2.3 Choose which treatment(s) and years to focus on ----
data1.all <- merge(x = data1.all, all.x = TRUE, y = trts)

data1.all <- data1.all[! data1.all$trt %in% c("Fence", "NPK+Fence"), ]
data1.all <- data1.all[data1.all$year_trt %in% c(0, 1, 2, 3), ]

#skip to include all sites:
data1.sub <- data1.all[data1.all$site_code %in% exptsites.sub$site_code, ] #63186 x 14

# 2.4 Site-specific adjustments - drop individual blocks, plots ----
#skip section to include all plots and blocks

#cbgb.us: six blocks. Verified with L. Biederman on 170324. Delete blocks 4-6
data1.sub <- with(data1.sub, data1.sub[!(site_code == "cbgb.us" & block %in% c(4, 5, 6)), ])

#cdcr.us: five blocks (delete 4-5)
data1.sub <- with(data1.sub, data1.sub[!(site_code == "cdcr.us" & block %in% c(4, 5)), ])

#cdpt.us: six blocks. On 170802, J. Knops suggested grouping blocks 1, 5, 6 (different) and 2, 3, 4 (similar). Delete blocks 1, 5, 6
data1.sub <- with(data1.sub, data1.sub[!(site_code == "cdpt.us" & block %in% c(1, 5, 6)), ])

#ethass.au: add four plot-years with no living plants (11_2013, 28_2013, 28_2014, 5_2014)
temp <- with(data1.sub, data1.sub[site_code == "ethass.au" & year == "2016" & plot %in% c(5, 11, 28) , ])
temp <- temp[c(1, 2, 3, 6), ]
temp$Taxon <- "dummy"
temp$year <- c(2013, 2014, 2014, 2013)
data1.sub <- rbind(data1.sub, temp)

#kbs.us: five blocks (delete 4, 5)
data1.sub <- with(data1.sub, data1.sub[!(site_code == "kbs.us" & block %in% c(4, 5)), ])

#kilp.fi: four blocks (delete 4)
data1.sub <- with(data1.sub, data1.sub[!(site_code == "kilp.fi" & block %in% c(4)), ])

#koffler.ca: 3 control plots in each block; different subplot for plot 19 in 2014
#focus on first control plot in each block. Verified with M. Cadotte on 170324.
data1.sub <- with(data1.sub, data1.sub[!(site_code == "koffler.ca" & plot %in% c(9, 11, 17, 21, 34, 36)), ])

#marc.ar: 3 control plots in blocks 1 and 2. Verified with J. Alberti on 170324.
data1.sub <- with(data1.sub, data1.sub[!(site_code == "marc.ar" & plot %in% c(6, 8, 11, 17)), ])

#mtca.au: four blocks (delete 4). Verified with S. Prober on 170324.
data1.sub <- with(data1.sub, data1.sub[!(site_code == "mtca.au" & block %in% c(4)), ])

#sedg.us: 2 control, 2 NPK (no fences) in each block
data1.sub <- with(data1.sub, data1.sub[!(site_code == "sedg.us" & plot %in% c(7, 10, 17, 18, 27, 28)), ])

#sevi.us: 3 of 5 'pseudo-blocks' (as contiguous as feasible) - using those farthest away from one another
data1.sub <- with(data1.sub, data1.sub[!(site_code == "sevi.us" & block %in% c(2, 4)), ])

#shps.us: four blocks (delete 4?)
data1.sub <- with(data1.sub, data1.sub[!(site_code == "shps.us" & block %in% c(4)), ])

#sier.us: five blocks (delete 4-5?); no Yr0 data for blocks 4-5
data1.sub <- with(data1.sub, data1.sub[!(site_code == "sier.us" & block %in% c(4, 5)), ])

#summ.za: 3 control (no fences) in each block; one control measured only in Yr0 and another not measured in Yr2
data1.sub <- with(data1.sub, data1.sub[!(site_code == "summ.za" & plot %in% c(1, 10, 15, 16, 21, 30)), ])

#temple.us: drop extra plots in block 2
data1.sub <- with(data1.sub, data1.sub[!(site_code == "temple.us" & plot %in% c(19, 20)), ])

#ukul.za: 2 control, 2 NPK (no fences) in each block
data1.sub <- with(data1.sub, data1.sub[!(site_code == "ukul.za" & plot %in% c(8, 10, 19, 20, 25, 30)), ])

#yarra.au: four blocks (delete 4)
data1.sub <- with(data1.sub, data1.sub[!(site_code == "yarra.au" & block %in% c(4)), ])


# 2.5 Final adjustments to site and plot data -------------------------
data1 <- data1.sub #set to data1.sub for balanced subset of sites

data1$site_code <- factor(data1$site_code)
data1$year_trt <- factor(data1$year_trt, ordered = TRUE, levels = c(0, 1, 2, 3))
data1$block <- factor(data1$block)
data1$plot <- factor(data1$plot)
data1$UBI <- do.call(paste, c(data1[c("site_code", "block")], sep = "_"))
data1$UPI <- do.call(paste, c(data1[c("site_code", "block", "plot")], sep = "_"))
data1$UPYI <- do.call(paste, c(data1[c("UPI", "year")], sep = "_"))
data1$STYI <- do.call(paste, c(data1[c("site_code", "trt", "year_trt")], sep = "_"))
data1$SNYI <- do.call(paste, c(data1[c("site_code", "N", "year_trt")], sep = "_"))

sites <- factor(unique(data1$site_code)) # 49 sites
UBIs <- unique(data1$UBI) # 147 blocks
UPIs <- unique(data1$UPI) # 1176 UPIs
UPYIs <- unique(data1$UPYI) # 4704 UPYIs
STYIs <- unique(data1$STYI) # 1568 site-trt-year combinations

#Commented code here creates a summary table for each site based just on plots included in this analysis
#sites.data1 <- unique(data1$site_code)
#for(i in 1:length(sites.data1)) {
#  write.csv(spplist(datafile = data1, site_code = sites.data1[i]), file = paste(as.character(sites.data1[i]), ".csv", sep = ""))
#}

# 3.0 PRE-PROCESS SITE-LEVEL COVARIATES -------------------------------------
# subset to selected sites
site.covars.orig$UPYI <- do.call(paste, c(site.covars.orig[c("site_code", "block", "plot", "year")], sep = "_"))
site.covars.sub <- site.covars.orig[site.covars.orig$UPYI %in% UPYIs, ]

# calculate management index
site.covars <- ddply(site.covars.sub[site.covars.sub$year_trt == 0, ], .(site_code), summarize, habitat = unique(habitat),
                     country = unique(country), managed = unique(managed), burned = unique(burned), grazed = unique(grazed),
                     anthropogenic = unique(anthropogenic), habitat = unique(habitat), elevation = unique(elevation))
site.covars$management <- with(site.covars, managed + burned + grazed + anthropogenic) #combine multiple types of management
site.covars$management[site.covars$management > 1] <- 1
site.covars$management <- as.factor(site.covars$management)

# calculate productivity covariate using control and pretreatment data (not all sites have pretreatment biomass data)
pre.con.mass <- ddply(site.covars.sub[(site.covars.sub$trt == "Control" | site.covars.sub$year_trt == 0 ) , ], .(site_code),
                      summarize, Number = length(total_mass[! is.na(total_mass)]),
                      total_mass_CV = CV(total_mass), total_mass = mean(total_mass, na.rm = TRUE))

# climate variables
bioclim.vars <- ddply(site.covars.sub, .(site_code), summarize,
                      MAT = unique(MAT_v2), MAT_RANGE = unique(MAT_RANGE_v2), ISO = unique(ISO_v2),
                      TEMP_VAR = unique(TEMP_VAR_v2), MAX_TEMP = unique(MAX_TEMP_v2), MIN_TEMP = unique(MIN_TEMP_v2),
                      ANN_TEMP_RANGE = unique(ANN_TEMP_RANGE_v2), TEMP_WET_Q = unique(TEMP_WET_Q_v2),
                      TEMP_DRY_Q = unique(TEMP_DRY_Q_v2), TEMP_WARM_Q = unique(TEMP_WARM_Q_v2),
                      TEMP_COLD_Q = unique(TEMP_COLD_Q_v2), MAP = unique(MAP_v2), MAP_WET_M = unique(MAP_WET_M_v2),
                      MAP_DRY_M = unique(MAP_DRY_M_v2), MAP_VAR = unique(MAP_VAR_v2), MAP_WET_Q = unique(MAP_WET_Q_v2),
                      MAP_DRY_Q = unique(MAP_DRY_Q_v2), MAP_WARM_Q = unique(MAP_WARM_Q_v2), MAP_COLD_Q = unique(MAP_COLD_Q_v2),
                      RAIN_PET = unique(RAIN_PET), AI = unique(AI), PET = unique(PET), N_Dep = unique(N_Dep))

# calculate site-level species richness (gamma diversity)
temp <- ddply(data1, .(site_code, Taxon), summarize, max_cover = sum(max_cover))
S <- ddply(temp, .(site_code), summarize, S = length(Taxon))

# combine all site-level covariates together
site.covars <- merge(x = site.covars, y = pre.con.mass) %>%
  merge(y = bioclim.vars[ , c("site_code", "MAT", "TEMP_VAR", "ANN_TEMP_RANGE", "TEMP_WET_Q", "MAP", "MAP_VAR", "N_Dep")]) %>%
  merge(y = S)
site.covars$log.S <- log(site.covars$S)


# 4.0 DATA PROCESSING -----------------------------------------------------------
# 4.1 Create Plot-Year (UPYI) x Species matrix ----
spp.data <- matrify(data1[ , c("UPYI", "Taxon", "max_cover")]) # 4704 x 1530; species data for analysis

#add dummy variable with low cover
spp.data$dummy <- 0.001

# calculate plot-level species richness
S.data <- (rowSums(spp.data > 0) - 1)

#order by plot and year
spp.data <- spp.data[ order(rownames(spp.data)), ]

data <- unique(data1[, c("STYI", "UPYI", "UPI", "UBI", "year", "site_code", "block", "plot", "year_trt", "trt", "N", "P", "K", "Num.Fert")])
data <- data[which(data$UPYI %in% rownames(spp.data)), ]
data <- data[ order(rownames(spp.data)), ]

# 4.2 Calculate various dissimilarity matrices and extract SS ----
data.sub <- data[data$site_code %in% exptsites.sub$site_code , ]
part.distances <- data.frame()
site.summary <- data.frame()
subset.summary <- data.frame()
for(i in 1:length(sites[sites %in% unique(exptsites.sub$site_code)])) {
  temp <- c()
  site.data <- data[data$site_code == sites[sites %in% unique(exptsites.sub$site_code)][i], ]
  site.data <- site.data[ order(site.data$UPYI), ]
  site.spp.data <- spp.data[rownames(spp.data) %in% unique(site.data$UPYI), ]
  site.spp.data <- site.spp.data[ order(rownames(site.spp.data)), colSums(site.spp.data) > 0]
  
  #Abundance-based dissimilarity
  temp <- bray.part(site.spp.data)
  bray.diss <- as.matrix(temp$bray)
  prop.bray.bal.diss <- with(temp, as.matrix(bray.bal / bray))
  prop.bray.bal.diss[is.na(prop.bray.bal.diss)] <- 0
  
  #Incidence-based dissimilarity
  temp.pa <- beta.pair(decostand(site.spp.data, "pa"))
  sor.pa.diss <- as.matrix(temp.pa$beta.sor)
  prop.sim.pa.diss <- with(temp.pa, as.matrix(beta.sim / beta.sor))
  prop.sim.pa.diss[is.na(prop.sim.pa.diss)] <- 0  
  
  permanova.bray <- adonis(bray.diss ~ (block + year_trt + trt)^2, data = site.data, permutations = 9, method = "euc")
  permanova.prop.bal <- adonis(prop.bray.bal.diss ~ (block + year_trt + trt)^2, data = site.data, permutations = 9, method = "euc")
  permanova.sor <- adonis(sor.pa.diss ~ (block + year_trt + trt)^2, data = site.data, permutations = 9, method = "euc")
  permanova.prop.sim <- adonis(prop.sim.pa.diss ~ (block + year_trt + trt)^2, data = site.data, permutations = 9, method = "euc")
  
  temp.bray <- data.frame(term = c(attr(permanova.bray$terms, "term.labels"), "residual", "total"),
                          df = as.numeric(permanova.bray$aov.tab$Df),
                          SS_bray = as.double(permanova.bray$aov.tab$SumsOfSqs),
                          SS_prop_bal = as.double(permanova.prop.bal$aov.tab$SumsOfSqs),
                          SS_sor = as.double(permanova.sor$aov.tab$SumsOfSqs),
                          SS_prop_sim = as.double(permanova.prop.sim$aov.tab$SumsOfSqs))
  temp.bray$site_code <- as.character(sites[i])
  part.distances <- rbind(part.distances, temp.bray)
  print(i); print(as.character(sites[i]))
  
  # site averages (all plot-years)
  temp.site <- data.frame(site_code = as.character(sites[i]),
                     bray = mean(temp$bray),
                     bray.bal = mean(temp$bray.bal),
                     bray.gra = mean(temp$bray.gra),
                     prop_bal = mean(as.dist(prop.bray.bal.diss)),
                     sor = mean(temp.pa$beta.sor),
                     sor.turn = mean(temp.pa$beta.sim),
                     sor.nest = mean(temp.pa$beta.sne),
                     prop_sim = mean(as.dist(prop.sim.pa.diss)))
  site.summary <- rbind(site.summary, temp.site)

  # site averages (subsets of plot-years)
  #index each control plot separately
  control.plots <- unique(site.data$UPI[site.data$trt == "Control"])
  control.plot.dists <- data.frame()
  for(j in 1:length(control.plots)) {
    temp.dists <- data.frame(site_code = as.character(sites[i]),
                             UPI = control.plots[j],
                             control.bray.dist = mean(as.dist(bray.diss[row.names(bray.diss) %in% site.data$UPYI[site.data$UPI == control.plots[j]],
                                                                         colnames(bray.diss) %in% site.data$UPYI[site.data$UPI == control.plots[j]]])),
                             control.prop_bal.dist = mean(as.dist(prop.bray.bal.diss[row.names(prop.bray.bal.diss) %in% site.data$UPYI[site.data$UPI == control.plots[j]],
                                                                             colnames(prop.bray.bal.diss) %in% site.data$UPYI[site.data$UPI == control.plots[j]]])),
                             control.sor.dist = mean(as.dist(sor.pa.diss[row.names(sor.pa.diss) %in% site.data$UPYI[site.data$UPI == control.plots[j]],
                                                                                 colnames(sor.pa.diss) %in% site.data$UPYI[site.data$UPI == control.plots[j]]])),
                             control.prop_sim.dist = mean(as.dist(prop.sim.pa.diss[row.names(prop.sim.pa.diss) %in% site.data$UPYI[site.data$UPI == control.plots[j]],
                                                                                      colnames(prop.sim.pa.diss) %in% site.data$UPYI[site.data$UPI == control.plots[j]]])))
    control.plot.dists <- rbind(control.plot.dists, temp.dists)
  }
  
  #combine overall, pretreatment, and control subsets
  temp.subset <- data.frame(site_code = as.character(sites[i]),
                          All.bray = mean(as.dist(temp$bray)),
                          Pretreatment.bray = mean(as.dist(bray.diss[row.names(bray.diss) %in% site.data$UPYI[site.data$year_trt == 0],
                                                             colnames(bray.diss) %in% site.data$UPYI[site.data$year_trt == 0]])),
                          Control.bray = mean(control.plot.dists$control.bray.dist),
                          All.prop_bal = mean(as.dist(prop.bray.bal.diss)),
                          Pretreatment.prop_bal = mean(as.dist(prop.bray.bal.diss[row.names(prop.bray.bal.diss) %in% site.data$UPYI[site.data$year_trt == 0],
                                                                                  colnames(prop.bray.bal.diss) %in% site.data$UPYI[site.data$year_trt == 0]])),
                          Control.prop_bal = mean(control.plot.dists$control.prop_bal.dist),
                          All.sor = mean(as.dist(sor.pa.diss)),
                          Pretreatment.sor = mean(as.dist(sor.pa.diss[row.names(sor.pa.diss) %in% site.data$UPYI[site.data$year_trt == 0],
                                                                      colnames(sor.pa.diss) %in% site.data$UPYI[site.data$year_trt == 0]])),
                          Control.sor = mean(control.plot.dists$control.sor.dist),
                          All.prop_sim = mean(as.dist(prop.sim.pa.diss)),
                          Pretreatment.prop_sim = mean(as.dist(prop.sim.pa.diss[row.names(prop.sim.pa.diss) %in% site.data$UPYI[site.data$year_trt == 0],
                                                                                colnames(prop.sim.pa.diss) %in% site.data$UPYI[site.data$year_trt == 0]])),
                          Control.prop_sim = mean(control.plot.dists$control.prop_sim.dist))
  subset.summary <- rbind(subset.summary, temp.subset)
}
rownames(site.summary) <- site.summary$site_code
site.summary <- merge(x = site.summary, y = site.covars)
subset.summary <- merge(x = subset.summary, y = site.covars)

# 5.0 ANALYSIS -------------------------------------------------
# 5.1 Compare subsets (All vs. Pretreatment vs. Control) ----
# total distance
bray.summary <- gather(subset.summary, All.bray, Pretreatment.bray, Control.bray, key = Subset, value = bray) %>%
  separate(Subset, into = c("Subset", "measure")) %>%
  select(c(site_code, Subset, bray)) %>%
  merge(y = data.frame(Subset = c("All", "Control", "Pretreatment"),
                       Type = c("All", "Temporal", "Spatial")))
summary(lmer(bray ~ Type + (1 | site_code), data = bray.summary))
CLD(emmeans(lmer(bray ~ Type + (1 | site_code), data = bray.summary), specs = pairwise ~ Type), adjust = "none", which = 1, Letters = letters)
#Every subset different
ddply(bray.summary, .(Type, Subset), summarize, min = min(bray), mean = mean(bray), max = max(bray))

subset.letters <- data.frame(Subset = c(0.75, 1.75, 2.75),
                             y = c(0.22, 0.22, 0.22),
                             letters = c("c", "b", "a"))
subset.plot(data = bray.summary, x = "Type", y = "bray", lab.x.title = "Subset", lab.y.title = "Bray-Curtis Dissimilarity") + lims(y = c(0, 0.8))
standard.ggsave(filename = "graphs/bray.subset.png")

# presence/absence data
sor.summary <- gather(subset.summary, All.sor, Pretreatment.sor, Control.sor, key = Subset, value = sor) %>%
  separate(Subset, into = c("Subset", "measure")) %>%
  select(c(site_code, Subset, sor))  %>%
  merge(y = data.frame(Subset = c("All", "Control", "Pretreatment"),
                       Type = c("All", "Temporal", "Spatial")))
summary(lmer(sor ~ Type + (1 | site_code), data = sor.summary))
CLD(emmeans(lmer(sor ~ Type + (1 | site_code), data = sor.summary), specs = pairwise ~ Type), adjust = "none", which = 1, Letters = letters)
#Every subset different
ddply(sor.summary, .(Type, Subset), summarize, min = min(sor), mean = mean(sor), max = max(sor))

subset.letters <- data.frame(Subset = c(0.75, 1.75, 2.75),
                             y = c(0.07, 0.07, 0.07),
                             letters = c("c", "b", "a"))
subset.plot(data = sor.summary, x = "Type", y = "sor", lab.x.title = "Subset", lab.y.title = "Sorensen Dissimilarity") + lims(y = c(0, 0.6))
standard.ggsave(filename = "graphs/sor.subset.png")

# proportion of dissimilarity due to balanced variation in abundance
prop_bal.summary <- gather(subset.summary, All.prop_bal, Pretreatment.prop_bal, Control.prop_bal, key = Subset, value = prop_bal) %>%
  separate(Subset, into = c("Subset", "measure"), extra = "merge") %>%
  select(site_code, Subset, prop_bal) %>%
  merge(y = data.frame(Subset = c("All", "Control", "Pretreatment"),
                       Type = c("All", "Temporal", "Spatial")))
summary(lmer(prop_bal ~ Type + (1 | site_code), data = prop_bal.summary))
CLD(emmeans(lmer(prop_bal ~ Type + (1 | site_code), data = prop_bal.summary), specs = pairwise ~ Type), adjust = "none", which = 1, Letters = letters)
#Temporal different than others
ddply(prop_bal.summary, .(Type, Subset), summarize, min = min(prop_bal), mean = mean(prop_bal), max = max(prop_bal))

subset.letters <- data.frame(Subset = c(0.75, 1.75, 2.75),
                             y = c(58, 58, 58),
                             letters = c("b", "b", "a"))
subset.plot(data = prop_bal.summary, x = "Type", y = "prop_bal * 100", lab.x.title = "Subset", lab.y.title = "Balanced Variation (%)") + lims(y = c(20, 100))
standard.ggsave(filename = "graphs/prop_bal.subset.png")

# proportion of presence/absence dissimilarity due to species turnover
prop_sim.summary <- gather(subset.summary, All.prop_sim, Pretreatment.prop_sim, Control.prop_sim, key = Subset, value = prop_sim) %>%
  separate(Subset, into = c("Subset", "measure"), extra = "merge") %>%
  select(site_code, Subset, prop_sim) %>%
  merge(y = data.frame(Subset = c("All", "Control", "Pretreatment"),
                       Type = c("All", "Temporal", "Spatial")))
summary(lmer(prop_sim ~ Type + (1 | site_code), data = prop_sim.summary))
CLD(emmeans(lmer(prop_sim ~ Type + (1 | site_code), data = prop_sim.summary), specs = pairwise ~ Type), adjust = "none", which = 1, Letters = letters)
#Temporal different than Spatial and All
ddply(prop_sim.summary, .(Type, Subset), summarize, min = min(prop_sim), mean = mean(prop_sim), max = max(prop_sim))

subset.letters <- data.frame(Subset = c(0.75, 1.75, 2.75),
                             y = c(28, 28, 28),
                             letters = c("b", "b", "a"))
subset.plot(data = prop_sim.summary, x = "Type", y = "prop_sim * 100", lab.x.title = "Subset", lab.y.title = "Species Turnover (%)") + lims(y = c(0, 100))
standard.ggsave(filename = "graphs/prop_sim.subset.png")

# compare different metrics to one another
summary.all <- merge(x = bray.summary, y = sor.summary) %>%
  merge(y = prop_bal.summary) %>%
  merge(y = prop_sim.summary)
ggplot(data = summary.all, aes(x = bray, y = sor)) +
  geom_point() +
  facet_wrap(facets = ~ Type) +
  geom_abline(slope = 1)
ggplot(data = summary.all, aes(x = prop_bal, y = prop_sim)) +
  geom_point() +
  facet_wrap(facets = ~ Type) +
  geom_abline(slope = 1)
ggplot(data = summary.all, aes(x = sor, y = prop_sim)) +
  geom_point() +
  facet_wrap(facets = ~ Type) +
  geom_abline(slope = 1)


# 5.2 Relate mean compositional dissimilarity to covariates ----

## Table S3 (Fit model for each explanatory variable separately)
#regression - everything but management
site.summary.gather <- gather(site.summary, log.S, total_mass, total_mass_CV, MAT, TEMP_VAR, ANN_TEMP_RANGE, TEMP_WET_Q, MAP, MAP_VAR, N_Dep, key = Variable, value = X) %>%
  gather(bray, prop_bal, sor, prop_sim, key = Response, value = Y) %>%
  select(site_code, Variable, X, Response, Y)
individ.models <- dlply(site.summary.gather, .(Response, Variable), function(df) { 
  lm(Y ~ X, data = df)
})
(individ.models.P <- ldply(individ.models, function(x) round(summary(x)$coefficients[2, 4], 4)) %>%
    spread(key = Response, value = V1))

#ANOVA - management
site.summary.gather.management <- gather(site.summary, bray, prop_bal, sor, prop_sim, key = Response, value = Y) %>%
  select(site_code, management, Response, Y)
management.models <- dlply(site.summary.gather.management, .(Response), function(df) { 
  aov(Y ~ management, data = df)
})
(management.models.P <- ldply(management.models, function(x) summary.aov(x)[[1]][["Pr(>F)"]][1]))

#combine all P-values in table
tableS3 <- rbind(individ.models.P,
                 c("management", round(t(management.models.P$V1), 4)))
tableS3$row.order <- c(7, 2, 9, 10, 5, 11, 6, 8, 3, 4, 1)
tableS3[order(tableS3$row.order), c(1, 2, 5, 3, 4)]

## Bray-Curtis dissimilarity
summary(bray.res <- stepAIC.function(data = site.summary, y = "bray")) #total_mass (r^2 = 0.1134)
total_mass.grid <- with(site.summary, expand.grid(
  total_mass = seq(min(total_mass), max(total_mass), length = 20)
))
total_mass.grid$bray <- stats::predict(lm(bray ~ total_mass, data = site.summary), newdata = total_mass.grid)
scatter.plot(data.actual = site.summary, data.pred = total_mass.grid, x = "total_mass", y = "bray", S = "solid",
             lab.x.title = expression(paste("Productivity (g/m"^2, ")")), lab.y.title = "", y.lims = c(0.2, 0.8))
standard.ggsave(filename = "graphs/bray.total_mass.png")

## Sorensen dissimilarity
summary(sor.res <- stepAIC.function(data = site.summary, y = "sor")) #log.S, ANN_TEMP_RANGE, total_mass_CV (r^2 = 0.5044)
S.ANN_TEMP_RANGE.grid <- with(site.summary, expand.grid(
  ANN_TEMP_RANGE = seq(min(ANN_TEMP_RANGE), max(ANN_TEMP_RANGE), length = 20),
  S = c(20, 60, 100)
))
S.ANN_TEMP_RANGE.grid$sor <- stats::predict(lm(sor ~ log(S) + ANN_TEMP_RANGE, data = site.summary), newdata = S.ANN_TEMP_RANGE.grid)
scatter.plot(data.actual = site.summary, data.pred = S.ANN_TEMP_RANGE.grid, x = "ANN_TEMP_RANGE", y = "sor", S = S,
             lab.x.title = "Annual Temperature Range (?C)", lab.y.title = "", y.lims = c(0.2, 0.6))
standard.ggsave(filename = "graphs/sor.ANN_TEMP_RANGE.logS.png")

S.total_mass_CV.grid <- with(site.summary, expand.grid(
  total_mass_CV = seq(min(total_mass_CV), max(total_mass_CV), length = 20),
  S = c(20, 60, 40)
))
S.total_mass_CV.grid$sor <- stats::predict(lm(sor ~ log(S) + total_mass_CV, data = site.summary), newdata = S.total_mass_CV.grid)
scatter.plot(data.actual = site.summary, data.pred = S.total_mass_CV.grid, x = "total_mass_CV", y = "sor", S = S,
             lab.x.title = "Local Variation in Productivity (%)", lab.y.title = "Sorensen Dissimilarity", y.lims = c(0.2, 0.6))
standard.ggsave(filename = "graphs/sor.total_mass_CV.logS.png")

## Balanced variation (%)
summary(prop_bal.res <- stepAIC.function(data = site.summary, y = "prop_bal")) #log.S, MAP (r^2 = 0.3695)
S.MAP.grid <- with(site.summary, expand.grid(
  MAP = seq(min(MAP), max(MAP), length = 20),
  S = c(20, 60, 100)
))
S.MAP.grid$prop_bal <- stats::predict(lm(prop_bal ~ log(S) + MAP, data = site.summary), newdata = S.MAP.grid)
scatter.plot(data.actual = site.summary, data.pred = S.MAP.grid, x = "MAP", y = "prop_bal * 100", S = S,
             lab.x.title = "Mean Annual Precipitation (mm)", lab.y.title = "", y.lims = c(50, 100))
standard.ggsave(filename = "graphs/prop_bal.MAP.logS.png")

## Species turnover (%)
summary(prop_sim.res <- stepAIC.function(data = site.summary, y = "prop_sim")) #log.S, total_mass, ANN_TEMP_RANGE (r^2 = 0.7051)
S.total_mass.grid <- with(site.summary, expand.grid(
  total_mass = seq(min(total_mass), max(total_mass), length = 20),
  S = c(20, 60, 100)
))
S.total_mass.grid$prop_sim <- stats::predict(lm(prop_sim ~ log(S) + total_mass, data = site.summary), newdata = S.total_mass.grid)
scatter.plot(data.actual = site.summary, data.pred = S.total_mass.grid, x = "total_mass", y = "prop_sim * 100", S = S,
             lab.x.title = expression(paste("Productivity (g/m"^2, ")")), lab.y.title = "", y.lims = c(25, 95))
standard.ggsave(filename = "graphs/prop_sim.logS.total_mass.png")

S.ANN_TEMP_RANGE.grid$prop_sim <- stats::predict(lm(prop_sim ~ log(S) + ANN_TEMP_RANGE, data = site.summary), newdata = S.ANN_TEMP_RANGE.grid)
scatter.plot(data.actual = site.summary, data.pred = S.ANN_TEMP_RANGE.grid, x = "ANN_TEMP_RANGE", y = "prop_sim * 100", S = S,
             lab.x.title = "Annual Temperature Range (?C)", lab.y.title = "Species Turnover (%)", y.lims = c(25, 95))
standard.ggsave(filename = "graphs/prop_sim.ANN_TEMP_RANGE.logS.png")


## build single graphic that shows all variables and significant predictors
#combine coefficients from final models
coefficients <- rbind(data.frame(Coef = bray.res$coefficients, Explan = row.names(data.frame(bray.res$coefficients)), Response = "bray"),
                      data.frame(Coef = prop_bal.res$coefficients, Explan = row.names(data.frame(prop_bal.res$coefficients)), Response = "prop_bal"),
                      data.frame(Coef = sor.res$coefficients, Explan = row.names(data.frame(sor.res$coefficients)), Response = "sor"),
                      data.frame(Coef = prop_sim.res$coefficients, Explan = row.names(data.frame(prop_sim.res$coefficients)), Response = "prop_sim"))

r.squared.values <- data.frame(R2 = c(summary(bray.res)$r.squared,
                                      summary(prop_bal.res)$r.squared,
                                      summary(sor.res)$r.squared,
                                      summary(prop_sim.res)$r.squared),
                               Response = c("bray", "prop_bal", "sor", "prop_sim"))

#merge files and flag significant terms
master <- site.summary %>%
  gather(bray, sor, prop_bal, prop_sim, key = Response, value = response) %>%
  gather(log.S, total_mass, total_mass_CV, ANN_TEMP_RANGE, MAP, key = Explan, value = explan) %>%
  dplyr::select(site_code, Response, response, Explan, explan) %>%
  mutate(Response = factor(Response, ordered = TRUE, levels = c("bray", "sor", "prop_bal", "prop_sim"))) %>%
  mutate(Explan = factor(Explan, ordered = TRUE, levels = c("log.S", "total_mass", "total_mass_CV", "ANN_TEMP_RANGE", "MAP")))
sig.codes <- unique(master[,c("Response", "Explan")]) %>%
  merge(y = coefficients, all.x = TRUE, all.y = FALSE)
sig.codes$Sig.AIC <- ifelse(! is.na(sig.codes$Coef), "Sig", "NS")
sig.codes$Sig.AIC <- factor(sig.codes$Sig.AIC, ordered = TRUE, levels = c("Sig", "NS"))
master <- merge(x = master, y = sig.codes)

#adjust based on individual tests 
master$Sig <- master$Sig.AIC
master$Sig[master$Response == "bray" & master$Explan == "total_mass_CV"] <- "Sig"
master$Sig[master$Response == "bray" & master$Explan == "MAP"] <- "Sig"
master$Sig[master$Response == "sor" & master$Explan == "total_mass_CV"] <- "NS"
master$Sig[master$Response == "prop_bal" & master$Explan == "total_mass_CV"] <- "Sig"
master$Sig[master$Response == "prop_sim" & master$Explan == "ANN_TEMP_RANGE"] <- "NS"
master$response[master$Response %in% c("prop_bal", "prop_sim")] <- master$response[master$Response %in% c("prop_bal", "prop_sim")] * 100
                        
response.metrics <- c(bray = "Bray-Curtis Dissimilarity", sor = "Sorensen Dissimilarity",
                      prop_bal = "Balanced Variation (%)", prop_sim = "Species Turnover (%)")
explanatory.variables <- c(log.S = "log(S)", total_mass = "Productivity", total_mass_CV = "Local Variation",
                           ANN_TEMP_RANGE = "Annual Temp. Range", MAP = "Mean Ann. Precip.")
ggplot(data = master, aes(x = explan, y = response)) +
  geom_rect(aes(fill = as.factor(Sig.AIC)), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  stat_smooth(aes(linetype = as.factor(Sig)), method = "lm", se = FALSE, size = 1, colour = "black") +
  geom_point() +
  scale_fill_manual(values = c("white", "light grey")) +
  facet_grid(facets = Response ~ Explan, scales = "free", labeller = labeller(Response = response.metrics, Explan = explanatory.variables)) +
  theme_bw(base_size = 10) +
  labs(x = "Explanatory Variable", y = "Response Variable") + guides(linetype = FALSE, fill = FALSE)

ggsave("graphs/Figure.S2.png", width = 6.5, height = 6.5, units = "in", dpi = 600)

stop()

# Alternative plot of partial residuals
bray_dat <- site.summary %>% dplyr::select(bray, log.S, total_mass, total_mass_CV, ANN_TEMP_RANGE, MAP)

bray_runs <- rbind(data.frame(resid = resid(lm(bray ~ . - log.S, bray_dat)), resp = "bray", var = bray_dat$log.S, varname = "log.S"),
      data.frame(resid = resid(lm(bray ~ . - total_mass, bray_dat)), resp = "bray", var = bray_dat$total_mass, varname = "total_mass"),
      data.frame(resid = resid(lm(bray ~ . - total_mass_CV, bray_dat)), resp = "bray", var = bray_dat$total_mass_CV, varname = "total_mass_cv"),
      data.frame(resid = resid(lm(bray ~ . - ANN_TEMP_RANGE, bray_dat)), resp = "bray", var = bray_dat$ANN_TEMP_RANGE, varname = "ann_temp_range"),
      data.frame(resid = resid(lm(bray ~ . - MAP, bray_dat)), resp = "bray", var = bray_dat$MAP, varname = "MAP"))

sor_dat <- site.summary %>% dplyr::select(sor, log.S, total_mass, total_mass_CV, ANN_TEMP_RANGE, MAP)

sor_runs <- rbind(data.frame(resid = resid(lm(sor ~ . - log.S, sor_dat)), resp = "sor", var = sor_dat$log.S, varname = "log.S"),
                   data.frame(resid = resid(lm(sor ~ . - total_mass, sor_dat)), resp = "sor", var = sor_dat$total_mass, varname = "total_mass"),
                   data.frame(resid = resid(lm(sor ~ . - total_mass_CV, sor_dat)), resp = "sor", var = sor_dat$total_mass_CV, varname = "total_mass_cv"),
                   data.frame(resid = resid(lm(sor ~ . - ANN_TEMP_RANGE, sor_dat)), resp = "sor", var = sor_dat$ANN_TEMP_RANGE, varname = "ann_temp_range"),
                   data.frame(resid = resid(lm(sor ~ . - MAP, sor_dat)), resp = "sor", var = sor_dat$MAP, varname = "MAP"))

prop_bal_dat <- site.summary %>% dplyr::select(prop_bal, log.S, total_mass, total_mass_CV, ANN_TEMP_RANGE, MAP)

prop_bal_runs <- rbind(data.frame(resid = resid(lm(prop_bal ~ . - log.S, prop_bal_dat)), resp = "prop_bal", var = prop_bal_dat$log.S, varname = "log.S"),
                   data.frame(resid = resid(lm(prop_bal ~ . - total_mass, prop_bal_dat)), resp = "prop_bal", var = prop_bal_dat$total_mass, varname = "total_mass"),
                   data.frame(resid = resid(lm(prop_bal ~ . - total_mass_CV, prop_bal_dat)), resp = "prop_bal", var = prop_bal_dat$total_mass_CV, varname = "total_mass_cv"),
                   data.frame(resid = resid(lm(prop_bal ~ . - ANN_TEMP_RANGE, prop_bal_dat)), resp = "prop_bal", var = prop_bal_dat$ANN_TEMP_RANGE, varname = "ann_temp_range"),
                   data.frame(resid = resid(lm(prop_bal ~ . - MAP, prop_bal_dat)), resp = "prop_bal", var = prop_bal_dat$MAP, varname = "MAP"))

prop_sim_dat <- site.summary %>% dplyr::select(prop_sim, log.S, total_mass, total_mass_CV, ANN_TEMP_RANGE, MAP)

prop_sim_runs <- rbind(data.frame(resid = resid(lm(prop_sim ~ . - log.S, prop_sim_dat)), resp = "prop_sim", var = prop_sim_dat$log.S, varname = "log.S"),
                   data.frame(resid = resid(lm(prop_sim ~ . - total_mass, prop_sim_dat)), resp = "prop_sim", var = prop_sim_dat$total_mass, varname = "total_mass"),
                   data.frame(resid = resid(lm(prop_sim ~ . - total_mass_CV, prop_sim_dat)), resp = "prop_sim", var = prop_sim_dat$total_mass_CV, varname = "total_mass_cv"),
                   data.frame(resid = resid(lm(prop_sim ~ . - ANN_TEMP_RANGE, prop_sim_dat)), resp = "prop_sim", var = prop_sim_dat$ANN_TEMP_RANGE, varname = "ann_temp_range"),
                   data.frame(resid = resid(lm(prop_sim ~ . - MAP, prop_sim_dat)), resp = "prop_sim", var = prop_sim_dat$MAP, varname = "MAP"))


sigvals <- master %>% dplyr::select(Response, Explan, "Sig.AIC") %>% distinct() %>%
  mutate(Explan = case_when(Explan == "ANN_TEMP_RANGE" ~ "ann_temp_range",
                            Explan == "total_mass_CV" ~ "total_mass_cv",
                            TRUE ~ as.character(Explan))) %>%
  rename("resp" = "Response",
         "varname" = "Explan",
         "Sig" = "Sig.AIC")

bind_rows(bray_runs, sor_runs, prop_bal_runs, prop_sim_runs) %>%
  left_join(sigvals) %>% 
  mutate(resp = case_when(resp == "bray" ~ "Bray-Curtis Dissimilarity",
                   resp == "sor" ~ "Sorensen Dissimilarity",
                   resp == "prop_bal" ~ "Balanced Variation (%)",
                   resp == "prop_sim" ~ "Species Turnover (%)")) %>%
  replace_na(list(Sig = "Sig")) %>%
  ggplot(aes(x = as.numeric(var),
             y = as.numeric(resid),
             linetype = Sig)) +
  geom_rect(aes(fill = as.factor(Sig)), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_point() +
  facet_grid(resp ~ varname, scales = "free") +
  stat_smooth(method = "lm", se = FALSE, color = "black", size = 1) + 
  theme_bw(base_size = 10) +
  scale_fill_manual(values = c("white", "light grey")) +
  labs(x = "Explanatory Variable", y = "Residuals") + 
  guides(linetype = FALSE, fill = FALSE) +
  ggtitle("Partial Correlation Plots")

# 5.3 PCA of variance components, and relate PCs to covariates ----
part.distances$term <- gsub(":", ".", part.distances$term)
part.distances$term.red <- revalue(part.distances$term, c(block = "Block", block.trt = "Block.Nutrient", block.year_trt = "Block.Year",
                                                          residual = "Block.Year.Nutrient", trt = "Nutrient", year_trt = "Year", year_trt.trt = "Year.Nutrient"))
part.distances$term.abbrev <- revalue(part.distances$term, c(block = "B", block.trt = "BxN", block.year_trt = "BxY",
                                                          residual = "BxYxN", trt = "N", year_trt = "Y", year_trt.trt = "YxN"))
part.distances[part.distances < 0 ] <- 0

#total B-C compositional variance
var.comp_bray <- part.distances %>% dplyr::select(site_code, term.abbrev, SS_bray) %>% matrify() %>% dplyr::select(-total) #49 x 7
var.prop_bray <- var.comp_bray / rowSums(var.comp_bray) # calculate each variance term as proportion of total; 49 x 7

var.prop.bray.triplet <- var.prop_bray
var.prop.bray.triplet$site_code <- row.names(var.prop.bray.triplet)
var.prop.bray.triplet <- gather(var.prop.bray.triplet, B, BxN, BxY, BxYxN, N, Y, YxN, key = term.abbrev, value = prop.SS) %>%
  merge(y = unique(part.distances[ , c("term", "term.red", "term.abbrev")]), by = "term.abbrev")

ggplot(data = var.prop.bray.triplet, aes(x = factor(1), y = prop.SS, fill = factor(term.red))) +
  geom_bar(width = 1, stat = "identity") + coord_polar(theta = "y") +
  facet_wrap(facets = ~ site_code, nrow = 6) + theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.border = element_blank()) +
  labs(x = "", y = "") +
  scale_y_continuous(breaks = NULL) + scale_x_discrete(breaks = NULL) + 
  scale_fill_manual(values = c("blue", "violet", "green", "grey", "red", "yellow", "orange"), 
                    guide = guide_legend(title = "Source"))
ggsave("graphs/bray.variance.piecharts.site_code.png", height = 7, width = 6, units = "in", dpi = 600)

summary(var.bray.pca <- princomp(var.prop_bray), loadings = TRUE, cutoff = 0)
biplot(var.bray.pca)

IGA <- envfit(var.bray.pca, 
       site.summary %>% dplyr::select(managed, burned, grazed, elevation, total_mass_CV, total_mass, MAT, TEMP_VAR, ANN_TEMP_RANGE, N_Dep, S, MAP, MAP_VAR),
       permutations = 9999)

pca.plot(var.bray.pca, x = "Comp.1", y = "Comp.2", expansion = 0.5, cutoff = 0,
         lab.x.title = "Temporal PC (58%)", lab.y.title = "Spatial PC (24%)") +
  theme_bw(base_size = 20) +
  geom_segment(data = data.frame(IGA$vectors$arrows),
               aes(xend = Comp.1 / 5,
                   yend = Comp.2 / 5,
                   x = 0,
                   y = 0,
                   linetype = IGA$vectors$pvals > 0.10),
               color = "red", arrow = arrow(length = unit(0.03, "npc"))) +
  guides(linetype = FALSE) +
  ggrepel::geom_text_repel(aes(x = Comp.1 / 5,
                y = Comp.2 / 5,
                label = rownames(data.frame(IGA$vectors$arrows))),
            data = data.frame(IGA$vectors$arrows),
            size = 3,
            color = "red")

ggsave(filename = "graphs/bray.PCA.large.png", width = 6, height = 4, units = "in", dpi = 600)

# graph to match design of cluster overlay below
PCAvalues <- var.bray.pca$scores
PCAloadings <- data.frame(var.bray.pca$loadings[ , c("Comp.1", "Comp.2")])
PCAloadings$Source <- row.names(PCAloadings)
PCAloadings$Plot <- c("Yes", "Yes", "No", "No", "Yes", "Yes", "No")
expansion <- 0.25
ggplot(data = PCAvalues, aes(x = Comp.1, y = Comp.2)) +
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = Comp.1 * expansion, yend = Comp.2 * expansion),
               arrow = arrow(length = unit(1/2, "picas")), color = "black") +
  geom_point(size = 1, shape = 21, colour = "black", fill = "dark grey") +
  annotate("text", x = (PCAloadings$Comp.1[PCAloadings$Plot == "Yes"] * expansion) * 1.05, 
           y = (PCAloadings$Comp.2[PCAloadings$Plot == "Yes"] * expansion) * 1.05,
           label = PCAloadings$Source[PCAloadings$Plot == "Yes"], size = 4) +
  theme_bw(base_size = 10) +
  labs(x = "Temporal PC (58%)", y = "Spatial PC (24%)") +
  coord_fixed()
ggsave("graphs/bray.PCA.png", height = 3, width = 3.25, units = "in", dpi = 600)

var.prop1_bray <- dematrify(var.prop_bray)
colnames(var.prop1_bray) <- c("site_code", "term.abbrev", "prop.var")
(summary.bray <- ddply(var.prop1_bray, .(term.abbrev), summarize, prop.var =  round(100 * sum(prop.var) / nrow(exptsites.sub), 1)))

site.covars.bray <- merge(x = site.covars, y = var.bray.pca$scores[,1:2], by.x = "site_code", by.y = "row.names")
site.covars.bray <- dplyr::rename(site.covars.bray, bray.Temporal.PC = Comp.1, bray.Spatial.PC = Comp.2)

summary(Temporal.bray.res <- stepAIC.function(data = site.covars.bray, y = "bray.Temporal.PC")) # log.S (R2 = 0.1118)
summary(Spatial.bray.res <- stepAIC.function(data = site.covars.bray, y = "bray.Spatial.PC")) # NS

#total Sorensen compositional variance
var.comp_sor <- part.distances %>% dplyr::select(site_code, term.abbrev, SS_sor) %>% matrify() %>% dplyr::select(-total) #49 x 7
var.prop_sor <- var.comp_sor / rowSums(var.comp_sor) # calculate each variance term as proportion of total; 49 x 7
summary(var.sor.pca <- princomp(var.prop_sor), loadings = TRUE, cutoff = 0)
biplot(var.sor.pca)

pca.plot(var.sor.pca, x = "Comp.1", y = "Comp.2", expansion = 0.4, cutoff = 0,
         lab.x.title = "PC1 (56%)", lab.y.title = "PC2 (32%)")
standard.ggsave(filename = "graphs/sor.PCA.png")

var.prop1_sor <- dematrify(var.prop_sor)
colnames(var.prop1_sor) <- c("site_code", "term.abbrev", "prop.sor")
(summary.sor <- ddply(var.prop1_sor, .(term.abbrev), summarize, prop.sor = round(100 * sum(prop.sor) / nrow(exptsites.sub), 1)))

site.covars.sor <- merge(x = site.covars, y = var.sor.pca$scores[,1:2], by.x = "site_code", by.y = "row.names")
site.covars.sor <- dplyr::rename(site.covars.sor, sor.Temporal.PC = Comp.1, sor.Spatial.PC = Comp.2)

summary(Temporal.sor.res <- stepAIC.function(data = site.covars.sor, y = "sor.Temporal.PC")) # log.S (R2 = 0.1243)
summary(Spatial.sor.res <- stepAIC.function(data = site.covars.sor, y = "sor.Spatial.PC")) # NS

#proportion of compositional variance due to balanced variation in abundance
var.comp_prop.bal <- part.distances %>% dplyr::select(site_code, term.abbrev, SS_prop_bal) %>% matrify() %>% dplyr::select(-total) #49 x 7
var.prop_prop.bal <- var.comp_prop.bal / rowSums(var.comp_prop.bal) # calculate each variance term as proportion of total; 49 x 7
summary(var.prop_prop.bal.pca <- princomp(var.prop_prop.bal), loadings = TRUE, cutoff = 0)
biplot(var.prop_prop.bal.pca)

pca.plot(var.prop_prop.bal.pca, x = "Comp.1", y = "Comp.2", expansion = 0.1, cutoff = 0,
         lab.x.title = "PC1 (68%)", lab.y.title = "PC2 (17%)")
standard.ggsave(filename = "graphs/prop_bal.PCA.png")

var.prop1_prop.bal <- dematrify(var.prop_prop.bal)
colnames(var.prop1_prop.bal) <- c("site_code", "term.abbrev", "prop.bal")
(summary.prop.bal <- ddply(var.prop1_prop.bal, .(term.abbrev), summarize, prop.bal =  round(100 * sum(prop.bal) / nrow(exptsites.sub), 1)))

site.covars.prop_bal <- merge(x = site.covars, y = var.prop_prop.bal.pca$scores[,1:2], by.x = "site_code", by.y = "row.names")
site.covars.prop_bal <- dplyr::rename(site.covars.prop_bal, Prop_bal.Spatial.PC = Comp.1, Prop_bal.Temporal.PC = Comp.2)

summary(Prop_bal.Spatial.res <- stepAIC.function(data = site.covars.prop_bal, y = "Prop_bal.Spatial.PC")) # total_mass (R2 = 0.0830)
summary(Prop_bal.Temporal.res <- stepAIC.function(data = site.covars.prop_bal, y = "Prop_bal.Temporal.PC")) # log.S (R2 = 0.0825)

#proportion of compositional variance due to species turnover
var.comp_prop.sim <- part.distances %>% dplyr::select(site_code, term.abbrev, SS_prop_sim) %>% matrify() %>% dplyr::select(-total) #49 x 7
var.prop_prop.sim <- var.comp_prop.sim / rowSums(var.comp_prop.sim) # calculate each variance term as proportion of total; 49 x 7
summary(var.prop_prop.sim.pca <- princomp(var.prop_prop.sim), loadings = TRUE, cutoff = 0)
biplot(var.prop_prop.sim.pca)

pca.plot(var.prop_prop.sim.pca, x = "Comp.1", y = "Comp.2", expansion = 0.2, cutoff = 0,
         lab.x.title = "PC1 (51%)", lab.y.title = "PC2 (28%)")
standard.ggsave(filename = "graphs/prop_sim.PCA.png")

var.prop1_prop.sim <- dematrify(var.prop_prop.sim)
colnames(var.prop1_prop.sim) <- c("site_code", "term.abbrev", "prop.sim")
(summary.prop.sim <- ddply(var.prop1_prop.sim, .(term.abbrev), summarize, prop.sim = round(100 * sum(prop.sim) / nrow(exptsites.sub), 1)))

site.covars.prop_sim <- merge(x = site.covars, y = var.prop_prop.sim.pca$scores[,1:2], by.x = "site_code", by.y = "row.names")
site.covars.prop_sim <- dplyr::rename(site.covars.prop_sim, prop_sim.Spatial.PC = Comp.1, prop_sim.Temporal.PC = Comp.2)

summary(Spatial.prop_sim.res <- stepAIC.function(data = site.covars.prop_sim, y = "prop_sim.Spatial.PC")) # NS
summary(Temporal.prop_sim.res <- stepAIC.function(data = site.covars.prop_sim, y = "prop_sim.Temporal.PC")) # log.S, total_mass, MAP (R2 = 0.3496)

#combine partitioning summary information from all four metrics
summary.partitions <- merge(x = summary.bray, y = summary.sor) %>%
  merge(y = summary.prop.bal) %>%
  merge(y = summary.prop.sim)

#combine PC scores
site.covars.all <- merge(x = site.covars.bray, y = site.covars.sor) %>%
  merge(y = site.covars.prop_bal) %>%
  merge(y = site.covars.prop_sim)


# 5.4 Cluster analysis of variance components (Bray-Curtis dissimilarity) ---- 
plot(prop.hclust.bray <- hclust(dist(var.prop_bray), method = "ward.D2"))
rect.hclust(prop.hclust.bray, k = 4)
prop.groups.bray <- data_frame(site_code = prop.hclust.bray$labels, Group = cutree(prop.hclust.bray, k = 4))
cutree(prop.hclust.bray, k = 4)
group_codes.bray <- data.frame(Group = unique(prop.groups.bray$Group), GroupCode = c("G1", "G2", "G3", "G4"))
prop.groups.bray <- merge(x = prop.groups.bray, y = group_codes.bray)
prop.groups.bray <- merge(x = prop.groups.bray, y = var.prop_bray, by.x = "site_code", by.y = "row.names")

(prop.groups.bray.sum <- ddply(prop.groups.bray, .(GroupCode), summarize, Sites = length(Group), Block = mean(B), Year = mean(Y), Nutrient = mean(N),
                          Block.Year = mean(BxY), Block.Nutrient = mean(BxN), Year.Nutrient = mean(YxN), Block.Year.Nutrient = mean(BxYxN)))

group_codes.bray <- data.frame(Group = unique(prop.groups.bray$Group),
                               GroupCode = c("G1", "G2", "G3", "G4"),
                               GroupName = c("Within Blocks", "Temporal", "Multiple", "Among Blocks"))
prop.groups.bray <- merge(x = prop.groups.bray, y = group_codes.bray)
(prop.groups.bray.sum <- ddply(prop.groups.bray, .(GroupName), summarize, Sites = length(Group), Block = mean(B), Year = mean(Y), Nutrient = mean(N),
                               Block.Year = mean(BxY), Block.Nutrient = mean(BxN), Year.Nutrient = mean(YxN), Block.Year.Nutrient = mean(BxYxN)))

site.covars.bray.groups <- merge(x = site.covars.bray, y = prop.groups.bray)

## draw dendrogram
bray.dendro <- dendro_data(prop.hclust.bray, type = "rectangle")
bray.dendro$labels <- merge(x = bray.dendro$labels, by.x = "label",
                            y = prop.groups.bray[ , c("site_code", "GroupName")], by.y = "site_code")
k <- 4
rect <- aggregate(x ~ GroupName, label(bray.dendro), range)
rect <- data.frame(GroupName = rect$GroupName, rect$x)
ymax <- mean(prop.hclust.bray$height[length(prop.hclust.bray$height)-((k-2):(k-1))])

ggplot(segment(bray.dendro)) +
  geom_rect(data=rect, aes(xmin = X1 - 0.35, xmax = X2 + 0.35, ymin = -0.35, ymax = ymax, colour = GroupName), 
            fill = "grey99", size = 1) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  coord_flip() + 
  scale_y_reverse(expand=c(0.4, 0)) +
  geom_text(data = bray.dendro$labels, 
            aes(x = x, y = y, label = label), size = 2.5, vjust = 0.3, hjust = 0) +
  scale_colour_discrete() +
  labs(colour = "Group") +
  theme_dendro()
ggsave("graphs/dendrogram.png", height = 8, width = 6.5, units = "in", dpi = 600)
  

## overlay groups onto PCA
chulls <- ddply(site.covars.bray.groups, .(GroupName), function(df) df[chull(df$bray.Temporal.PC, df$bray.Spatial.PC), ])
groups.text <- ddply(site.covars.bray.groups, .(GroupName), summarize, bray.Temporal.PC = mean(bray.Temporal.PC), bray.Spatial.PC = mean(bray.Spatial.PC) + 0.015)
ggplot(data = site.covars.bray.groups, aes(x = bray.Temporal.PC, y = bray.Spatial.PC)) +
  geom_polygon(data = chulls, aes(group = GroupName), alpha = 0.15) +
  geom_point(aes(fill = GroupName), shape = 21, colour = "black", size = 1) +
  annotate("text", x = groups.text$bray.Temporal.PC, y = groups.text$bray.Spatial.PC, label = groups.text$GroupName, size = 4) +
  scale_fill_discrete(guide = FALSE) +
  theme_bw(base_size = 10) + labs(x = "Temporal PC (58%)", y = "") +
  coord_fixed()
ggsave("graphs/bray.PCA.cluster.png", height = 3, width = 3.25, units = "in", dpi = 600)

## graph each site separately, within cluster groups
site.covars.bray.groups.triplet <- site.covars.bray.groups
#site.covars.bray.groups$site_code <- row.names(site.covars.bray.groups)
site.covars.bray.groups.triplet <- gather(site.covars.bray.groups.triplet, B, BxN, BxY, BxYxN, N, Y, YxN, key = term.abbrev, value = prop.SS) %>%
  merge(y = unique(part.distances[ , c("term", "term.red", "term.abbrev")]), by = "term.abbrev")

ggplot(data = site.covars.bray.groups.triplet[site.covars.bray.groups.triplet$GroupName == "Among Blocks" , ],
       aes(x = factor(1), y = prop.SS, fill = factor(term.red))) +
  geom_bar(width = 1, stat = "identity") + coord_polar(theta = "y") +
  facet_wrap(facets = ~ site_code, ncol = 6) + theme_bw(base_size = 10) +
  theme(legend.position = NULL, panel.border = element_blank()) +
  labs(x = "", y = "") +
  scale_y_continuous(breaks = NULL) + scale_x_discrete(breaks = NULL) + 
  scale_fill_manual(values = c("blue", "violet", "green", "grey", "red", "yellow", "orange"), 
                    guide = FALSE)
ggsave("graphs/bray.variance.piecharts.site_code_AmongBlocks.png", height = 1.5, width = 6, units = "in", dpi = 600)

ggplot(data = site.covars.bray.groups.triplet[site.covars.bray.groups.triplet$GroupName == "Multiple" , ],
       aes(x = factor(1), y = prop.SS, fill = factor(term.red))) +
  geom_bar(width = 1, stat = "identity") + coord_polar(theta = "y") +
  facet_wrap(facets = ~ site_code, ncol = 6) + theme_bw(base_size = 10) +
  theme(legend.position = NULL, panel.border = element_blank()) +
  labs(x = "", y = "") +
  scale_y_continuous(breaks = NULL) + scale_x_discrete(breaks = NULL) + 
  scale_fill_manual(values = c("blue", "violet", "green", "grey", "red", "yellow", "orange"), 
                    guide = FALSE)
ggsave("graphs/bray.variance.piecharts.site_code_Multiple.png", height = 6, width = 6, units = "in", dpi = 600)

ggplot(data = site.covars.bray.groups.triplet[site.covars.bray.groups.triplet$GroupName == "Temporal" , ],
       aes(x = factor(1), y = prop.SS, fill = factor(term.red))) +
  geom_bar(width = 1, stat = "identity") + coord_polar(theta = "y") +
  facet_wrap(facets = ~ site_code, ncol = 6) + theme_bw(base_size = 10) +
  theme(legend.position = NULL, panel.border = element_blank()) +
  labs(x = "", y = "") +
  scale_y_continuous(breaks = NULL) + scale_x_discrete(breaks = NULL) + 
  scale_fill_manual(values = c("blue", "violet", "green", "grey", "red", "yellow", "orange"), 
                    guide = FALSE)
ggsave("graphs/bray.variance.piecharts.site_code_Temporal.png", height = 4.5, width = 6, units = "in", dpi = 600)

ggplot(data = site.covars.bray.groups.triplet[site.covars.bray.groups.triplet$GroupName == "Within Blocks" , ],
       aes(x = factor(1), y = prop.SS, fill = factor(term.red))) +
  geom_bar(width = 1, stat = "identity") + coord_polar(theta = "y") +
  facet_wrap(facets = ~ site_code, ncol = 6) + theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.border = element_blank()) +
  labs(x = "", y = "") +
  scale_y_continuous(breaks = NULL) + scale_x_discrete(breaks = NULL) + 
  scale_fill_manual(values = c("blue", "violet", "green", "grey", "red", "yellow", "orange"), 
                    guide = guide_legend(title = "Source"))
ggsave("graphs/bray.variance.piecharts.site_code_WithinBlocks.png", height = 4.5, width = 6, units = "in", dpi = 600)


# 5.5 PCA of dissimilarity components ----
temp <- site.summary %>%
  dplyr::rename("BC" = bray, "Sor" = sor, "Bal. Var." = prop_bal, "Spp Turnover" = prop_sim)
summary(dissim.pca <- princomp(temp[ , c("BC", "Sor", "Bal. Var.", "Spp Turnover")]), loadings = TRUE, cutoff = 0)
biplot(dissim.pca)
site.summary.pca <- data.frame(site.summary, dissim.pca$scores)
site.summary.pca <- dplyr::rename(site.summary.pca, Dissim.PC1 = Comp.1, Dissim.PC2 = Comp.2, Dissim.PC3 = Comp.3, Dissim.PC4 = Comp.4)

autoplot(dissim.pca, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 4, size = 2, alpha = 0.5) +
  theme_bw(base_size = 10) +
  labs(x = "PC1 (59%)", y = "PC2 (22%)")
ggsave(filename = "graphs/Dissim.PCA.png", height = 3, width = 5, units = "in", dpi = 600)

plot(hclust(dist(temp[ , c("BC", "Sor", "Bal. Var.", "Spp Turnover")]))) #cowi.ca different from all others

# test individual terms
#regression - everything but management
site.summary.pca2<- gather(site.summary.pca, log.S, total_mass, total_mass_CV, MAT, TEMP_VAR, ANN_TEMP_RANGE, TEMP_WET_Q, MAP, MAP_VAR, N_Dep, key = Variable, value = X) %>%
  gather(Dissim.PC1, Dissim.PC2, Dissim.PC3, Dissim.PC4, key = Response, value = Y) %>%
  select(site_code, Variable, X, Response, Y)
individ.models <- dlply(site.summary.pca2, .(Response, Variable), function(df) { 
  lm(Y ~ X, data = df)
})
(individ.models.P <- ldply(individ.models, function(x) round(summary(x)$coefficients[2, 4], 4)) %>%
    spread(key = Response, value = V1))

#ANOVA - management
site.summary.pca2.management <- gather(site.summary.pca, Dissim.PC1, Dissim.PC2, Dissim.PC3, Dissim.PC4, key = Response, value = Y) %>%
  select(site_code, management, Response, Y)
management.models <- dlply(site.summary.pca2.management, .(Response), function(df) { 
  aov(Y ~ management, data = df)
})
(management.models.P <- ldply(management.models, function(x) summary.aov(x)[[1]][["Pr(>F)"]][1]))

#combine all P-values in table
appendixS2 <- rbind(individ.models.P,
                 c("management", round(t(management.models.P$V1), 4)))
appendixS2$row.order <- c(7, 2, 9, 10, 5, 11, 6, 8, 3, 4, 1)
appendixS2[order(appendixS2$row.order), ]

#model selection
summary(Dissim.PC1.res <- stepAIC.function(data = site.summary.pca, y = "Dissim.PC1")) # log.S, ANN_TEMP_RANGE, total_mass (R2 = 0.5743)
summary(Dissim.PC2.res <- stepAIC.function(data = site.summary.pca, y = "Dissim.PC2")) # total_mass_CV, MAP, log.S (R2 = 0.3955)
summary(Dissim.PC3.res <- stepAIC.function(data = site.summary.pca, y = "Dissim.PC3")) # log.S (R2 = 0.0596)
summary(Dissim.PC4.res <- stepAIC.function(data = site.summary.pca, y = "Dissim.PC4")) # total_mass (R2 = 0.2317)


# 6.0 Supplementary Material
# 6.1 Table S1 ----
tableS1 <- site.covars[ , c("site_code", "country", "elevation", "management", "S", "total_mass", "total_mass_CV",
                            "MAT", "TEMP_VAR", "ANN_TEMP_RANGE", "TEMP_WET_Q", "MAP", "MAP_VAR")]
tableS1$elevation <- round(tableS1$elevation, 0)
tableS1$total_mass <- round(tableS1$total_mass, 0)
tableS1$total_mass_CV <- round(tableS1$total_mass_CV, 0)
tableS1$MAT <- round(tableS1$MAT, 1)
tableS1$TEMP_VAR <- round(tableS1$TEMP_VAR, 0)
tableS1$ANN_TEMP_RANGE <- round(tableS1$ANN_TEMP_RANGE, 0)
tableS1$TEMP_WET_Q <- round(tableS1$TEMP_WET_Q, 1)
tableS1$MAP <- round(tableS1$MAP, 0)
tableS1$MAP_VAR <- round(tableS1$MAP_VAR, 0)
tableS1
write.csv(tableS1, "graphs/tableS1.csv")




# 6.2 WORLD MAP (Figure S1) ----
#Create world map of sites, coded by cluster group
map.dat.red <- map.dat[map.dat$site_code %in% exptsites.sub$site_code, ]

map.dat.red <- merge(x = map.dat.red, y = site.covars.bray.groups[ , c("site_code", "GroupName")], all.y = TRUE)

WorldData <- map_data('world')
WorldData %>% filter(region != "Antarctica") -> WorldData
WorldData <- fortify(WorldData)

ggplot(data = WorldData, aes(x = long, y = lat)) +
  geom_map(map = WorldData, aes(map_id = region), size=0.5, fill = "grey", color = "grey") +
  geom_point(data = map.dat.red, aes(x = longitude, y = latitude, fill = GroupName, shape = GroupName),
             size = 2, alpha = 0.5) +
  scale_fill_discrete() +
  scale_shape_manual(values = c(22, 21, 23, 24)) + 
  #coord_map("rectangular", lat0=0, xlim=c(-180,180), ylim=c(-60, 90)) +
  labs(x = NULL, y = NULL, fill = "Group", shape = "Group") +
  scale_y_continuous(breaks=c()) +
  scale_x_continuous(breaks=c()) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
ggsave("graphs/Figure.S1.png", height = 4, width = 6.5, units = "in", dpi = 600)


# 6.3 Appendix illustrating partitioning
appendix <- data.frame(Block = c(1,1,2,2),
                       Nutrient = c("N", "C", "N", "C"),
                       Plot = c(2,6,19, 11),
                       S = c(18, 20, 19, 28))
adonis(dist(appendix$S) ~ Block + Nutrient, data = appendix, permutations = 0)
summary(aov(appendix$S ~ Block + Nutrient, data = appendix))

