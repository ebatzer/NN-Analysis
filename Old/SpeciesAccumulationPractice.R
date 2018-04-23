# Generating species accumulation curves:
# How many new species are we seeing in each dataset over time?
# If we compare the total number of species between control and nutrient-added datasets across single and multiple years,
# Do we end up seeing a decrease in the effect size of nutrient enrichment treatment?
# Is this the same if we look at total richness and abundance-weighted diversity?
# Generic species accumulation function:
# Sequentially within the dataset, do we see more species in nutrient-added vs. control communities?
# I.e. do they experience more species over the duration of the experiment?

setwd("C:/Users/ebatzer/Box Sync/Eviner lab shared/Evan/Research Projects/NutNet/NutNetData")

library(vegan)

# Reading in all living species dataset:

livedata <- read.csv("nn_livedata.csv", stringsAsFactors = F, header=T)

head(livedata)

# What years and treatments to look at?

livedata <- livedata[livedata$year > 2008, ]
livedata <- livedata[livedata$trt %in% c("Control", "N", "P", "K", "NPK"),]


# What else to think about:

# Abundance-weighted diversity loss?
# Proportional correction to diversity value?

idcols <- 6

collector_SAR <- list()

counter <- 1

plotindex <- 1

# For each site in the dataset{
for(siteindex in unique(livedata$site_code)){
  
  site_subset <- livedata[livedata$site_code == siteindex,]
  
  for(plotindex in unique(site_subset$plot)){
    
    if(sum(site_subset$plot == plotindex) > 1){
      
      specmat <- site_subset[site_subset$plot == plotindex, c((idcols + 1):ncol(site_subset))]
      
      
      collector_SAR[[counter]] <- list(
        
        siteid = data.frame(
        site_name = siteindex,
        block = unique(site_subset$block[site_subset$plot == plotindex]),
        trt = unique(site_subset$trt[site_subset$plot == plotindex]),
        plot = plotindex),
        
        accum = data.frame(
          sites = specaccum(specmat,
                            method="collector")$sites,
          richness = specaccum(specmat,
                               method="collector")$richness)
      )
      
      counter <- counter + 1
      
    }

  }

}


organizeSamples = function(x){
  siteid = x$siteid
  accum = x$accum
  output = data.frame(key=siteid,richness=accum[,2],sites=accum[,1])
  return(output)
}

l_new = lapply(collector_SAR, organizeSamples)

samples = dplyr::bind_rows(l_new)

library(ggplot2)

p1 <- ggplot(aes(x = sites, y= richness, color = key.trt, group = key.plot), data=samples)
p1 + geom_point() + geom_line() + facet_wrap(~key.site_name)

library(dplyr)

cntrlmeans <- data.frame(group_by(samples, key.site_name, key.block, key.trt, sites) %>% filter(key.trt == "Control") %>% summarise(cntrlmean = mean(richness)))

findeffect = function(x){
  trtrich <- x[5]
  trtsites <- x[6]
  cntrlrich <- cntrlmeans$cntrlmean[cntrlmeans$key.site_name == x[1] &
                                    cntrlmeans$key.block == x[2] &
                                    cntrlmeans$sites == x[6]]
  
  output <- data.frame(site = x[1],
             block = x[2],
             sites = trtsites,
             trt = x[3],
             plot = x[4],
             trtrich = trtrich,
             cntrlrich = cntrlrich,
             totalloss = cntrlrich - as.numeric(trtrich),
             proploss = (cntrlrich - as.numeric(trtrich)) / cntrlrich)
  
  return(output)
}

new_list <- apply(samples, 1, findeffect)
samples2 = dplyr::bind_rows(new_list)

p1 <- ggplot(aes(x = as.numeric(sites), y= -totalloss, color = trt, group = plot), data=samples2)
p1 + geom_point() + geom_line() + facet_wrap(~site)
p1 <- ggplot(aes(x = as.numeric(sites), y= -proploss, color = trt, group = plot), data=samples2)
p1 + geom_point() + geom_line() + facet_wrap(~site)


