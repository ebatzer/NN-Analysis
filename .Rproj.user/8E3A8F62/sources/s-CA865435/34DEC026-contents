
# Unconstrained ordination

```{r}
# Calculating the shape of the ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

pretty_NMDS <- function(sp.mat, env.mat, group){
  
  nmds_run = metaMDS(sp.mat,
                     try = 50,
                     trymax = 100,
                     trace = FALSE,
                     distance = "horn"
  )
  
  groupvar <- subset(env.mat, select = group)
  names(groupvar) = "group"
  NMDS = bind_cols(data.frame(nmds_run$points), groupvar)
  NMDS$group = as.factor(NMDS$group)
  NMDS.mean=aggregate(NMDS[,1:2], list(group=NMDS$group), mean)
  
  df_ell <- data.frame()
  for(g in levels(NMDS$group)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                     veganCovEllipse(cov.wt(cbind(MDS1,MDS2),
                                                                            wt=rep(1/length(MDS1),length(MDS1)))$cov,
                                                                     center=c(mean(MDS1),mean(MDS2)))))
                                  ,group=g))
  }
  
  # Adding other details for facetting
  NMDS = bind_cols(NMDS, env.mat)
  
  NMDS.mean = left_join(NMDS.mean, bind_cols(group = groupvar, env.mat[,1:7]),
                        by = "group")
  
  NMDS.mean = NMDS.mean[!duplicated(NMDS.mean[,2:2]),]
  
  return(list(NMDS, NMDS.mean, df_ell, nmds_run))
}

```

```{r}
siteName <- "hopl.us"
distmethod <- "bray"
trts <- c("Control", "N", "P", "K")

ord.sp.mat <- plot.sp.mat[plot.sp.mat$trt %in% trts & plot.sp.mat$site_code == siteName,]

# Combine climate data with community dataset
ord.sp.mat$year <- as.numeric(ord.sp.mat$year)
ord.sp.mat <- right_join(climate_agg[!is.na(climate_agg$tot_precip),], ord.sp.mat, by = c("year", "site_code"))

sp.mat <- ord.sp.mat[ord.sp.mat$site_code == siteName, -(1:10)]
sp.mat <- wisconsin(sp.mat)
env.mat <- ord.sp.mat[ord.sp.mat$site_code == siteName, 1:10]

NMDS_output <- pretty_NMDS(sp.mat, env.mat, "year_trt")

NMDS = NMDS_output[[1]]
NMDS.mean = NMDS_output[[2]]
df_ell = NMDS_output[[3]]
nmds_run = NMDS_output[[4]]

# Producing figure
ggplot(data = NMDS, aes(MDS1, MDS2)) + 
  geom_point(aes(color = NMDS$trt), size = 3, alpha = 1) + 
  geom_path(data = df_ell, aes(x=MDS1, y=MDS2, group = group), size=1, linetype=1) +
  geom_text(data = NMDS.mean, aes(label = NMDS.mean$group)) + 
  theme_bw()+ 
  geom_text(data = data.frame(MDS1 = .5, MDS2 = -1.05), 
            aes(label = paste("Stress:", round(nmds_run$stress, 2))),
            fontface = "bold")
```



