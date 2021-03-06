# Loading necessary packages
library("ggplot2") # Figure plotting
library("tidyverse") # Data manipulation
library("lme4") # Linear mixed effects models
library("vegan") # General ecology functions
library("data.table") # Fast reading of data tables
library("testthat") # Unit testing
library("dtplyr") # Data.table dplyr functions
library("gridExtra")
library("JostDiv") # Exponentiated diversity indices (Hill Numbers).
library("betapart")
library("permute")
library("lmerTest")
library("RRPP")
library("emmeans")
library("gridExtra")

### Reading in Datasets

# Biomass data
biomass_data <- fread("Data/comb-by-plot-02-August-2019.csv",
                      stringsAsFactors = FALSE,
                      na.strings = c('NA','NULL'))

# Loading in dataset (Updated to Feb 22st 2019)
cover <- fread('Data/full-cover-02-August-2019.csv',
               stringsAsFactors = FALSE,
               na.strings = c('NA','NULL'))

# Plot descriptions
plot.descriptions <-  read.csv("Data/comb-by-plot-clim-soil-diversity-02-Aug-2019.csv",
                               stringsAsFactors = FALSE)

# Soil characteristics
soil.chars <- read.csv("Data/soil-nutrients-22-February-2019.csv",
                       stringsAsFactors = FALSE)

###################################################################################################
# Temporal Beta-Diversity
###################################################################################################

# Filtering to sites that have 5 years of observational data
site_selection <- plot.descriptions %>% 
  filter(trt == "Control" & year_trt < 5) %>%
  filter(!is.na(proportion_par)) %>%
  group_by(site_code, plot) %>%
  summarise(nplots = n()) %>%
  filter(nplots == 5)

com_temporal <- cover %>% 
  filter(trt == "Control" & year_trt < 5) %>%
  filter(site_code %in% site_selection$site_code) %>% 
  filter(live == 1)

cover.wide <- data.frame(com_temporal) %>%
  arrange(site_code, plot, year_trt) %>%
  select(site_code, block, plot, Taxon, max_cover, year_trt) %>%
  dcast(site_code+block+plot+year_trt ~ Taxon,
        value.var='max_cover', 
        fun.aggregate = sum,drop=T,fill=0)

temporal_turnover <- cover.wide %>% 
  select(-block, -plot, -site_code, -year_trt) %>%
  mutate_all(as.numeric)

temporal_turnover[temporal_turnover > 0] = 1

plot_desc_filtered <- plot.descriptions %>% 
  filter(trt == "Control" & year_trt < 5) %>%
  filter(!is.na(proportion_par)) %>%
  filter(site_code %in% site_selection$site_code)

com.list <- split(temporal_turnover, paste(cover.wide$site_code, cover.wide$plot))

envtl.list <- split(plot_desc_filtered$proportion_par, paste(plot_desc_filtered$site_code, plot_desc_filtered$plot))

time.list <- split(plot_desc_filtered$year_trt, paste(plot_desc_filtered$site_code, plot_desc_filtered$plot))

bmass.list <- split(plot_desc_filtered$total_mass, paste(plot_desc_filtered$site_code, plot_desc_filtered$plot))

com_beta <- lapply(com.list, beta.pair, index.family = "jaccard")

# Turnover component
jac_turn <- lapply(com_beta, `[[`, 1)

jac_turn_mat <- do.call(rbind, jac_turn) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(key = "site_code",
         value = "jac_tur",
         -rowname)

jac_turn_mantel <- lapply(FUN = as.matrix, jac_turn) %>%
  do.call(rbind,.)

# Nestedness component
jac_nest<- lapply(com_beta, `[[`, 2)

jac_nest_mat <- do.call(rbind, jac_nest) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(key = "site_code",
         value = "jac_nes",
         -rowname)

jac_nest_mantel <- lapply(FUN = as.matrix, jac_nest) %>%
  do.call(rbind,.)

# Total turnover component

total_jac <- lapply(com_beta, `[[`, 3)

jac_mat_total <- do.call(rbind, total_jac) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(key = "site_code",
         value = "jac_dis",
         -rowname)

jac_tot_mantel <- lapply(FUN = as.matrix, total_jac) %>%
  do.call(rbind,.)

# Light heterogeneity
light_het <- lapply(envtl.list, dist, "euclidean")

light_mat <- do.call(rbind, light_het) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(key = "site_code",
         value = "lighthet",
         -rowname)

# Light heterogeneity
time_dif <- lapply(time.list, dist, "euclidean")

time_mat <- do.call(rbind, time_dif) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(key = "site_code",
         value = "timedif",
         -rowname)


# Biomass heterogeneity
bmass_dif <- lapply(bmass.list, dist, "euclidean")

bmass_mat <- do.call(rbind, bmass_dif) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(key = "site_code",
         value = "bmass_dif",
         -rowname)

# Plotting components
jac_df <- left_join(jac_turn_mat, light_mat) %>%
  left_join(jac_mat_total) %>%
  left_join(jac_nest_mat) %>%
  left_join(time_mat) %>%
  left_join(bmass_mat) %>%
  select(-site_code) %>%
  rename("site_code" = "rowname") %>%
  arrange(site_code)

jac_df %>%
  ggplot(aes(x = sqrt(lighthet),
             color = as.factor(timedif))) +
  stat_smooth(aes(y = jac_tur),
              method = "lm",
              color = "black",
              linetype = 2) +
  stat_smooth(aes(y = jac_dis), se = FALSE, method = "lm") +
  geom_point(aes(y = jac_tur)) + 
  stat_smooth(aes(y = jac_nes),
              method = "lm",
              color = "black",
              linetype = 3) +
  stat_smooth(aes(y = jac_dis),
              method = "lm",
              color = "black",
              linetype = 1) +
  xlim(0,.75) +
  ylim(0,1) +
  xlab("Delta Light Availability") +
  ylab("Jaccard Turnover") 

CTRL <- how(blocks = paste(jac_df$site_code, jac_df$plot))
rows <- c(1:nrow(jac_df))

perm_Mantel <- function(com, env, iter, rows, CTRL, type = "two-sided"){
  cor_obs <- cor(com, env)
  cor_stat <- c()
  for(i in 1:iter){
    perm <- shuffle(rows, control = CTRL)
    cor_stat[i] <- cor(com, env[perm])
  }
  
  if(type == "two-sided"){
    p_val <- 1 - (sum(abs(cor_obs) >= abs(cor_stat)) / (iter + 1))
  }else{
    print("I only do two-sided tests right now...\n Specify that type = 'two-sided'")
  }
  
  return(p_val)
}

perm_Mantel(com = jac_df$jac_dis,
            env = sqrt(jac_df$lighthet),
            iter = 999, 
            rows = rows, 
            CTRL = CTRL)

perm_Mantel(com = jac_df$jac_tur,
            env = sqrt(jac_df$lighthet),
            iter = 999, 
            rows = rows, 
            CTRL = CTRL)

perm_Mantel(com = jac_df$jac_nes,
            env = sqrt(jac_df$lighthet),
            iter = 999, 
            rows = rows, 
            CTRL = CTRL)

newdf <- expand.grid(lighthet = sqrt(seq(0,1, by = .01)),
                     timedif = factor(c(1,2,3,4)),
                     site_code = unique(jac_df$site_code))

mod_re_jtot <- lmer(jac_dis ~ sqrt(lighthet) + (1 + sqrt(lighthet) | site_code), data = jac_df)
mod_re_jtur <- lmer(jac_tur ~ sqrt(lighthet) + (1 + sqrt(lighthet) | site_code), data = jac_df)
mod_re_jnes <- lmer(jac_nes ~ sqrt(lighthet) + (1 + sqrt(lighthet) | site_code), data = jac_df)

anova(mod_re_jtur)

pred_dat <- newdf %>%
  bind_cols(list(jac_tot = predict(mod_re_jtot, newdf),
                 jac_tur = predict(mod_re_jtur, newdf),
                 jac_nes = predict(mod_re_jnes, newdf))) %>%
  gather(key = "source",
         value = "dis", 
         -site_code, - lighthet, -timedif)

pred_dat %>%
  ggplot(aes(x = lighthet ^ 2,
             y = dis,
             color = site_code)) +
  geom_line() + 
  guides(color = FALSE) +
  facet_wrap(~source, nrow = 1)

# Generating predicted mean responses

tot_pred <- emmeans(mod_re_jtot, var = "~lighthet", specs = "lighthet", at = list(lighthet = seq(0,1, by = .01))) %>%
  as.data.frame(.) %>%
  rename("est_tot" = "emmean",
         "tot_lower" = "lower.CL" ,
         "tot_higher" = "upper.CL") %>%
  select(-SE, -df)

tur_pred <- emmeans(mod_re_jtur, var = "~lighthet", specs = "lighthet", at = list(lighthet = seq(0,1, by = .01))) %>%
  as.data.frame(.) %>%
  rename("est_tur" = "emmean",
         "tur_lower" = "lower.CL" ,
         "tur_higher" = "upper.CL") %>%
  select(-SE, -df)

nes_pred <- emmeans(mod_re_jnes, var = "~lighthet", specs = "lighthet", at = list(lighthet = seq(0,1, by = .01)))%>%
  as.data.frame(.) %>%
  rename("est_nes" = "emmean",
         "nes_lower" = "lower.CL" ,
         "nes_higher" = "upper.CL") %>%
  select(-SE, -df)

# Plotting predicted values

temp_jac_plot <- left_join(tot_pred, tur_pred, by = "lighthet") %>%
  left_join(nes_pred, by = "lighthet") %>%
  ggplot(aes(x = lighthet)) +
  geom_line(aes(y = est_tot, linetype = factor(1)), size = 1.5) +
  geom_ribbon(aes(ymin = tot_lower,
                  ymax = tot_higher),
              alpha = .2) +
  geom_line(aes(y = est_tur, linetype = factor(2)), size = 1.5) +
  geom_ribbon(aes(ymin = tur_lower,
                  ymax = tur_higher),
              alpha = .2) +
  geom_line(aes(y = est_nes, linetype = factor(4)), size = 1.5) +
  geom_ribbon(aes(ymin = nes_lower,
                  ymax = nes_higher),
              alpha = .2) +
  xlim(0, 1) +
  ylim(0, .75) +
  xlab("Change in Proportion Light Interception") +
  ylab("Total Dissimilarity") + 
  scale_linetype_manual(values = c(1, 2, 4),
                        labels = c("Total", "Turnover", "Nestedness")) +
  ggtitle("Temporal Jaccard Dissimilarity") +
  labs(linetype = "Jaccard Dissimilarity") +
  theme_bw() + 
  theme(legend.position="bottom")

# Saving everything together

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(temp_jac_plot)

p3 <- arrangeGrob(arrangeGrob(spat_jac_plot + theme(legend.position="none"),
                              temp_jac_plot + theme(legend.position="none"),
                              nrow=1),
                  mylegend, 
                  nrow=2,heights=c(10, 1))

ggsave("../../Figures/DivTurnover_JacPlot.pdf",
       p3,
       height = 6,
       width = 8)
