
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

# Filtering to sites that have 30 plots, 
site_selection <- plot.descriptions %>% 
  filter(year_trt == 0) %>%
  filter(!is.na(proportion_par)) %>%
  group_by(site_code) %>%
  summarise(nplots = n()) %>%
  filter(nplots == 30) %>%
  filter(!(site_code == "ethass.au" | site_code == "hopl.us"))

com_spatial <- cover %>% 
  filter(year_trt == 0) %>%
  filter(site_code %in% site_selection$site_code) %>% 
  filter(live == 1)

cover.wide <- data.frame(com_spatial) %>%
  select(site_code, block, plot, Taxon, max_cover) %>%
  dcast(site_code+block+plot ~ Taxon,
        value.var='max_cover', 
        fun.aggregate = sum,drop=T,fill=0)

temporal_turnover <- cover.wide %>% 
  select(-block, -plot, -site_code) %>%
  mutate_all(as.numeric)

temporal_turnover[temporal_turnover > 0] = 1

plot_desc_filtered <- plot.descriptions %>% 
  filter(year_trt == 0) %>%
  filter(!is.na(proportion_par)) %>%
  filter(site_code %in% site_selection$site_code)

# Cleaning errors with Barta Brothers site
plot_desc_filtered$Ground_PAR[plot_desc_filtered$site_code == "barta.us" &
                                plot_desc_filtered$proportion_par == max(plot_desc_filtered$proportion_par)] = 493.5

plot_desc_filtered$proportion_par[plot_desc_filtered$site_code == "barta.us" &
                                    plot_desc_filtered$proportion_par == max(plot_desc_filtered$proportion_par)] =
  plot_desc_filtered$Ground_PAR[plot_desc_filtered$site_code == "barta.us" &
                                  plot_desc_filtered$proportion_par == max(plot_desc_filtered$proportion_par)]/ 
  plot_desc_filtered$Ambient_PAR[plot_desc_filtered$site_code == "barta.us" &
                                   plot_desc_filtered$proportion_par == max(plot_desc_filtered$proportion_par)]

com.list <- split(temporal_turnover, cover.wide$site_code)

envtl.list <- split(plot_desc_filtered$proportion_par, plot_desc_filtered$site_code)

com_beta <- lapply(com.list, beta.pair, index.family = "jaccard")

names(com_beta[[1]]) # Jaccard due to turnover, Jaccard due to nestedness, Total Jaccard

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

light_mat_mantel <- lapply(FUN = as.matrix, light_het) %>%
  do.call(rbind,.)

# Plotting components
jac_df <- left_join(jac_turn_mat, light_mat) %>%
  left_join(jac_mat_total) %>%
  left_join(jac_nest_mat) %>%
  select(-site_code) %>%
  rename("site_code" = "rowname") %>%
  arrange(site_code)
  
jac_df %>%
  ggplot(aes(x = sqrt(lighthet),
             color = site_code)) +
  stat_smooth(aes(y = jac_tur),
              method = "lm",
              color = "black",
              linetype = 2) +
  stat_smooth(aes(y = jac_nes),
              method = "lm",
              color = "black",
              linetype = 3) +
  stat_smooth(aes(y = jac_dis),
              method = "lm",
              color = "black",
              linetype = 1) +
  xlim(0,1) +
  ylim(0,1) +
  xlab("Delta Light Availability") +
  ylab("Jaccard Turnover") + 
  guides(color = FALSE)

# Running a bootleg Mantel test
CTRL <- how(blocks = jac_df$site_code)
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
            env = jac_df$lighthet,
            iter = 999, 
            rows = rows, 
            CTRL = CTRL)

perm_Mantel(com = jac_df$jac_tur,
            env = jac_df$lighthet,
            iter = 999, 
            rows = rows, 
            CTRL = CTRL)

perm_Mantel(com = jac_df$jac_nes,
            env = jac_df$lighthet,
            iter = 999, 
            rows = rows, 
            CTRL = CTRL)

# Modeling Distance-Dissimilarity

newdf <- expand.grid(lighthet = sqrt(seq(0,1, by = .01)),
                     site_code = unique(jac_df$site_code))

mod_re_jtot <- lmer(jac_dis ~ sqrt(lighthet) + (1 + sqrt(lighthet) | site_code), data = jac_df)
mod_re_jtur <- lmer(jac_tur ~ sqrt(lighthet) + (1 + sqrt(lighthet) | site_code), data = jac_df)
mod_re_jnes <- lmer(jac_nes ~ sqrt(lighthet) + (1 + sqrt(lighthet) | site_code), data = jac_df)

pred_dat <- newdf %>%
  bind_cols(list(jac_tot = predict(mod_re_jtot, newdf),
                 jac_tur = predict(mod_re_jtur, newdf),
                 jac_nes = predict(mod_re_jnes, newdf))) %>%
  gather(key = "source",
         value = "dis", 
         -site_code, - lighthet)

pred_dat %>%
  ggplot(aes(x = lighthet ^ 2,
             y = dis,
             color = site_code)) +
  geom_line() + 
  guides(color = FALSE) +
  facet_wrap(~source, nrow = 1)

tot_pred <- emmeans(mod_re_jtot, var = "~lighthet", specs = "lighthet", at = list(lighthet = seq(0,1, by = .01))) %>%
  as.data.frame(.) %>%
  rename("est_tot" = "emmean",
         "tot_lower" = "asymp.LCL" ,
         "tot_higher" = "asymp.UCL") %>%
  select(-SE, -df)

tur_pred <- emmeans(mod_re_jtur, var = "~lighthet", specs = "lighthet", at = list(lighthet = seq(0,1, by = .01))) %>%
  as.data.frame(.) %>%
  rename("est_tur" = "emmean",
         "tur_lower" = "asymp.LCL" ,
         "tur_higher" = "asymp.UCL") %>%
  select(-SE, -df)

nes_pred <- emmeans(mod_re_jnes, var = "~lighthet", specs = "lighthet", at = list(lighthet = seq(0,1, by = .01)))%>%
  as.data.frame(.) %>%
  rename("est_nes" = "emmean",
         "nes_lower" = "asymp.LCL" ,
         "nes_higher" = "asymp.UCL") %>%
  select(-SE, -df)

spat_jac_plot <- left_join(tot_pred, tur_pred, by = "lighthet") %>%
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
  ggtitle("Spatial Jaccard Dissimilarity") +
  labs(linetype = "Jaccard Dissimilarity") +
  theme_bw()
