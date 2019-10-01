
library(lavaan)
library(lavaan.survey)
library(ggplot2) # Figure plotting
library(tidyverse) # Data manipulation
library(lme4) # Linear mixed effects models
library(vegan) # General ecology functions
library(data.table) # Fast reading of data tables
library(testthat) # Unit testing
library(dtplyr) # Data.table dplyr functions
library(JostDiv) # Exponentiated diversity indices (Hill Numbers)
library(gridExtra) # Multiple GGplots
library(corrplot) # Correlation plots
library(sjPlot) # Model coefficient plotting
library(tidyverse)

# Explanatory variables
explanatory_var <- read.csv("Data/divturnover_explanatory.csv",
                            stringsAsFactors = FALSE,
                            na.strings = c('NA','NULL'))

# Temporal diversity
temporal_div <- read.csv("Data/diversity_temporal_full.csv",
                         stringsAsFactors = FALSE,
                         na.strings = c('NA','NULL'))

# diversity
spatial_div <- read.csv("Data/diversity_spatial_full.csv",
                        stringsAsFactors = FALSE,
                        na.strings = c('NA','NULL'))

# General plot descriptions
plot.descriptions <- read.csv("Data/comb-by-plot-clim-soil-diversity-02-Aug-2019.csv",
                              stringsAsFactors = FALSE,
                              na.strings = c('NA','NULL'))

temporal_table <- temporal_div %>% 
  select(site_code, block, alpha_q0, beta_q0, gamma_q0) %>% 
  group_by(site_code, block) %>%
  summarise(beta_q0 = mean(beta_q0),
            alpha_q0 = mean(alpha_q0),
            gamma_q0 = mean(gamma_q0)) %>%
  rename("temporal_beta_q0" = "beta_q0",
         "temporal_alpha_q0" = "alpha_q0",
         "temporal_gamma_q0" = "gamma_q0")

spatial_table <- spatial_div %>% 
  select(site_code, block, alpha_q0, beta_q0, gamma_q0) %>% 
  group_by(site_code, block) %>%
  summarise(beta_q0 = mean(beta_q0),
            alpha_q0 = mean(alpha_q0),
            gamma_q0 = mean(gamma_q0)) %>%
  rename("spatial_beta_q0" = "beta_q0",
         "spatial_alpha_q0" = "alpha_q0",
         "spatial_gamma_q0" = "gamma_q0")

biomass_summary_spatial <- plot.descriptions %>% 
  filter(year_trt == 0) %>% 
  select(site_code, block, total_mass, proportion_par, Ambient_PAR, live_mass) %>% 
  group_by(site_code, block) %>%
  mutate(proportion_par = proportion_par * 100) %>%
  summarise(spatial_mean_biomass = mean(na.omit(total_mass)),
            spatial_cv_biomass = sd(na.omit(total_mass)) /  mean(na.omit(total_mass)),
            spatial_mean_prod = mean(na.omit(live_mass)),
            spatial_cv_prod = sd(na.omit(live_mass)) / mean(na.omit(live_mass)),
            spatial_mean_par = mean(100 - proportion_par),
            spatial_sd_par = sd(100 - proportion_par),
            spatial_ambient_par = mean(Ambient_PAR))

biomass_summary_temporal <- plot.descriptions %>% 
  filter(year_trt < 5 & trt == "Control" & !is.na(first_nutrient_year)) %>%
  select(site_code, block, total_mass, proportion_par, Ambient_PAR, live_mass) %>%
  group_by(site_code, block) %>%
  mutate(proportion_par = proportion_par * 100) %>%
  summarise(temporal_mean_biomass = mean(na.omit(total_mass)),
            temporal_cv_biomass = sd(na.omit(total_mass)) / mean(na.omit(total_mass)),
            temporal_mean_prod = mean(na.omit(live_mass)),
            temporal_cv_prod = sd(na.omit(live_mass)) / mean(na.omit(live_mass)),
            temporal_mean_par = mean(100 - proportion_par),
            temporal_sd_par = sd(100 - proportion_par),
            temporal_ambient_par = mean(Ambient_PAR)) %>%
  # group_by(site_code, block) %>%
  # summarise(temporal_mean_biomass = mean(na.omit(temporal_mean_biomass)),
  #           temporal_sd_biomass = mean(na.omit(temporal_sd_biomass)),
  #           temporal_mean_par = mean(na.omit(temporal_mean_par)),
  #           temporal_sd_par = sd(na.omit(temporal_sd_par)),
  #           temporal_ambient_par = mean(na.omit(temporal_ambient_par))) %>%
  mutate(temporal_cv_biomass = ifelse(site_code == "ethass.au", NA, temporal_cv_biomass))


explanatory_cvs <- explanatory_var %>%
  transmute(site_code = site_code,
            block = block, 
            cv_N = sd_N / mean_N,
            cv_C = sd_C / mean_C,
            cv_P = sd_P / mean_P,
            cv_K = sd_K / mean_K,
            cv_ppt = sd_ppt / mean_ppt,
            cv_pet = sd_pet / mean_pet)

SEM_table <- spatial_table %>% 
  full_join(temporal_table, by = c("site_code", "block")) %>%
  left_join(biomass_summary_spatial, by = c("site_code", "block")) %>% 
  left_join(biomass_summary_temporal, by = c("site_code", "block")) %>%
  left_join(explanatory_var, by = c("site_code", "block")) %>%
  left_join(explanatory_cvs, by = c("site_code", "block")) %>%
  filter(block < 4)

write.csv(x = SEM_table,
          "Data/SEM_table.csv", 
          row.names = FALSE)

# distr_plot <- SEM_table %>%
#   select(-anthropogenic, - burned, -grazed, -managed) %>%
#   gather(key = "sd", value = "val", -site_code) %>%
#   ggplot(aes(x = val)) +
#   geom_histogram(bins = 30, fill = "white", color = "black") +
#   facet_wrap(~sd, scales = "free")
# 
# ggsave(plot = distr_plot, 
#        "../../Figures/check_distributions.pdf", height = 15, width = 25)


# Spatial dataset SEM

spatial_SEM_dat = SEM_table %>% select(site_code, block,
                                       spatial_beta_q0, spatial_gamma_q0, site_richness,
                                       spatial_cv_biomass, spatial_cv_prod, spatial_sd_par, spatial_mean_biomass,
                                       managed, burned, grazed, 
                                       cv_N, cv_C, cv_K, cv_P, sd_pH) %>%
  na.omit(.) %>%
  filter(cv_P != 0) %>%
  ungroup() %>%
  mutate(managed = factor(managed, levels = c(0,1)),
         burned = factor(burned, levels = c(0,1)),
         grazed = factor(grazed, levels = c(0,1)),
         block = factor(block, levels = c(1,2,3,4,5,6))) %>%
  mutate_if(is.numeric, scale)


soil_composite_model <- lm(spatial_cv_prod ~ cv_N + cv_C + cv_K + cv_P + sd_pH, data = spatial_SEM_dat)
summary(soil_composite_model)

soil_composite_coefs <- coef(soil_composite_model)[2:6]
soil_composite <- soil_composite_coefs[1] * spatial_SEM_dat$cv_N +
  soil_composite_coefs[2] * spatial_SEM_dat$cv_C + 
  soil_composite_coefs[3] * spatial_SEM_dat$cv_K +
  soil_composite_coefs[4] * spatial_SEM_dat$cv_P +
  soil_composite_coefs[5] * spatial_SEM_dat$sd_pH

plot(spatial_SEM_dat$spatial_cv_prod ~ soil_composite)

biomass_composite_model <- lm(spatial_cv_biomass ~ spatial_cv_prod + managed + grazed + burned, data = spatial_SEM_dat)
summary(biomass_composite_model)
biomass_composite_coefs <- coef(biomass_composite_model)[3:5]
biomass_composite <- biomass_composite_coefs[1] * (as.integer(spatial_SEM_dat$managed) - 1) +
  biomass_composite_coefs[2] * (as.integer(spatial_SEM_dat$grazed) - 1) +
  biomass_composite_coefs[3] * (as.integer(spatial_SEM_dat$burned) - 1)

spatial_SEM_dat = bind_cols(spatial_SEM_dat, soil_composite = soil_composite, biomass_composite = biomass_composite)

sitemod <- '# Composite variables
              spatial_cv_prod ~ soil_composite
              spatial_cv_biomass ~ spatial_cv_prod + biomass_composite
              spatial_beta_q0 ~ spatial_cv_biomass + site_richness

            # Correlated error
              spatial_cv_prod ~~ spatial_cv_biomass
'

sitemod.fit <- sem(sitemod, data=spatial_SEM_dat, fixed.x=FALSE, estimator = "MLM")  # run model
fit_Ctrl <- sem(sitemod, data = spatial_SEM_dat, estimator="MLM")
exp.design_Ctrl <- svydesign(ids=~site_code, nest=TRUE, data=spatial_SEM_dat)
survey.fit.plot_Ctrl <- lavaan.survey(lavaan.fit=fit_Ctrl, survey.design=exp.design_Ctrl)
summary(survey.fit.plot_Ctrl, fit.measures=T, standardized=T, rsquare=T) # This seems to work!

# Repeating the same for temporal data

temporal_SEM_dat = SEM_table %>% select(site_code, block, 
                                       temporal_beta_q0, temporal_gamma_q0, temporal_alpha_q0, site_richness,
                                       temporal_cv_biomass, temporal_cv_prod,
                                       managed, burned, grazed, ann_frac, anthropogenic,
                                       cv_ppt, cv_pet, mean_wtrdef, sd_ppt) %>%
  na.omit(.) %>%
  ungroup() %>%
  mutate(managed = factor(managed, levels = c(0,1)),
         burned = factor(burned, levels = c(0,1)),
         grazed = factor(grazed, levels = c(0,1)),
         block = factor(block, levels = c(1,2,3,4,5,6))) %>%
  mutate_if(is.numeric, scale) %>%
  filter(site_code != "marc.ar")


clim_composite_model <- lm(temporal_cv_prod ~ sd_ppt*mean_wtrdef, data = temporal_SEM_dat)
summary(clim_composite_model)

clim_composite_coefs <- coef(clim_composite_model)[2:4]
clim_composite <- clim_composite_coefs[1] * temporal_SEM_dat$cv_ppt +
                  clim_composite_coefs[2] * temporal_SEM_dat$mean_wtrdef + 
                  clim_composite_coefs[3] * (temporal_SEM_dat$cv_ppt *  temporal_SEM_dat$mean_wtrdef) 

biomass_composite_model <- lm(temporal_cv_biomass ~ temporal_cv_prod + managed + grazed + burned, data = temporal_SEM_dat)
summary(biomass_composite_model)

biomass_composite_coefs <- coef(biomass_composite_model)[3:5]
biomass_composite <- biomass_composite_coefs[1] * (as.integer(temporal_SEM_dat$managed) - 1) +
  biomass_composite_coefs[2] * (as.integer(temporal_SEM_dat$grazed) - 1) +
  biomass_composite_coefs[3] * (as.integer(temporal_SEM_dat$burned) - 1)

temporal_SEM_dat = bind_cols(temporal_SEM_dat, clim_composite = clim_composite, biomass_composite = biomass_composite)

p <- temporal_SEM_dat %>%
  ggplot(aes(x = temporal_cv_biomass,
             y = temporal_beta_q0,
             color = site_code)) +
  geom_point()

ggplotly(p)



summary(lm(temporal_beta_q0 ~ site_richness, data= temporal_SEM_dat))

summary(lm(temporal_beta_q0 ~  temporal_cv_biomass, data= temporal_SEM_dat))




sitemod <- '# Composite variables
              temporal_cv_prod ~ clim_composite
              temporal_cv_biomass ~ temporal_cv_prod + biomass_composite
              temporal_beta_q0 ~ temporal_cv_biomass + temporal_gamma_q0

            # Correlated error
              temporal_cv_prod ~~ temporal_cv_biomass'


sitemod.fit <- sem(sitemod, data=temporal_SEM_dat, fixed.x=FALSE, estimator = "MLM")  # run model
summary(sitemod.fit, fit.measures = TRUE, standardized = TRUE)

fit_Ctrl <- sem(sitemod, data = temporal_SEM_dat, estimator="MLM")
exp.design_Ctrl <- svydesign(ids=~site_code, nest=TRUE, data=temporal_SEM_dat)
survey.fit.plot_Ctrl <- lavaan.survey(lavaan.fit=fit_Ctrl, survey.design=exp.design_Ctrl)
summary(survey.fit.plot_Ctrl, fit.measures=T, standardized=T, rsquare=T) # This seems to work!



# Combining these together?

temporal_SEM_dat <- temporal_SEM_dat %>% rename("temp_biomass_composite" = "biomass_composite") %>% select(-site_richness)
spatial_SEM_dat <- spatial_SEM_dat %>% rename("spat_biomass_composite" = "biomass_composite")%>% select(-site_richness)

all_SEM_dat <- full_join(temporal_SEM_dat, spatial_SEM_dat)

p2 <- all_SEM_dat %>%
  ggplot(aes(x= spatial_cv_biomass,
      y = temporal_cv_biomass,
      color = site_code)) +
  geom_point()
ggplotly(p2)


summary(lm(temporal_beta_q0 ~ temporal_cv_biomass + spatial_cv_biomass, data= all_SEM_dat))


sitemod <- '# Composite variables
              temporal_cv_prod ~ clim_composite
              temporal_cv_biomass ~ temporal_cv_prod + temp_biomass_composite
              
              spatial_cv_prod ~ soil_composite
              spatial_cv_biomass ~ spatial_cv_prod + spat_biomass_composite


            # Standard regressions
              temporal_beta_q0 ~ temporal_cv_biomass + temporal_gamma_q0
              spatial_beta_q0 ~ spatial_cv_biomass  + spatial_gamma_q0

            # Correlated error
              spatial_cv_prod ~~ spatial_cv_biomass
              temporal_cv_prod ~~ temporal_cv_biomass
              spatial_beta_q0 ~~ temporal_beta_q0
'

bigmod.fit <- sem(sitemod, data=na.omit(all_SEM_dat), fixed.x=FALSE, estimator = "MLM")  # run model
summary(bigmod.fit, fit.measures = TRUE, rsquare = T)

View(na.omit(all_SEM_dat))
exp.design_Ctrl <- svydesign(ids=~site_code, nest=TRUE, data=na.omit(all_SEM_dat))
survey.fit.plot_Ctrl <- lavaan.survey(lavaan.fit=bigmod.fit, survey.design=exp.design_Ctrl, estimator = "MLMV")
summary(survey.fit.plot_Ctrl, fit.measures=T, standardized=T, rsquare=T) # This seems to work!

