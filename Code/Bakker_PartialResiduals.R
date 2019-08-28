standard.ggsave <- function(filename = filename) {
  ggsave(filename, height = 6, width = 10, units = "in", dpi = 600)
}

# Alternative plot of partial residuals
bray_dat <- site.summary %>% dplyr::select(bray, log.S, total_mass, total_mass_CV, ANN_TEMP_RANGE, MAP)

bray_runs <- rbind(data.frame(resid = (resid(lm(bray ~ . - log.S, bray_dat))), resp = "bray", var = bray_dat$log.S, varname = "log.S"),
                   data.frame(resid = (resid(lm(bray ~ . - total_mass, bray_dat))), resp = "bray", var = bray_dat$total_mass, varname = "total_mass"),
                   data.frame(resid = (resid(lm(bray ~ . - total_mass_CV, bray_dat))), resp = "bray", var = bray_dat$total_mass_CV, varname = "total_mass_cv"),
                   data.frame(resid = (resid(lm(bray ~ . - ANN_TEMP_RANGE, bray_dat))), resp = "bray", var = bray_dat$ANN_TEMP_RANGE, varname = "ann_temp_range"),
                   data.frame(resid = (resid(lm(bray ~ . - MAP, bray_dat))), resp = "bray", var = bray_dat$MAP, varname = "MAP"))

sor_dat <- site.summary %>% dplyr::select(sor, log.S, total_mass, total_mass_CV, ANN_TEMP_RANGE, MAP)

sor_runs <- rbind(data.frame(resid = (resid(lm(sor ~ . - log.S, sor_dat))), resp = "sor", var = sor_dat$log.S, varname = "log.S"),
                  data.frame(resid = (resid(lm(sor ~ . - total_mass, sor_dat))), resp = "sor", var = sor_dat$total_mass, varname = "total_mass"),
                  data.frame(resid = (resid(lm(sor ~ . - total_mass_CV, sor_dat))), resp = "sor", var = sor_dat$total_mass_CV, varname = "total_mass_cv"),
                  data.frame(resid = (resid(lm(sor ~ . - ANN_TEMP_RANGE, sor_dat))), resp = "sor", var = sor_dat$ANN_TEMP_RANGE, varname = "ann_temp_range"),
                  data.frame(resid = (resid(lm(sor ~ . - MAP, sor_dat))), resp = "sor", var = sor_dat$MAP, varname = "MAP"))

prop_bal_dat <- site.summary %>% dplyr::select(prop_bal, log.S, total_mass, total_mass_CV, ANN_TEMP_RANGE, MAP)

prop_bal_runs <- rbind(data.frame(resid = (resid(lm(prop_bal ~ . - log.S, prop_bal_dat))), resp = "prop_bal", var = prop_bal_dat$log.S, varname = "log.S"),
                       data.frame(resid = (resid(lm(prop_bal ~ . - total_mass, prop_bal_dat))), resp = "prop_bal", var = prop_bal_dat$total_mass, varname = "total_mass"),
                       data.frame(resid = (resid(lm(prop_bal ~ . - total_mass_CV, prop_bal_dat))), resp = "prop_bal", var = prop_bal_dat$total_mass_CV, varname = "total_mass_cv"),
                       data.frame(resid = (resid(lm(prop_bal ~ . - ANN_TEMP_RANGE, prop_bal_dat))), resp = "prop_bal", var = prop_bal_dat$ANN_TEMP_RANGE, varname = "ann_temp_range"),
                       data.frame(resid = (resid(lm(prop_bal ~ . - MAP, prop_bal_dat))), resp = "prop_bal", var = prop_bal_dat$MAP, varname = "MAP"))

prop_sim_dat <- site.summary %>% dplyr::select(prop_sim, log.S, total_mass, total_mass_CV, ANN_TEMP_RANGE, MAP)

prop_sim_runs <- rbind(data.frame(resid = (resid(lm(prop_sim ~ . - log.S, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$log.S, varname = "log.S"),
                       data.frame(resid = (resid(lm(prop_sim ~ . - total_mass, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$total_mass, varname = "total_mass"),
                       data.frame(resid = (resid(lm(prop_sim ~ . - total_mass_CV, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$total_mass_CV, varname = "total_mass_cv"),
                       data.frame(resid = (resid(lm(prop_sim ~ . - ANN_TEMP_RANGE, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$ANN_TEMP_RANGE, varname = "ann_temp_range"),
                       data.frame(resid = (resid(lm(prop_sim ~ . - MAP, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$MAP, varname = "MAP"))


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

standard.ggsave(filename = "../Figures/Bakker_PartialCor_All.png")


# Alternative plot of partial residuals
bray_dat <- site.summary %>% dplyr::select(bray, total_mass)

bray_runs <- rbind(data.frame(resid = (resid(lm(bray ~ . - total_mass, bray_dat))), resp = "bray", var = bray_dat$total_mass, varname = "total_mass"))

sor_dat <- site.summary %>% dplyr::select(sor, log.S,total_mass_CV, ANN_TEMP_RANGE)

sor_runs <- rbind(data.frame(resid = (resid(lm(sor ~ . - log.S, sor_dat))), resp = "sor", var = sor_dat$log.S, varname = "log.S"),
                  data.frame(resid = (resid(lm(sor ~ . - total_mass_CV, sor_dat))), resp = "sor", var = sor_dat$total_mass_CV, varname = "total_mass_cv"),
                  data.frame(resid = (resid(lm(sor ~ . - ANN_TEMP_RANGE, sor_dat))), resp = "sor", var = sor_dat$ANN_TEMP_RANGE, varname = "ann_temp_range"))

prop_bal_dat <- site.summary %>% dplyr::select(prop_bal, log.S, MAP)

prop_bal_runs <- rbind(data.frame(resid = (resid(lm(prop_bal ~ . - log.S, prop_bal_dat))), resp = "prop_bal", var = prop_bal_dat$log.S, varname = "log.S"),
                       data.frame(resid = (resid(lm(prop_bal ~ . - MAP, prop_bal_dat))), resp = "prop_bal", var = prop_bal_dat$MAP, varname = "MAP"))

prop_sim_dat <- site.summary %>% dplyr::select(prop_sim, log.S, total_mass, ANN_TEMP_RANGE)

prop_sim_runs <- rbind(data.frame(resid = (resid(lm(prop_sim ~ . - log.S, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$log.S, varname = "log.S"),
                       data.frame(resid = (resid(lm(prop_sim ~ . - total_mass, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$total_mass, varname = "total_mass"),
                       data.frame(resid = (resid(lm(prop_sim ~ . - ANN_TEMP_RANGE, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$ANN_TEMP_RANGE, varname = "ann_temp_range"))

sigvals <- master %>% dplyr::select(Response, Explan, "Sig.AIC") %>% distinct() %>%
  mutate(Explan = case_when(Explan == "ANN_TEMP_RANGE" ~ "ann_temp_range",
                            Explan == "total_mass_CV" ~ "total_mass_cv",
                            TRUE ~ as.character(Explan))) %>%
  dplyr::rename("resp" = "Response",
         "varname" = "Explan",
         "Sig" = "Sig.AIC")

bind_rows(bray_runs, sor_runs, prop_bal_runs, prop_sim_runs) %>%
  right_join(sigvals) %>% 
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
  ggtitle("Partial Correlation Plots -- significant variables")

standard.ggsave(filename = "../Figures/Bakker_PartialCor_Sigonly.png")

# Alternative plot of partial residuals
bray_dat <- site.summary %>% dplyr::select(bray, total_mass)

bray_runs <- rbind(data.frame(resid = (resid(lm(bray ~ . - total_mass, bray_dat))), resp = "bray", var = bray_dat$total_mass, varname = "total_mass", order = 1))

sor_dat <- site.summary %>% dplyr::select(sor, log.S,total_mass_CV, ANN_TEMP_RANGE)

sor_runs <- rbind(data.frame(resid = (resid(lm(sor ~ 1, sor_dat))), resp = "sor", var = sor_dat$log.S, varname = "log.S", order = 1),
                  data.frame(resid = (resid(lm(sor ~ log.S, sor_dat))), resp = "sor", var = sor_dat$total_mass_CV, varname = "total_mass_cv", order = 2),
                  data.frame(resid = (resid(lm(sor ~ log.S + total_mass_CV, sor_dat))), resp = "sor", var = sor_dat$ANN_TEMP_RANGE, varname = "ann_temp_range", order = 3))

prop_bal_dat <- site.summary %>% dplyr::select(prop_bal, log.S, MAP)

prop_bal_runs <- rbind(data.frame(resid = (resid(lm(prop_bal ~ 1, prop_bal_dat))), resp = "prop_bal", var = prop_bal_dat$log.S, varname = "log.S", order = 1),
                       data.frame(resid = (resid(lm(prop_bal ~ 1 + log.S, prop_bal_dat))), resp = "prop_bal", var = prop_bal_dat$MAP, varname = "MAP", order = 2))

prop_sim_dat <- site.summary %>% dplyr::select(prop_sim, log.S, total_mass, ANN_TEMP_RANGE)

prop_sim_runs <- rbind(data.frame(resid = (resid(lm(prop_sim ~ 1, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$log.S, varname = "log.S", order = 1),
                       data.frame(resid = (resid(lm(prop_sim ~ 1 + log.S, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$total_mass, varname = "total_mass", order = 2),
                       data.frame(resid = (resid(lm(prop_sim ~ 1 + log.S + total_mass, prop_sim_dat))), resp = "prop_sim", var = prop_sim_dat$ANN_TEMP_RANGE, varname = "ann_temp_range", order = 3))

sigvals <- master %>% dplyr::select(Response, Explan, "Sig.AIC") %>% distinct() %>%
  mutate(Explan = case_when(Explan == "ANN_TEMP_RANGE" ~ "ann_temp_range",
                            Explan == "total_mass_CV" ~ "total_mass_cv",
                            TRUE ~ as.character(Explan))) %>%
  dplyr::rename("resp" = "Response",
                "varname" = "Explan",
                "Sig" = "Sig.AIC")

pc_dat <- bind_rows(bray_runs, sor_runs, prop_bal_runs, prop_sim_runs) %>%
  right_join(sigvals) %>% 
  mutate(resp = case_when(resp == "bray" ~ "Bray-Curtis Dissimilarity",
                          resp == "sor" ~ "Sorensen Dissimilarity",
                          resp == "prop_bal" ~ "Balanced Variation (%)",
                          resp == "prop_sim" ~ "Species Turnover (%)")) %>%
  replace_na(list(Sig = "Sig")) %>% 
  group_by(resp) %>%
  mutate(resid = scale(resid)) %>%
  mutate(varname = case_when(varname == "log.S" ~ "Log(S)",
                         varname == "MAP" ~ "Mean Ann Precip",
                         varname == "total_mass" ~ "Total Biomass",
                         varname == "total_mass_cv" ~ "CV Biomass",
                         varname == "ann_temp_range" ~ "Ann Temp Range")) %>%
  mutate(varname = factor(varname,
                       levels = c("Log(S)",
                                  "Mean Ann Precip",
                                  "Total Biomass",
                                  "CV Biomass",
                                  "Ann Temp Range"))) %>%
  mutate(resp = factor(resp,
                       levels = c("Bray-Curtis Dissimilarity",
                                  "Sorensen Dissimilarity",
                                  "Balanced Variation (%)",
                                  "Species Turnover (%)")))  %>% 
  mutate(resid = case_when(is.na(resid) ~ 0,
                           TRUE ~ as.numeric(resid))) %>%
  mutate(siglabel = case_when(Sig == "NS" ~ "N.S.",
                              Sig == "Sig" ~ ""))

labdat = pc_dat %>%   
  ungroup() %>%
  dplyr::group_by(varname) %>%
  dplyr::summarise(meanval = mean(range(na.omit(var))))

left_join(pc_dat, labdat) %>%
  ggplot() +
  geom_rect(aes(fill = as.factor(Sig)), 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
  geom_point(aes(x = as.numeric(var),
                 y = as.numeric(resid))) +
  # geom_text(aes(label = as.character(order)), x = -Inf, y = Inf, hjust = -0.4, vjust = 1.1) +
  geom_text(aes(label = siglabel,
                x = meanval,
                y = 0)) +
  facet_grid(resp ~ varname, scales = "free_x") +
  stat_smooth(method = "lm", se = FALSE, color = "black", size = 1,
              aes(x = as.numeric(var),
                  y = as.numeric(resid),
                  linetype = Sig)) + 
  theme_bw(base_size = 10) +
  ylim(-3, 3) +
  scale_fill_manual(values = c("white", "grey")) +
  labs(x = "Explanatory Variable", y = "Scaled Residuals") + 
  guides(linetype = FALSE, fill = FALSE) +
  ggtitle("Partial Correlation of Explanatory Variables")

standard.ggsave(filename = "../Figures/Bakker_PartialCor_Sequential_ordered.png")

