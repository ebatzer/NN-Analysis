# Nutrient Co-Limitation and Grassland Variance, 2019.
#
# (Press 'alt-o' (windows) / 'alt-command-o' (mac) to collapse script into neat list)
#
# ----
# 1.1 BEGIN: Load Packages ----

#install.packages("")

library(arm)
library(tidyverse)
library(readxl)
library(gridExtra) # !masks dplyr 'combine'!
library(lme4) # !masks tidyr 'expand'
library(nlme) # !masks lmList, collapse
library(emmeans)
library(broom)
library(visreg)
library(car)
library(vegan)
library(sjPlot)
library(openxlsx)
library(RColorBrewer)
library(lmerTest)

options(scipen=999)

# 1.2 Load Data ----

#     Upload is currently 25th January 2019 version, available NutNet.
Nut <- read_excel("comb-by-plot-clim-soil-diversity-25-Jan-2019.xlsx", na = "NA")

# & Species cover:
#     Upload is currently 25th January 2019 version, available NutNet.
Cover <- read_excel("full-cover-25-January-2019.xlsx", na = "NA")

head(Cover)

# ----

# Objects ----

YearFilter <- Nut %>% 
  group_by(site_code) %>% 
  summarise(MaxTrtYr = max(year_trt)) 

Nut2 <- Nut %>% 
  select(site_code, continent, habitat, block, plot, region, year_trt, 
         trt, live_mass, dead_mass, richness_vegan, MAP_v2, pct_N) %>%
  merge(YearFilter, all.x=T) %>%
  filter(MaxTrtYr >=7 ) %>%
  filter(year_trt <= 7) %>%
  filter(trt!='NPK+Fence', trt!='Fence')

rm(YearFilter)


Nut3 <- Nut2 %>%
  mutate(Outliers = ifelse(site_code == 'elliot.us' & block == 2 & year_trt == 6 & trt == 'N', "Outlier",
                           ifelse(site_code == 'elliot.us' & block == 2 & year_trt == 6 & trt == 'K', "Outlier",
                                  ifelse(site_code == 'mcla.us' & block == 2 & year_trt == 1 & trt == 'K', "Outlier",
                                         ifelse(site_code == 'hall.us' & block == 1 & year_trt == 3 & trt == 'NK', "Outlier",
                                                ifelse(site_code == 'elliot.us' & block == 1 & year_trt == 6 & trt == 'PK', "Outlier",
                                                       "Retain")))))) %>%
  filter(Outliers != 'Outlier') %>%
  filter(site_code != 'sevi.us') %>%
  select(site_code, block, plot, region, year_trt, trt, live_mass, dead_mass, richness_vegan, MAP_v2, pct_N)

# ----
# ANPP ~ trt*year ----
# LRR ----
# Log-Log ----

NutLog <- Nut2 %>% mutate(log_live = log(live_mass + 1)) %>% 
  mutate(log_year = log(year_trt + 1))

logmod_quad <- lmer(log_live ~ log_year + trt + log_year:trt + I(log_year^2) + I(log_year^2):trt
              +(1|site_code/block), data=NutLog, REML = FALSE, na.action = na.exclude) 

logmod_lin <- lmer(log_live ~ log_year + trt + log_year:trt +
                   (1|site_code/block), 
                   data = NutLog, REML = FALSE, na.action = na.exclude) 

anova(logmod_quad, logmod_lin)

plot(logmod_quad)
plot(logmod_noint)

summary(logmod_quad)
tab_model(logmod_quad)
anova(logmod_quad)
qqnorm(Mod7)
hist(residuals(logmod_quad))
plot(fitted(logmod_quad),residuals(logmod_quad))

summary(logmod_lin)
tab_model(logmod_lin)
anova(logmod_lin)
qqnorm(logmod_lin)
hist(residuals(logmod_lin))
plot(fitted(logmod_lin),residuals(logmod_lin))

Exp_Figure_quad <- data.frame(year_trt = sort(rep(c(0:7), length(unique(NutLog$trt)))),
                         trt = rep(unique(NutLog$trt), 8)) %>%
  mutate(log_year = log(year_trt + 1)) %>%
  mutate(predval = predict(logmod_quad, newdata = ., re.form = NA)) %>%
  mutate(exp_year_trt = (exp(log_year) -1)) %>%
  mutate(exp_live = (exp(predval)-1))

ggplot(data = Exp_Figure_quad, 
       aes(x = exp_year_trt,
           y = exp_live,
           color = trt)) +
  geom_point() +
  geom_smooth(formula = y ~ x + x^2, se = FALSE)

Exp_Figure_lin <- data.frame(year_trt = sort(rep(c(0:7), length(unique(NutLog$trt)))),
                              trt = rep(unique(NutLog$trt), 8)) %>%
  mutate(log_year = log(year_trt + 1)) %>%
  mutate(predval = predict(logmod_lin, newdata = ., re.form = NA)) %>%
  mutate(exp_year_trt = (exp(log_year) -1)) %>%
  mutate(exp_live = (exp(predval)-1))

ggplot(data = Exp_Figure_lin, 
       aes(x = exp_year_trt,
           y = exp_live,
           color = trt)) +
  geom_point() +
  geom_smooth(formula = y ~ x, se = FALSE)

calcresid2 <- cbind(NutLog, predval = predict(logmod_lin, type = "response")) %>%
  group_by(trt, year_trt) %>%
  mutate(SSres = sqrt((exp(predval)-1) - live_mass) ^2) %>%
  summarise(mean_sd = mean(na.omit(SSres))) 

ggplot(data = calcresid2,
       aes(x = year_trt,
           y = mean_sd,
           color = trt)) +
  geom_point() +
  stat_smooth(formula = y ~ x + x^2, se = FALSE)

# ---- 
# 7 year Mean ----

yr7M <- Nut3 %>%
  select(site_code,block,plot,year_trt,trt,live_mass) %>%
  filter(year_trt != 0) %>%
  group_by(site_code, block, trt) %>%
  summarise(MEAN = mean(live_mass, na.rm = T)) %>%
  spread(trt, MEAN, fill=NA, convert = FALSE, drop = TRUE) %>%
  mutate(N=log(N/Control),P=log(P/Control),
         K=log(K/Control),NP=log(NP/Control),
         NK=log(NK/Control),PK=log(PK/Control),
         NPK=log(NPK/Control)) %>%
  select(-Control) %>%
  gather(trt, LRRMean, -site_code, -block)

ggplot(data = (yr7M %>%
                 group_by(trt) %>%
                 summarise(MEAN = mean(LRRMean))), 
       aes(x=trt, y = MEAN)) +
  geom_bar(stat = "identity")

Mod7M_1 <- lmer(LRRMean ~ trt + (1|site_code/block),
                data=yr7M, REML = FALSE, na.action=na.exclude)
Mod7M_2 <- lmer(LRRMean ~ trt + (1|site_code),
                data=yr7M, REML = FALSE, na.action=na.exclude)
Mod7M_3 <- lmer(LRRMean ~ 1 + (1|site_code/block),
                data=yr7M, REML = FALSE, na.action=na.exclude)
Mod7M_4 <- lmer(LRRMean ~ trt + (1|site_code),
                data=yr7M, REML = FALSE, na.action=na.exclude)

anova(Mod7M_1, Mod7M_2, Mod7M_3,Mod7M_4)

rm(Mod7M_2, Mod7M_3,Mod7M_4)

qqnorm(Mod7M_1)
hist(residuals(Mod7M_1))
plot(fitted(Mod7M_1),residuals(Mod7M_1))

summary(Mod7M_1)
tab_model(Mod7M_1)
pairs(emmeans(Mod7M_1, "trt"))
plot(emmeans(Mod7M_1, "trt"), comparisons = TRUE)

# TAKE CIs FROM TAB MOD

# Model matrix of your fitted mixed effects model
mod_matrix <- unique(model.matrix(Mod7M_1))
# Calculate the confidence intervals of fixed and random effects
confints <- matrix(confint(Mod7M_1_lmer, nsim = 9999), ncol = 2)

# Dataframe of treatment levels
Mod7M_1_DF <- data.frame(trt = rep(unique(yr7M$trt),1))
  # Calculate overall predicted values from model
Mod7M_1_DF <- Mod7M_1_DF %>%
  mutate(predval = predict(Mod7M_1, newdata = Mod7M_1_DF, re.form = NA)) 
  # Matrix multiply by confidence interval matrix to generate upper and lower bounds
  mutate(ciLow = mod_matrix %*% confints[4:10,1],
         ciHigh = mod_matrix %*% confints[4:10,2]) 

Mod7M_1_DF <- new_tibble(Mod7M_1_DF)

ggplot(data = (Mod7M_1_DF %>%
                 mutate(trt = factor(trt, levels = 
                                       c("N", "P", "K", "NP", "NK", "PK","NPK"))) %>%
                 mutate(predval = round(predval, digits=2))),
       aes(x=trt, y = predval)) +
  geom_bar(stat = "identity", colour ="black", fill ="grey80") +
  #geom_errorbar(aes(ymin=ciLow, ymax=ciHigh,width=.3), colour = 'black') +
  geom_hline(yintercept = 0, color = 'red', lty = 1, size = 1) +
  geom_text(aes(label=predval), vjust=1.5,size=3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill= NA, colour = "black")) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)) +
  labs(x="Nutrient treatment", y="Change in 7 year mean ANPP (LRR.Mean)") +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))

# 7 year SD ----

yr7SD <- Nut3 %>%
  select(site_code,block,plot,year_trt,trt,live_mass) %>%
  filter(year_trt != 0) %>%
  group_by(site_code, block, trt) %>%
  summarise(SD = sd(live_mass, na.rm = T)) %>%
  spread(trt, SD, fill=NA, convert = FALSE, drop = TRUE) %>%
  mutate(N=log(N/Control),P=log(P/Control),
         K=log(K/Control),NP=log(NP/Control),
         NK=log(NK/Control),PK=log(PK/Control),
         NPK=log(NPK/Control)) %>%
  select(-Control) %>%
  gather(trt, LRRsd, -site_code, -block)

ggplot(data = (yr7SD %>%
                 group_by(trt) %>%
                 summarise(MEAN = mean(LRRsd))), 
       aes(x=trt, y = MEAN)) +
  geom_bar(stat = "identity")

Mod7sd_1 <- lmer(LRRsd ~ trt + (1|site_code/block),
                data=yr7SD, REML = FALSE, na.action=na.exclude)
Mod7sd_2 <- lmer(LRRsd ~ trt + (1|site_code),
                 data=yr7SD, REML = FALSE, na.action=na.exclude)
Mod7sd_3 <- lmer(LRRsd ~ 1 + (1|site_code/block),
                 data=yr7SD, REML = FALSE, na.action=na.exclude)
Mod7sd_4 <- lmer(LRRsd ~ 1 + (1|site_code),
                 data=yr7SD, REML = FALSE, na.action=na.exclude)

anova(Mod7sd_1, Mod7sd_2, Mod7sd_3,Mod7sd_4)

rm(Mod7sd_2, Mod7sd_3,Mod7sd_4)

qqnorm(Mod7sd_1)
hist(residuals(Mod7sd_1))
plot(fitted(Mod7sd_1),residuals(Mod7sd_1))

summary(Mod7sd_1)
tab_model(Mod7sd_1)
pairs(emmeans(Mod7sd_1, "trt"))
plot(emmeans(Mod7sd_1, "trt"), comparisons = TRUE)

Mod7sd_1_DF <- data.frame(trt = rep(unique(yr7SD$trt),1))
Mod7sd_1_DF <- Mod7sd_1_DF %>%
  mutate(predval = predict(Mod7sd_1, newdata = Mod7sd_1_DF, re.form = NA))

ggplot(data = (Mod7sd_1_DF %>%
                 mutate(trt = factor(trt, levels = 
                                       c("N", "P", "K", "NP", "NK", "PK","NPK"))) %>%
                 mutate(predval = round(predval, digits=2))),
       aes(x=trt, y = predval)) +
  geom_bar(stat = "identity", colour ="black", fill ="grey80") +
  geom_text(aes(label=predval), vjust=1.5,size=3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill= NA, colour = "black")) +
  expand_limits(y = c(0,0.8)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)) +
  labs(x="Nutrient treatment", y="Change in 7 year mean SD (LRR.SD)") +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))

# 7 year S ----

yr7S <- yr7M %>%
  left_join(yr7SD, by = c("site_code", "block", "trt")) %>%
  mutate(LRRs = LRRMean-LRRsd)
  

ggplot(data = (yr7S %>%
                 group_by(trt) %>%
                 summarise(MEAN = mean(LRRs))), 
       aes(x=trt, y = MEAN)) +
  geom_bar(stat = "identity")

Mod7s_1 <- lmer(LRRs ~ trt + (1|site_code/block),
                 data=yr7S, REML = FALSE, na.action=na.exclude)
Mod7s_2 <- lmer(LRRs ~ trt + (1|site_code),
                 data=yr7S, REML = FALSE, na.action=na.exclude)
Mod7s_3 <- lmer(LRRs ~ 1 + (1|site_code/block),
                 data=yr7S, REML = FALSE, na.action=na.exclude)
Mod7s_4 <- lmer(LRRs ~ 1 + (1|site_code),
                 data=yr7S, REML = FALSE, na.action=na.exclude)

anova(Mod7s_1, Mod7s_2, Mod7s_3,Mod7s_4)

rm(Mod7s_2, Mod7s_3,Mod7s_4)

qqnorm(Mod7s_1)
hist(residuals(Mod7s_1))
plot(fitted(Mod7s_1),residuals(Mod7s_1))

summary(Mod7s_1)
tab_model(Mod7s_1)
pairs(emmeans(Mod7s_1, "trt"))
plot(emmeans(Mod7s_1, "trt"), comparisons = TRUE)

Mod7s_1_DF <- data.frame(trt = rep(unique(yr7S$trt),1))
Mod7s_1_DF <- Mod7s_1_DF %>%
  mutate(predval = predict(Mod7s_1, newdata = Mod7s_1_DF, re.form = NA)) 

ggplot(data = (Mod7s_1_DF %>%
                 mutate(trt = factor(trt, levels = 
                                       c("N", "P", "K", "NP", "NK", "PK","NPK"))) %>%
                 mutate(predval = round(predval, digits=2))),
       aes(x=trt, y = predval)) +
  geom_bar(stat = "identity", colour ="black", fill ="grey80") +
  geom_hline(yintercept = 0, color = 'red', lty = 1, size = 1) +
  geom_text(aes(label=predval), vjust=-1,size=3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill= NA, colour = "black")) +
  labs(x="Nutrient treatment", y="Change in 7 year mean S (LRR.Mean - LRR.SD)") +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))

Mod7MSD <- Mod7M_1_DF %>%
  mutate(MEAN = predval) %>%
  select(-predval) %>%
  left_join(Mod7sd_1_DF, by = "trt") %>%
  mutate(SD = predval) %>%
  select(-predval) %>%
  mutate(Diff = SD-MEAN) %>%
  select(-SD) %>%
  gather(MSD, Pred, -trt)

ggplot(data = (Mod7MSD %>%
                 mutate(trt = factor(trt, levels = 
                                       c("K", "P", "N", "PK", "NK", "NP","NPK")))),
       aes(x=trt, y = Pred, colour = MSD)) +
  geom_bar(stat = "identity", position = "stack", colour ="black", aes(fill = MSD)) +
  geom_hline(yintercept = 0, color = 'blue', lty = 1, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill= NA, colour = "black")) +
  scale_fill_manual(name="", 
                    values = c("Diff"="red", "MEAN"="grey80"),
                    labels = c("LRR.SD", "LRR.Mean")) +
  labs(x="", y="Change in 7 year ANPP mean (LRR.Mean) and SD (LRR.SD)") +
  expand_limits(y = c(0,0.8)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)) +
  theme(legend.position = c(0.2, 0.88), legend.key=element_blank(),
        legend.text=element_text(size=7,color = "black")) +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))


# Figure2_MSD ----

Figure2_MSD <- list()
Figure2_MSD[[1]] <- ggplot(data = (Mod7M_1_DF %>%
                                     mutate(trt = factor(trt, levels = 
                                                           c("N", "P", "K", "NP", "NK", "PK","NPK"))) %>%
                                     mutate(predval = round(predval, digits=2))),
                           aes(x=trt, y = predval)) +
  geom_bar(stat = "identity", colour ="black", fill ="grey80") +
  geom_hline(yintercept = 0, color = 'red', lty = 1, size = 1) +
  geom_text(aes(label=predval), vjust=1.5,size=3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill= NA, colour = "black")) +
  expand_limits(y = c(0,0.6)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6)) +
  labs(x="", y="Change in 7 year mean ANPP (LRR.Mean)") +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))

Figure2_MSD[[2]] <- ggplot(data = (Mod7sd_1_DF %>%
                                     mutate(trt = factor(trt, levels = 
                                                           c("N", "P", "K", "NP", "NK", "PK","NPK"))) %>%
                                     mutate(predval = round(predval, digits=2))),
                           aes(x=trt, y = predval)) +
  geom_bar(stat = "identity", colour ="black", fill ="grey80") +
  geom_hline(yintercept = 0, color = 'red', lty = 1, size = 1) +
  geom_text(aes(label=predval), vjust=1.5,size=3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill= NA, colour = "black")) +
  expand_limits(y = c(0,0.6)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6)) +
  labs(x="", y="Change in 7 year mean SD (LRR.SD)") +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))

Figure2_MSD[[3]] <- ggplot(data = (Mod7s_1_DF %>%
                                     mutate(trt = factor(trt, levels = 
                                                           c("N", "P", "K", "NP", "NK", "PK","NPK"))) %>%
                                     mutate(predval = round(predval, digits=2))),
                           aes(x=trt, y = predval)) +
  geom_bar(stat = "identity", colour ="black", fill ="grey80") +
  geom_hline(yintercept = 0, color = 'red', lty = 1, size = 1) +
  geom_text(aes(label=predval), vjust=-1,size=3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill= NA, colour = "black")) +
  labs(x="", y="Change in 7 year mean S (LRR.Mean - LRR.SD)") +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))

ggsave("2019_08_Fig_2_MSD.png", marrangeGrob(grobs = Figure2_MSD, nrow=1, ncol=3, top=NULL), 
       width = 45, height = 15, unit = 'cm')

ggsave("2019_08_Fig_2_MSDV.png", marrangeGrob(grobs = Figure2_MSD, nrow=3, ncol=1, top=NULL), 
       width = 15, height = 45, unit = 'cm')

Figure2_MSD_2 <- list()
Figure2_MSD_2[[1]] <-  ggplot(data = (Mod7MSD %>%
                                        mutate(trt = factor(trt, levels = 
                                                              c("K", "P", "N", "PK", "NK", "NP","NPK")))),
                              aes(x=trt, y = Pred, colour = MSD)) +
  geom_bar(stat = "identity", position = "stack", colour ="black", aes(fill = MSD)) +
  geom_hline(yintercept = 0, color = 'blue', lty = 1, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill= NA, colour = "black")) +
  scale_fill_manual(name="", 
                    values = c("Diff"="red", "MEAN"="grey80"),
                    labels = c("LRR.SD", "LRR.Mean")) +
  labs(x="", y="Change in 7 year ANPP mean (LRR.Mean) and SD (LRR.SD)") +
  expand_limits(y = c(0,0.8)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)) +
  theme(legend.position = c(0.2, 0.88), legend.key=element_blank(),
        legend.text=element_text(size=7,color = "black")) +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))

Figure2_MSD_2[[2]] <- ggplot(data = (Mod7s_1_DF %>%
                                       mutate(trt = factor(trt, levels = 
                                                             c("K", "P", "N", "PK", "NK", "NP","NPK")))),
                             aes(x=trt, y = predval)) +
  geom_bar(stat = "identity", colour ="black", fill ="grey80") +
  geom_hline(yintercept = 0, color = 'blue', lty = 1, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill= NA, colour = "black")) +
  labs(x="", y="Change in 7 year ANPP stability (LRR.Mean - LRR.SD)") +
  theme(axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))

ggsave("2019_08_Fig_2_MSD_2.png", marrangeGrob(grobs = Figure2_MSD_2, nrow=1, ncol=2, top=NULL), 
       width = 30, height = 15, unit = 'cm')

# ----
# Functional ----

YearFilter2 <- Cover %>% 
  group_by(site_code) %>% 
  summarise(MaxTrtYr = max(year_trt)) 

FuncGroup <- Cover %>% 
  merge(YearFilter2, all.x=T) %>%
  filter(year_trt <=7) %>%
  filter(MaxTrtYr >=7 ) %>%
  filter(trt!='NPK+Fence', trt!='Fence') %>%
  filter(site_code != 'hogone.us') %>%
  select(year_trt, trt, block, plot, Family, Taxon, functional_group)

rm(YearFilter2)

FuncGroup <- FuncGroup %>%
  group_by(year_trt, trt, block) %>%
  summarise(Func = count(functional_group))




# Installing stan
install.packages("rstan")
library("rstan")

data{
  int<lower=0> N;
  int<lower = 0> K;
  vector[N] y;
  matrix[N, K] x;
}

parameters{
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
}












live_formula <- alist(
  live_mass ~ dnorm(mu, sigma),
  mu <- a + b[trt]*year_trt,
  sigma <- a_sigma + b_sigma[trt]*year_trt,
  
  a ~ dnorm(0, 100),
  b[trt] ~ dnorm(b_bar, var_mean),
  b_bar ~ dnorm(0, 10),
  var_mean ~ dnorm(0, 10),
  
  a_sigma ~ dnorm(0, 100),
  b_sigma[trt] ~ dnorm(sig_bar, var_sig),
  sig_bar ~ dnorm(0, 10),
  var_sig ~ dnorm(0, 10)
)
starting_values <- list(
  a <- mean(Nut2$live_mass),
  b_bar = 0,
  a_sigma = sd(Nut2$live_mass),
  sig_bar = 0,
  var_mean = 1,
  var_sig = 1
)

install.packages("rlang")
devtools::install_github("rmcelreath/rethinking")
library(rethinking)

stan_attempt <- rethinking::map2stan(live_formula,
                         data = Nut2[!is.na(Nut2$live_mass),],
                         iter = 2000,
                         warmup = 1000, 
                         chains = 1,
                         refresh = 500)

