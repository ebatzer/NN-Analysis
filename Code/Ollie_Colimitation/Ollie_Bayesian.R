# Messing with Ollie's Data
nut2 <- read.csv("../Data/Nut2.csv")
model_string <- "model{

# Likelihood
for(i in 1:n){
Y[i]   ~ dnorm(mu[i],sigma)
mu[i] <- beta[1] + beta[trt[i]] * year_trt[i] + sitemean[siteval[i]]
}

# Prior for beta
for(j in 1:C){
  beta[j] ~ dnorm(0,.1)
}

# Prior for site mean
for(k in 1:S){
  sitemean[k] ~ dnorm(0, 100)
}

# Prior for the inverse variance
inv.var   ~ dgamma(0.01, 0.01)
sigma     <- 1/sqrt(inv.var)

}"

dat <- nut2[!is.na(nut2$live_mass),]
Y <- log(dat$live_mass)
year_trt <- log(dat$year_trt + 1)
trt <- as.integer(dat$trt) + 1
n <- nrow(dat)
C <- 1 + length(unique(trt))
siteval <- as.integer(dat$site_code)
S <- length(unique(dat$site_code))

library(rjags)

model <- jags.model(textConnection(model_string), 
                    data = list(Y=Y,n=n,year_trt = year_trt,trt=trt,C=C, siteval = siteval,
                                S = S))

update(model, 10000, progress.bar="none"); # Burnin for 10000 samples

samp <- coda.samples(model, 
                     variable.names=c("beta","sigma", "sitemean"), 
                     n.iter=20000, progress.bar="none")


# Messing with Ollie's Data
nut2 <- read.csv("../Data/Nut2.csv")
model_string <- "model{

# Likelihood
for(i in 1:n){
  Y[i]   ~ dnorm(mu[i], sigma[i])
  mu[i] <- beta[1] + beta[trt[i]] * year_trt[i] + sitemean[siteval[i]]
  sigma[i] <- sig_mean + beta_sig[trt[i]] * year_trt[i]
}

# Prior for beta
for(j in 1:C){
  beta[j] ~ dnorm(0,.1)
}

# Prior for site mean
for(k in 1:S){
  sitemean[k] ~ dnorm(0, 100)
}

# Prior for the sigma change
for(l in 2:C){
  beta_sig[l] ~ dnorm(0, .1)
}

# Prior for the inverse variance
inv.var   ~ dgamma(0.01, 0.01)
sig_mean     <- 1/sqrt(inv.var)

}"

dat <- nut2[!is.na(nut2$live_mass),]
Y <- log(dat$live_mass)
year_trt <- log(dat$year_trt + 1)
trt <- as.integer(dat$trt) + 1
n <- nrow(dat)
C <- 1 + length(unique(trt))
siteval <- as.integer(dat$site_code)
S <- length(unique(dat$site_code))

library(rjags)

model <- jags.model(textConnection(model_string), 
                    data = list(Y=Y,n=n,year_trt = year_trt,trt=trt,C=C, siteval = siteval,
                                S = S))

update(model, 1000); # Burnin for 10000 samples

samp <- coda.samples(model, 
                     variable.names=c("beta","sig_mean", "beta_sig"), 
                     n.iter=2000)

bayescoefs <- summary(samp)
str(bayescoefs)

library(dplyr)

mod_df <- data.frame(Y, trt = as.factor(dat$trt), year_trt) %>% distinct()
modmat <- model.matrix(lm(Y ~ trt : year_trt), data = mod_df)
modmat <- as.matrix(data.frame(modmat))

nrow(mod_df)
nrow(modmat)


# Predicted mean values
predvals <- modmat %*% matrix(bayescoefs$quantiles[1:9,3], nrow = 9)
mod_df <- bind_cols(pred = predvals, mod_df)

mod_df %>% filter(trt == "Control")
mod_df %>% filter(trt == "K")
mod_df %>% filter(trt == "N")

mod_df %>% 
  ggplot(aes(x = exp(year_trt) - 1,
                      y = exp(pred),
                      color = trt)) +
  stat_smooth() +
  geom_point()

est_coefs <- bayescoefs$quantiles[1:9,3]

names(est_coefs) <- c("Intercept", levels(mod_df$trt))

# Predicted changes in variance
sigma_coefs <- bayescoefs$quantiles[c(18, 10:17),3]
sigma_ests <- modmat %*% sigma_coefs
mod_df <- bind_cols(mod_df, sigma_ests = c(sigma_ests))
mod_df %>%
  ggplot(aes(x = exp(year_trt) - 1,
             y = sigma_ests,
             color = trt)) +
  stat_smooth() +
  geom_point()

mod_df %>%
  ggplot(aes(x = exp(year_trt) - 1,
             y = sqrt((exp(pred)^2) * (exp(sigma_ests ^ 2) - 1)),
             color = trt)) +
  stat_smooth() +
  geom_point()

mod_df <- mod_df %>% mutate(sd_modeled = c(sqrt((exp(pred)^2) * (exp(sigma_ests ^ 2) - 1))))

mod_df %>%
  ggplot(aes(x = exp(year_trt) - 1,
             y = pred / sd_modeled,
             color = trt)) +
  stat_smooth() +
  geom_point() +
  xlab("Year Trt") +
  ylab("Inverse CV")



require("lmvar"); require(MASS)

dat <- dat[!is.na(log(dat$live_mass))]
live_mass <- log(dat$live_mass)
year_trt <- log(dat$year_trt + 1)
site_code <- as.factor(dat$site_code)

varvar <- lmvar(y = live_mass,
      X_mu = year_trt + as.integer(dat$trt), 
      X_sigma = year_trt)

summary(varvar)
