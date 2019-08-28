library(tidyverse)

norm <- rnorm(n = 101, sd = 1, mean = 10)
increasingmean <- rnorm(n = 101, sd = seq(1, 2, by = .01), mean = seq(10, 20, by = .1))
constantlarger <- rnorm(n = 101, sd = 2, mean = 20)
increasingmean_samevar <- rnorm(n = 101, sd = seq(1, .5, by = -.01), mean = seq(10, 20, by = .1))

data.frame(i = c(1:101),
           norm,
           increasingmean, constantlarger, increasingmean_samevar) %>%
  gather(key = "source", value = "val", -i) %>%
  ggplot(aes(x = i, y = val, color = source)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  facet_wrap(~source)

mean(norm) / sd(norm)
mean(increasingmean) / sd(increasingmean)
mean(constantlarger) / sd(constantlarger)


norm <- rnorm(n = 101, sd = 1, mean = 10)
increasingmean <- rnorm(n = 101, sd = seq(1, 2, by = .01), mean = seq(10, 20, by = .1))

norm <- rnorm(n = 1001, sd = 5, mean = 10)
increasingmean_1fold <- rnorm(n = 1001, sd = 5, mean = seq(10, 20, by = .01))
increasingmean_2fold <- rnorm(n = 1001, sd = 5, mean = seq(10, 30, by = .02))
increasingmean_3fold <- rnorm(n = 1001, sd = 5, mean = seq(10, 50, by = .04))

log(sd(norm))
log(sd(increasingmean_1fold))
log(sd(increasingmean_2fold))
log(sd(increasingmean_3fold))






mean(increasingmean) / sd(increasingmean - seq(10, 20, by = .1))

sd(increasingmean_samevar)
sd(increasingmean_samevar - seq(10, 20, by = .1))

mean(increasingmean_samevar) / sd(increasingmean_samevar)
mean(increasingmean_samevar) / sd(increasingmean_samevar - seq(10, 20, by = .1))

