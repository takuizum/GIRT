library(tidyverse)
library(rstan)
library(mirt)

# sim data
set.seed(1234)
a <- rlnorm(10)
d <- rnorm(10, sd = 2)
data <- mirt::simdata(a, d, N = 100, itemtype = "dich")
# Set data
mod <- stan_model("src/bin_girt_cor.stan")
lst <- list(y = data, N = 100, M = 10)
# Inference
smp <- vb(mod, lst)
options(mc.cores=parallel::detectCores())
smp <- sampling(mod, lst, chains = 2)
shinystan::launch_shinystan(smp)