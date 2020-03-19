library(tidyverse)
library(rstan)
library(mirt)

# sim data
set.seed(1234)
a <- rlnorm(10)
d <- rnorm(10, sd = 2)
data <- mirt::simdata(a, d, N = 100, itemtype = "dich")
# Set data
mod <- stan_model("src/girtmodel.stan")
lst <- list(y = data, N = 100, M = 10)
# Inference
smp <- vb(mod, lst)
options(mc.cores=parallel::detectCores())
smp <- sampling(mod, lst, chains = 4, iter = 5000, thin = 4)
shinystan::launch_shinystan(smp)

# marginal model
mod2 <- stan_model("src/marginal_model.stan")
lst2 <- list(y = data, N = 100, P = 10, K = 17)
smp2 <- vb(mod2, lst2)
smp <- sampling(mod2, lst2, chains = 4, iter = 2000, thin = 1)
shinystan::launch_shinystan(smp)

# correlated model
mod3 <- stan_model("src/bin_girt_cor.stan")
lst <- list(y = data, N = 100, M = 10)
# Inference
smp3 <- vb(mod3, lst)
options(mc.cores = parallel::detectCores())
smp3 <- sampling(mod3, lst, chains = 4, iter = 2000, thin = 4)
shinystan::launch_shinystan(smp3)