# GIRT simulation
# Condition----
# 1, Assume 2PLM
# 2, Assume GIRT MODEL
# 3, Assume 2 dimentional IRT (Non Compensatory)

# using pkg----
library(tidyverse)
library(irtfun2)
library(mirt)

# Simdat----

t <- 0 # initialize count
model <- "2PL"
mod <- mirt.model("F1 = 1-30")
seed1 <- numeric(100)
while(t < 100){
  t <- t + 1
  seed1[t] <- round(runif(1) * 10000)
  set.seed(seed1[t])
  at <- rlnorm(30, -0.2, 0.5) # non negetive real
  bt <- rnorm(30, 0, 1.5) # real
  tt <- rnorm(5000, 0, 1) # real 
  data <- simdata(a = at, d = -bt*at, Theta = as.matrix(tt), itemtype = "dich") # b = -d/a then d = -a*b
  # 2PL estimation
  fit_2pl <- mirt(data, mod, itemtype = model, quadpts = 31)
  par_2pl <- coef(fit_2pl, IRTpars = T, simplify = T)$items[,c("a", "b")] %>% as.data.frame
  # GIRT estimation
  fit_gpl <- estGip(data, fc = 1, IDc = 1, method = "Fisher_Scoring", phi_dist = "lognormal", min_ph = 0.1, max_ph = 2, sigma_ph = 1, maxiter_em = 500, eEM = 1e-6, eDIST = 1e-6)
  par_gpl <- fit_gpl$item %>% transmute(a = 1.702 * a, b = b)
  
plot(at, par_gpl$a)
  plot(at, par_2pl$a)
  
}
