# GIRT simulation
# Condition----
# 1, Assume 2PLM
# 2, Assume GIRT MODEL
# 3, Assume 2 dimentional IRT (Non Compensatory)

devtools::install_github("takuizum/irtfun2")
# using pkg----
library(tidyverse)
library(irtfun2)
library(mirt)
library(extraDistr)

# using function ----
rmse <- function(x, y = list()){
  n <- length(as.vector(x))
  z <- numeric(length(y))
  for(i in 1:length(y)){
    yy <- y[[i]]
    z[i] <- sqrt(sum((yy - x)^2)/n)
  }
  return(z)
}
r2norm <- function(n, mu, sigma, rho) {
  tmp <- rnorm(n)
  x   <- mu+sigma*tmp
  y   <- rho*x + sqrt(1-rho^2)*rnorm(n)
  return(as.matrix(data.frame(th1=x,th2=y)))
}

# phi prior distribution ----
# log normal
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dlnorm, args = list(meanlog = -0.5, sdlog = 0.5))
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dlnorm, args = list(meanlog = 0, sdlog = 0.25))
# invchi 
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dinvchi, args = list(v = 0.1, tau = 0.1))
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dinvchi, args = list(v = 1, tau = 1))
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dinvchi, args = list(v = 5, tau = 3))
# inv chi sq
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dinvchisq, args = list(nu = 1, tau = 1))
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dinvchisq, args = list(nu = 5, tau = 3))
# chi square
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dchisq, args = list(df = 1, ncp = 1))
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dchisq, args = list(df = 2, ncp = 2))
# beta dist
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dbeta, args = list(shape1 = 1, shape2 = 3))
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dbeta, args = list(shape1 = 1, shape2 = 2))
# beta prime 
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dbetapr, args = list(shape1 = 1, shape2 = 5))
tibble(x = 0.001:4) %>% ggplot(aes(x = x)) + stat_function(fun = dbetapr, args = list(shape1 = 2, shape2 = 3))


estGip(sim_gen(rnorm(5000), rlnorm(30, 0, 0.25), rnorm(30, 0, 1.5)), fc = 2, phi_dist = "gbeta")
estGip(sim_gen(rnorm(5000), rlnorm(30, 0, 0.25), rnorm(30, 0, 1.5), phi = rlnorm(5000, 0, 0.25)), fc = 2)

# Simdat----

# the data follows 2 Parameter Logistic Model ----
t <- 0 # initialize count
model <- "2PL"
mod <- mirt.model("F1 = 1-30")
seed1 <- numeric(100)
girtpar1 <- list_along(numeric(100))
irtpar1 <- list_along(numeric(100))
result_obj1 <- list_along(numeric(100))
while(t < 100){
  t <- t + 1
  cat(t, "time simulation!\n")
  seed1[t] <- round(runif(1) * 10000)
  set.seed(seed1[t])
  at <- rlnorm(30, -0.5, 0.5) # non negetive real
  bt <- rnorm(30, 0, 1.5) # real
  tt <- rnorm(5000, 0, 1) # real
  data <- simdata(a = at, d = -bt*at, Theta = as.matrix(tt), itemtype = "dich") # b = -d/a then d = -a*b
  # GIRT estimation
  fit_gpl <- try(estGip(data, fc = 1, IDc = 0, phi_dist = "gbeta", 
                        min_ph = 0, max_ph = 4, eEM = 0.0001, Nphi = 10, D = 1.0, maxiter_em = 200))
  if(class(fit_gpl) == "try-error"){
    t <- t - 1
    next
  }
  par_gpl <- fit_gpl$item %>% transmute(a = 1.0 * a, b = b)
  # 2PL estimation
  fit_2pl <- mirt(data, mod, itemtype = model, quadpts = 31, TOL = 0.0001, accelerate = "none")
  par_2pl <- coef(fit_2pl, IRTpars = T, simplify = T)$items[,c("a", "b")] %>% as.data.frame
  
  # plot
  # discrimination
  par(mfrow=c(2,2)) 
  plot(at, par_gpl$a, ylim = c(0, 3), xlim = c(0, 3), xlab = "GIRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  plot(at, par_2pl$a, ylim = c(0, 3), xlim = c(0, 3), xlab = "IRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  # difficulty
  plot(bt, par_gpl$b, ylim = c(-4, 4), xlim = c(-4, 4), xlab = "GIRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  plot(bt, par_2pl$b, ylim = c(-4, 4), xlim = c(-4, 4), xlab = "IRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  # rmse calc
  rmse_a <- rmse(x = at, y = list(par_gpl$a, par_2pl$a))
  names(rmse_a) <- c("girt", "irt")
  rmse_b <- rmse(x = bt, y = list(par_gpl$b, par_2pl$b))
  names(rmse_b) <- c("girt", "irt")
  # AIC
  aic <- c("girt" = -2*last(fit_gpl$mll_history) + 2*30*2, "irt" = fit_2pl@Fit$AIC)
  # n of iter
  iter <- c("girt" = length(fit_gpl$mll_history) - 1, "irt" = fit_2pl@OptimInfo$converged)
  # contain the estiamte data in list
  result_obj1[[t]] <- list("rmse_a" = rmse_a, "rmse_b" = rmse_b, "girtpar" = par_gpl, "irtpar" = par_2pl, "aic" = aic, "iter" = iter)
}

# the data follows General IRT model ----
t <- 0 # initialize count
model <- "2PL"
mod <- mirt.model("F1 = 1-30")
seed3 <- numeric(100)
girtpar3 <- list_along(numeric(100))
irtpar3 <- list_along(numeric(100))
result_obj3 <- list_along(numeric(100))
while(t < 100){
  t <- t + 1
  cat(t, "time simulation!\n")
  seed3[t] <- round(runif(1) * 10000)
  set.seed(seed3[t])
  at <- rlnorm(30, -0.5, 0.5) # non negetive real
  bt <- rnorm(30, 0, 1.5) # real
  tt <- rnorm(5000, 0, 1) # real
  pt <- rbetapr(5000, shape1 = 1, shape2 = 3)
  data <- sim_gen(a = at, b = bt, theta = as.matrix(tt), phi = pt, D = 1.0)[,-1] # b = -d/a then d = -a*b
  # GIRT estimation
  fit_gpl <- try(estGip(data, fc = 1, IDc = 0, phi_dist = "gbeta", 
                        min_ph = 0, max_ph = 4, eEM = 0.0001, Nphi = 10, D = 1.0, maxiter_em = 200))
  if(class(fit_gpl) == "try-error"){
    t <- t - 1
    next
  }
  par_gpl <- fit_gpl$item %>% transmute(a = 1.0 * a, b = b)
  # 2PL estimation
  fit_2pl <- mirt(data, mod, itemtype = model, quadpts = 31, TOL = 0.0001, accelerate = "none")
  par_2pl <- coef(fit_2pl, IRTpars = T, simplify = T)$items[,c("a", "b")] %>% as.data.frame
  
  # plot
  # discrimination
  par(mfrow=c(2,2)) 
  plot(at, par_gpl$a, ylim = c(0, 3), xlim = c(0, 3), xlab = "GIRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  plot(at, par_2pl$a, ylim = c(0, 3), xlim = c(0, 3), xlab = "IRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  # difficulty
  plot(bt, par_gpl$b, ylim = c(-4, 4), xlim = c(-4, 4), xlab = "GIRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  plot(bt, par_2pl$b, ylim = c(-4, 4), xlim = c(-4, 4), xlab = "IRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  # rmse calc
  rmse_a <- rmse(x = at, y = list(par_gpl$a, par_2pl$a))
  names(rmse_a) <- c("girt", "irt")
  rmse_b <- rmse(x = bt, y = list(par_gpl$b, par_2pl$b))
  names(rmse_b) <- c("girt", "irt")
  # AIC
  aic <- c("girt" = -2*last(fit_gpl$mll_history) + 2*30*2, "irt" = fit_2pl@Fit$AIC)
  # n of iter
  iter <- c("girt" = length(fit_gpl$mll_history) - 1, "irt" = fit_2pl@OptimInfo$converged)
  # contain the estiamte data in list
  result_obj3[[t]] <- list("rmse_a" = rmse_a, "rmse_b" = rmse_b, "girtpar" = par_gpl, "irtpar" = par_2pl, "aic" = aic, "iter" = iter)
}

# the data follows 2 dimensional IRT model ----
t <- 0 # initialize count
model <- "2PL"
mod <- mirt.model("F1 = 1-30")
modm <- mirt.model("F1 = 1-30
                    F2 = 1-30
                    COV = F1*F2")
seed2 <- numeric(100)
girtpar2 <- list_along(numeric(100))
irtpar2 <- list_along(numeric(100))
result_obj2 <- list_along(numeric(100))
while(t < 100){
  t <- t + 1
  cat(t, "time simulation!\n")
  seed2[t] <- round(runif(1) * 10000)
  set.seed(seed2[t])
  a1t <- rlnorm(30, -0.5, 0.5) # non negetive real
  a2t <- runif(30, min = 0.5, max = 0.8) * a1t
  at <- cbind(a1 = a1t, a2 = a2t)
  bt <- rnorm(30, 0, 1.5) # real
  tt <- r2norm(5000, 0, 1, rho = 0.7) # real
  data <- simdata(a = at, d = -bt*(a1t^2+a2t^2), Theta = tt, itemtype = "dich") # b = -d/a then d = -a*b
  # GIRT estimation
  fit_gpl <- try(estGip(data, fc = 1, IDc = 0, phi_dist = "gbeta", 
                        min_ph = 0, max_ph = 4, eEM = 0.0001, Nphi = 10, D = 1.0, maxiter_em = 200))
  if(class(fit_gpl) == "try-error"){
    t <- t - 1
    next
  }
  par_gpl <- fit_gpl$item %>% transmute(a = 1.0 * a, b = b)
  # 2PL estimation
  fit_2pl <- mirt(data, mod, itemtype = model, quadpts = 31, TOL = 0.0001, iter = 200)
  par_2pl <- coef(fit_2pl, IRTpars = T, simplify = T)$items[,c("a", "b")] %>% as.data.frame
  # 2 dim IRT estimation
  fit_mpl <- mirt(data, modm, itemtype = model, quadpts = 31, TOL = 0.0001, accelerate = "none", rotate = "none")
  par_mpl <- coef(fit_mpl, IRTpars = F, simplify = T)$items[,c("a1", "a2", "d")] %>% as.data.frame %>% mutate(b = -d/(a1^2+a2^2))
  
  # plot
  # discrimination
  par(mfrow=c(2,3)) 
  plot(a1t, par_gpl$a, ylim = c(0, 3), xlim = c(0, 3), xlab = "GIRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  plot(a1t, par_2pl$a, ylim = c(0, 3), xlim = c(0, 3), xlab = "IRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  plot(a1t, par_mpl$a1, ylim = c(0, 3), xlim = c(0, 3), xlab = "2DMIRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  
  # difficulty
  plot(bt, par_gpl$b, ylim = c(-4, 4), xlim = c(-4, 4), xlab = "GIRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  plot(bt, par_2pl$b, ylim = c(-4, 4), xlim = c(-4, 4), xlab = "IRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  plot(bt, par_mpl$b, ylim = c(-4, 4), xlim = c(-4, 4), xlab = "2DMIRT_estimates", ylab = "true");abline(a=0, b = 1, col='red')
  
  # rmse calc
  rmse_a <- rmse(x = at, y = list(par_gpl$a, par_2pl$a))
  names(rmse_a) <- c("girt", "irt")
  rmse_b <- rmse(x = bt, y = list(par_gpl$b, par_2pl$b))
  names(rmse_b) <- c("girt", "irt")
  # AIC
  aic <- c("girt" = -2*last(fit_gpl$mll_history) + 2*30*2, "irt" = fit_2pl@Fit$AIC)
  # n of iter
  iter <- c("girt" = length(fit_gpl$mll_history) - 1, "irt" = fit_2pl@OptimInfo$converged)
  # contain the estiamte data in list
  result_obj1[[t]] <- list("rmse_a" = rmse_a, "rmse_b" = rmse_b, "girtpar" = par_gpl, "irtpar" = par_2pl, "aic" = aic, "iter" = iter)
}