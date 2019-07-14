# GIRTmodel study

# 現状のGIRTモデルの問題点
# 1, The interpretation of phi parameter is difficult (has not been solved).
# 2, Unstability of phi parameters estimati

library(irtfun2)
library(tidyverse)
library(extraDistr)

set.seed(0204)
theta <- rnorm(3000)
phi <- rinvchi(3000, max = 2)
a <- rlnorm(1000, sdlog = 0.25)
b <- rnorm(1000)
dat1 <- sim_gen(theta=theta, phi=phi, a=a, b=b)

fit1 <- estGip(dat1, fc = 2, Nphi = 21, esteap = T, min_ph = 0.001, max_ph = 4, phi_dist = "gbeta")
fit1$item

fit1$person$phi %>% hist()
cor(fit1$person$phi, phi)
plot(fit1$person$phi, phi)


# fix phi parameter to 1
set.seed(0204)
theta <- rnorm(3000)
phi <- rep(1, 3000)
a <- rlnorm(30, sdlog = 0.25)
b <- rnorm(30)
dat2 <- sim_gen(theta=theta, phi=phi, a=a, b=b)
fit2 <- estGip(dat2, fc = 2, Nphi = 21, esteap = T, min_ph = 0.001, max_ph = 4, phi_dist = "gbeta")
fit2$item
plot(fit2$item$a, a)
plot(fit2$item$b, b)
plot(fit2$person$theta, theta)
plot(fit2$person$phi, phi)
hist(fit2$person$phi)
mean(abs(fit2$person$phi -  phi))


# compare theta estimate
set.seed(0204)
theta <- rnorm(3000)
# phi <- rep(1, 3000)
a <- rlnorm(30, sdlog = 0.25)
b <- rnorm(30)
dat3 <- sim_gen(theta=theta, phi=phi, a=a, b=b)
fit3 <- estGip(dat3, fc = 2, Nphi = 21, esteap = T, min_ph = 0.001, max_ph = 4, phi_dist = "gbeta")
fit4 <- estip2(dat3, fc = 2)
fit4$preson <- estheta(dat3, fit4$para, est = "EAP", fc = 2, gc = 0)

plot(theta, fit3$person$theta); cor(theta, fit3$person$theta) # true vs GIRT
plot(theta, fit4$preson$res$EAP); cor(theta, fit4$preson$res$EAP) # true vs IRT
plot(fit3$person$theta, fit4$preson$res$EAP) # GIRT vs IRT
# item parameter 
plot(a, fit3$item$a, xlim = c(0.5, 1.5), ylim = c(0.5, 1.5)); cor(a, fit3$item$a) # GIRT
plot(a, fit4$para$a, xlim = c(0.5, 1.5), ylim = c(0.5, 1.5)); cor(a, fit4$para$a) # IRT
plot(fit3$item$a, fit4$para$a, xlim = c(0.5, 1.5), ylim = c(0.5, 1.5))


# comprehension of phi parameter
tibble(phi = 0:100) %>% ggplot(aes(x = phi)) + stat_function(fun = gptheta, args = list(theta = -1, b = 0, a = 1, D = 1.702))


LLgirt <- function(para, u, a, b, D, Bayes = FALSE){
  theta <- para[1]
  phi <- para[2]
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  sum(log(p)*u + log(1-p)*(1-u)) + Bayes*(log(dnorm(theta,0,1)) + log(dbetapr(phi, shape1 = 2, shape2 = 2)))
}
LLgirt2 <- function(theta, phi, u, a, b, D){
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  sum(log(p)*u + log(1-p)*(1-u))
}
LLgirt_apply <- function(theta, phi, u, a, b, D){
  mapply(FUN = LLgirt2, theta, phi, MoreArgs = list(u = u, a = a, b = b, D = D), SIMPLIFY = T) %>% as.vector
}

LLgrit_gr <- function(para, u, a, b, D){ # incorrect
  theta <- para[1]
  phi <- para[2]
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  A <- D*a/(sqrt(1+phi^2*a^2))
  B <- D*a^3*phi*(theta-b) / (1+phi^2*a^2)^(3/2)
  thetagr <- sum(A*(u-p))
  phigr <- sum(B*(u-p))
  c(thetagr, phigr)
}

optim(par = c(0, 1), fn = LLgirt, u = c(1,0,1,1), a = c(1,1,1,1), b = c(-1,0,1,2), D = 1.0, method = "BFGS", hessian = T, control = list(fnscale = -1))

optim(par = c(0, 1), fn = LLgirt, gr = LLgrit_gr, u = c(1,0,1,1), a = c(1,1,1,1), b = c(-1,0,1,2), D = 1.0, method = "BFGS", hessian = T, control = list(fnscale = -1))

optim(par = c(0, 1), fn = LLgirt, gr = LLgrit_gr, u = c(1,1,1,0,1), a = c(1,1,1,1,1), b = c(-2,-1,0,1,2), D = 1.0, method = "BFGS", hessian = T, control = list(fnscale = -1))
optim(par = c(0, 1), fn = LLgirt, u = c(1,0,1,1,1), a = c(1,1,1,1,1), b = c(-2,-1,0,1,2), D = 1.0, hessian = T, control = list(fnscale = -1))

# visualization
test <- tibble(theta = seq(-4, 4, length.out = 100) %>% rep(100), phi = apply(matrix(seq(0, 4, length.out = 100)), 1, rep, 100) %>% as.vector) %>% 
  mutate(LL = LLgirt_apply(theta = theta, phi = phi, 
                           u = c(0,0,1,1,1), 
                           a = c(1,1,1,1,1), b = c(-2,-1,0,1,2), D = 1.0))

test %>% ggplot(aes(x = theta, y = phi, z = LL)) + geom_contour(binwidth = 0.05) + xlim(1, 4) + ylim(0, 4)
test %>% ggplot(aes(x = theta, y = phi, z = LL)) + geom_raster(aes(fill = LL), hjust = 0, vjust = 0, interpolate = F)

girt_llplot <- function(u, a, b, D, Bayes = FALSE, binwidth = 0.1){
  # estimate pars
  opt <- optim(par = c(0, 1), fn = LLgirt, 
              # gr = LLgrit_gr, 
               u = u, a = a, b = b, D = D, Bayes = Bayes, method = "BFGS", hessian = T, control = list(fnscale = -1))
  cat(opt$par)
  if(opt$par[2] < 0) warning("Failuer to optimise.")
  par <- tibble(theta = opt$par[1], phi = opt$par[2], LL = opt$value)
  dat <- tibble(theta = seq(-4, 4, length.out = 100) %>% rep(100), phi = apply(matrix(seq(0, 4, length.out = 100)), 1, rep, 100) %>% as.vector) %>% 
    mutate(LL = LLgirt_apply(theta = theta, phi = phi, u = u, a = a, b = b, D = D))
  plt <- dat %>% ggplot(aes(x = theta, y = phi, z = LL)) + geom_contour(binwidth = binwidth) + xlim(-4, 4) + ylim(0, 4)
  plt + geom_point(data = par, aes(x = theta, y = phi))
}
girt_llplot(u = c(0,1,1,0,1,1,0,1,0,0), 
            a = c(1,1,1,1,1,1,1,1,1,1), 
            b = c(-2,-1,0,1,2,-2,-1,0,1,2), D = 1.0, binwidth = 0.1, Bayes = TRUE)
