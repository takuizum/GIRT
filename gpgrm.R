# General IRT + GRM


# response probability----
gpgrm <- function(theta, phi, a, b, k, D = 1.0){
  K <- length(b) + 1 # n of category(n of b parameter + 1)
  if(k == 1){
    A <- sqrt(1+phi^2*a^2)
    e <- exp(-D*a/A*(theta-b[k]))
    p1 <- 1/(1+e)
    p0 <- 1 # category 0
  } else if (k == K){
    p1 <- 0 # category K + 1
    A <- sqrt(1+phi^2*a^2)
    e <- exp(-D*a/A*(theta-b[k-1]))
    p0 <- 1/(1+e)
  } else if (k < K && k > 1) {
    A <- sqrt(1+phi^2*a^2)
    e1 <- exp(-D*a/A*(theta-b[k]))
    e0 <- exp(-D*a/A*(theta-b[k-1]))
    p1 <- 1/(1+e1)
    p0 <- 1/(1+e0)
  } else { # NA
    p0 <- 1
    p1 <- 0
  }
  p0-p1
}

bt <- c(-1.5, 0, 1)
# pgrm(0, a = 2, b = bt, k = 1)
tibble(theta = c(-4:4)) %>%
  ggplot(aes(x = theta))+
  stat_function(fun = gpgrm, args = list(a = 2, b = bt, k = 1, D = 1, phi = 100), colour = 2)+
  stat_function(fun = gpgrm, args = list(a = 2, b = bt, k = 2, D = 1, phi = 100), colour = 3)+
  stat_function(fun = gpgrm, args = list(a = 2, b = bt, k = 3, D = 1, phi = 100), colour = 4)+
  stat_function(fun = gpgrm, args = list(a = 2, b = bt, k = 4, D = 1, phi = 100), colour = 6)
