# InDiscPackage.R

install.packages("InDisc")
library(InDisc)
library(mirt)
library(tidyverse)

data <- expand.table(LSAT6)
fit1 <- InDisc(data, model = "graded")
fit1$INDIES %>% head(30)

devtools::install_github("takuizum/irtfun2")
library(irtfun2)
fit2 <- estGip(data, fc = 1, esteap = TRUE)
fit2$person