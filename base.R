# GIRTmodel study

# 現状のGIRTモデルの問題点
# 1, The interpretation of phi parameter is difficult (has not been solved).
# 2, Unstability of phi parameters estimati

library(irtfun2)
library(tidyverse)

res <- estGip(sim_dat_girt, fc = 2, Nphi = 21, esteap = T, min_ph = 0.001, max_ph = 4, phi_dist = "gbeta")
res$item

res$person$phi %>% hist()
