#test out make food
source("R/makefood.R")
num_sources <- 3
popsize <- c(10,100,1000)
mu_carb <- c(1, 100, 1000)
var_carb <- c(0,1,10)
makefood(3, popsize, mu_carb, var_carb)
