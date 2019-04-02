source("R/makefood.R")
source("R/eatfood.R")

num_sources <- 3
popsize <- c(10,100,1000)
mu_carb <- c(1, 100, 1000)
var_carb <- c(0,1,10)
diet_prop <- c(.1, .25, .65)

food <- makefood(3, popsize, mu_carb, var_carb)

ind_eatfood(num_sources, diet_prop, food)
