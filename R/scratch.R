source("R/makefood.R")
source("R/eatfood.R")
source("R/gettissueiso.R")

num_sources <- 3
popsize <- c(1000,1000,1000)
mu_carb <- c(1, 100, 1000)
var_carb <- c(0,1,10)
diet_prop <- c(.1, .25, .65)
steps <- 100
num_individ <- 10

food <- makefood(3, popsize, mu_carb, var_carb)

specimens <- eatfood(num_individ, steps, num_sources, diet_prop, food)

getiso(specimens, time = 100, window =4)
