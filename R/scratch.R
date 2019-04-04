source("R/makefood.R")
source("R/eatfood.R")
source("R/gettissueiso.R")

num_sources <- 3
popsize <- c(1000,1000,1000)
mu_carb <- c(1, 100, 1000)
sd_carb <- c(0,1,10)
mu_nit <- c(1, 100, 1000)
sd_nit <- c(0,1,10)


diet_prop <- c(.1, .25, .65)
steps <- 100
num_individ <- 100

#make the food sources
food <- makefood(3, popsize, mu_carb, sd_carb, mu_nit, sd_nit)

#simulate agents over time
specimens <- eatfood(num_individ, steps, num_sources, diet_prop, food)

#collect "observed" isotop values
obs_iso <- getiso(specimens, time = steps, window = 10)

#format to data frame ans save as csv for import to MixSIAR


#sample analysis
library(MixSIAR)
