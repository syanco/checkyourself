iteratesourcesamp <- function (folder, num_sources, popsize, mu_carb, sd_carb,
                               mu_nit, sd_nit, sourcesamp, num_individ, steps,
                               diet_prop, window){
  wd <- getwd()#get current wd
  #change working directory to designated data folder
  setwd(folder)

  #simulate prey sources
  preysource <- makefood(num_sources = num_sources, popsize = popsize,
                   mu_carb = mu_carb, sd_carb = sd_carb, mu_nit = mu_nit,
                   sd_nit = sd_nit)

  #sample from pery source populations walking across parameter space of sample
  #sizes.  Writes source iso estimates to file
  for(i in sourcesamp){
    samplesources(num_samples = i, food = preysource,
                  filename = paste0("simulated_source_iso_n", i,
                                    ".csv"))
  }
  #
  #
  # sapply(sourcesamp, FUN = samplesources, food = preysource,
  #        filename = paste0("simulated_source_iso_n", sourcesamp, ".csv"))

  #save isotopic discrimination factors to file
  savediscrimination(food = preysource, filename="simulated_discrimination.csv")

  #simulate agents over time
  specimens <- eatfood(num_individ =  num_individ, steps = steps,
                       num_sources = num_sources, diet_prop = diet_prop,
                       food = preysource)

  #collect "observed" isotope values
  obs_iso <- getiso(food = specimens, time = steps, window = window)

  #format to data frame and save as csv for import to MixSIAR
  formatiso(iso = obs_iso, filename = "simulated_consumer_iso.csv")

  #return wd to previous
  setwd(wd)
}


