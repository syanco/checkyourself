#' iteratesourcesamp
#'
#'Function to run entire trophic iso-ecology simulation iterating across sample
#'sizes for sample-based estimation of prey source isotopic distribution
#'estimation.  Requires the designation of a folder for all data output. Calls
#'all necessary functions to run a complete simulation and saves outputs into
#'format ready for analysis with `MixSIAR`.
#'
#' @param folder String, path to folder where simulated data should be saved
#' @param num_sources integer, number of potential food sources to simulate
#' @param popsize vector whose length = 'num_sources', each entry of which is
#' the size of the diet source population to simulate
#' @param mu_carb vector whose length = 'num_sources', each entry of which is
#' the mean 13C value of the diet source population to simulate
#' @param sd_carb vector whose length = 'num_sources', each entry of which is
#' the standard deviation of the 13C of the diet source population to simulate
#'#' @param mu_nit vector whose length = 'num_sources', each entry of which is
#' the mean 15N value of the diet source population to simulate
#' @param sd_nit vector whose length = 'num_sources', each entry of which is
#' the standard deviation of the 15N of the diet source population to simulate
#' @param sourcesamp Vector of sample sizes through which the source sampling
#' walks
#' @param num_individ Integer, the number of individuals to simulate
#' @param steps Integer, the number of time steps to simulate. Passed to
#' `eatfoodsteps`.
#' @param diet_prop vector, the relative probability of selecting a source.
#' Must have an entry for each source. Passed to `eatfoodsteps` and then to
#' `eatfoodonce`.
#' @param window Integer, the consumer isotope integration window representing
#' the number of steps across which prey isotope values are integrated in an
#' individual consumer to generate consumer tissue isotope values. Window must
#' be less than or equal to `time`.
#'
#' @return Writes simulated consumer mixture distribution file, sample-based
#' prey source estimate files for each value of `sourcesamp`, and a trophic
#' discrimination factor file all in format that can be used by `MixSIAR`.
#' @export
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

  #sample from prey source populations walking across parameter space of sample
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
