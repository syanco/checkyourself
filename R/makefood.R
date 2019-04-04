#' Make food
#'
#' Simulates the 13C isotope values for an arbitray number of potential diet
#' endpoints. Each endpoint population is simulated as a normally distributed
#' variable
#'
#' @param num_sources integer, number of potential food sources to simulate
#' @param popsize vector whose length = 'num_sources', each entry of which is
#' the size of the diet source population to simulate
#' @param mu_carb vector whose length = 'num_sources', each entry of which is
#' the mean 13C value of the diet source population to simulate
#' @param sd_carb vector whose length = 'num_sources', each entry of which is
#' the standard deviation of the 13C of the diet source population to simulate
#'#' @param mu_nit vector whose length = 'num_sources', each entry of which is
#' the mean 14N value of the diet source population to simulate
#' @param sd_nit vector whose length = 'num_sources', each entry of which is
#' the standard deviation of the 14N of the diet source population to simulate
#'
#' @return A list, each element of which is a vector of isotope vales for each
#' simulated diet endpoint individual
#' @export
#'
#' @examples
makefood <- function(num_sources, popsize, mu_carb, sd_carb, mu_nit, sd_nit) {
  food <- list()
  for(i in 1:num_sources){
    food[[i]] <- list(rnorm(n = popsize[i], mean = mu_carb[i], sd = sd_carb[i]),
                   rnorm(n = popsize[i], mean = mu_nit[i], sd = sd_nit[i]))
    }
  return(food)
  }

