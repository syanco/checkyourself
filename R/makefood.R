
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
#' @param var_carb vector whose length = 'num_sources', each entry of which is
#' the standard deviation of the diet source population to simulate
#'
#' @return A list, each element of which is a vector of isotope vales for each
#' simulated diet endpoint individual
#' @export
#'
#' @examples
makefood <- function(num_sources, popsize, mu_carb, var_carb) {
  food <- list()
  for(i in 1:num_sources){
    food[[i]] <- rnorm(n = popsize[i], mean = mu_carb[i], sd = var_carb[i])
  }
  return(food)
}

