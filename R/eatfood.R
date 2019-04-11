#' eatfoodonce
#'
#'Simulates a single consumption event by a consumer and returns the carbon and
#' nitrogen values for the consumed prey.  Meant to be call by `eatfoodsteps`
#' and `eatfood`.
#'
#' @param num_sources integer, the number of food sources to choose from.
#' Can not exceed the actual number of available food sources.
#' @param diet_prop vector, the relative probability of selecting a source.
#' Must have an entry for each source.
#' @param food list, the simulated food sources from which to select.
#'
#' @return A vector containing the carbon, then nitrogen values for the prey
#' item selected.
#' @export
#'
#' @examples
eatfoodonce <- function(num_sources, diet_prop, food) {
  #check that num_sources doesn't exceed the simulated available sources
  if(num_sources > length(food))
    stop("num_sources can not exceed the number of simulated food sources")

  #check that there is a probability for each source
  if(num_sources != length(diet_prop))
    stop("Must supply a source selection probability for each source (length of diet_prop must equal num_sources")

  #choose a diet source population
  source <- sample(x = c(num_sources), size = 1, prob = diet_prop)
  #choose an individual from that population
  eaten <- c(sample(food[[source]][[1]], size = 1), #pull a carbon value
                sample(food[[source]][[2]], size = 1)) #pull a nitrogen value
  return(eaten)
}


#' eatfoodsteps
#'
#' Wrapper function to interate `eatfoodonce` across multiple time steps (for
#' a sinlge individual). Returns a matrix of C and N values sampled by the
#' individual over `steps` number of steps.
#'
#' @param steps Integer, the number of time steps to simulate
#' @param num_sources integer, the number of food sources to choose from.
#' Can not exceed the actual number of available food sources. Passed to
#' `eatfoodonce`.
#' @param diet_prop vector, the relative probability of selecting a source.
#' Must have an entry for each source. Passed to `eatfoodonce`.
#' @param food list, the simulated food sources from which to select. Passed to `eatfoodonce`.
#'
#' @return A list of consumer-sampled prey isotope values. Each element of the
#' list is a step containing the C and N values for the prey item sampled.
#' @export
#'
#' @examples
eatfoodsteps <- function(steps, num_sources, diet_prop, food){
  foodtrack <- replicate(n = steps, eatfoodonce(num_sources = num_sources,
                                                diet_prop = diet_prop,
                                                food = food), simplify = F)
  return(foodtrack)
}

#' eatfood
#'
#'Wrapper function to iterate `eatfoodsteps` across multiple individuals.
#' @param num_individ Integer, the number of individuals to simulate
#' @param steps Integer, the number of time steps to simulate. Passed to
#' `eatfoodsteps`.
#' @param num_sources integer, the number of food sources to choose from.
#' Can not exceed the actual number of available food sources. Passed to
#' `eatfoodsteps` and then to `eatfoodonce`.
#' @param diet_prop vector, the relative probability of selecting a source.
#' Must have an entry for each source. Passed to `eatfoodsteps` and then to
#' `eatfoodonce`.
#' @param food list, the simulated food sources from which to select. Passed to
#' `eatfoodsteps` and then to `eatfoodonce`.
#'
#' @return A list, each element of which represents an individual consumer and
#'  contains the time step matrix supplied by `eatfoodsteps`.
#' @export
#'
#' @examples
eatfood <- function(num_individ, steps, num_sources, diet_prop, food) {
  diethistory <- replicate(n = num_individ, eatfoodsteps(steps = steps,
                                                     num_sources = num_sources,
                                                     diet_prop = diet_prop,
                                                     food = food),
                           simplify = F) #force the result to be a list
  return(diethistory)
}

