#' eatfoodonce
#'
#' @param num_sources
#' @param diet_prop
#' @param food
#'
#' @return
#' @export
#'
#' @examples
eatfoodonce <- function(num_sources, diet_prop, food) {
  #choose a diet source population
  source <- sample(x = c(num_sources), size = 1, prob = diet_prop)
  #choose an individual from that population
  eaten <- sample(food[[source]], size = 1)
  return(eaten)
}

#' eatfood
#'
#' @param num_individ
#' @param steps
#' @param num_sources
#' @param diet_prop
#' @param food
#'
#' @return
#' @export
#'
#' @examples
eatfood <- function(num_individ, steps, num_sources, diet_prop, food) {

  diethistory <- replicate(num_individ,
                           replicate(n = steps,
                                     eatfoodonce(num_sources = num_sources,
                                                  diet_prop = diet_prop,
                                                  food = food)))
  return(diethistory)
}
