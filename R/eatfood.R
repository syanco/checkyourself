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
  eaten <- c(sample(food[[source]][[1]], size = 1), #pull a carbon value
                sample(food[[source]][[2]], size = 1)) #pull a nitrogen value
  return(eaten)
}

#vaue returned is a mtrix with each step being a column and each element a row
eatfoodsteps <- function(steps, num_sources, diet_prop, food){
  foodtrack <- replicate(n = steps, eatfoodonce(num_sources = num_sources,
                                                diet_prop = diet_prop,
                                                food = food),
                         simplify = F)
  return(foodtrack)
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
  diethistory <- replicate(n = num_individ, eatfoodsteps(steps = steps,
                                                     num_sources = num_sources,
                                                     diet_prop = diet_prop,
                                                     food = food),
                           simplify = F) #force the result to be a list
  return(diethistory)
}
