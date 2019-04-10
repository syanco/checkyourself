iterateiso <- function (num_individ, steps, num_sources, diet_prop, food,
                        popsize, mu_carb = NULL, sd_carb = NULL, mu_nit = NULL,
                        sd_nit = NULL){


  }

makeparameterspace <- function(carbon = T, nitrogen = T, sources){
  range1 <- c(-20, -1)
  interval1 <- 1
  one <- as.factor(seq(range1[1], range1[2], by = interval1))

  range2 <- c(-10, -1)
  interval2 <- 1
  two <- seq(range2[1], range2[2], by = interval2)

  range3 <- c(-100, -50)
  interval3 <- 1
  three <- seq(range3[1], range3[2], by = interval3)
}
