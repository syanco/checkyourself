eatfoodonce <- function(num_sources, diet_prop, food) {
  source <- sample(x = c(num_sources), size = 1, prob = diet_prop)
  eaten <- sample(food[[source]], size = 1)
}

individualeatfood <- function() {

}
