#overall location selection function
choose.loc <- function(r, pd.rate, steps, lambda, coef.d, coef.r, blank.rast, ...) {
  move.list <- list(c(0,0)) #establish list to store locations
  for (i in 2:steps) {
    if (move.list[[i-1]][1] == 0 & move.list[[i-1]][2] == 0)  {
      #check last location, if at nest - allow hab-based movement
      list.tmp <- list() #create temp list to store each location per iteration
      loc <- probsel(r, move.list = move.list, coef.d = coef.d, coef.r = coef.r, blank.rast = blank.rast)
      list.tmp[1] <- loc@coords[1]
      list.tmp[2] <- loc@coords[2]
      move.list[[i]] <- list.tmp
    } else { #if not at nest, return to nest probabilistically (set by pd.rate)
      if (runif(1) <= pd.rate) {
        list.tmp <- list()
        list.tmp[1] <- 0
        list.tmp[2] <- 0
        move.list[[i]] <- list.tmp
      } else {
        list.tmp <- list() #create temp list to store each location per iteration
        loc <- probsel(r, move.list = move.list, coef.d = coef.d, coef.r = coef.r, blank.rast = blank.rast)
        list.tmp[1] <- loc@coords[1]
        list.tmp[2] <- loc@coords[2]
        move.list[[i]] <- list.tmp
      }
    }
  }
  hr.kde.90 <- hrFun(move.list)
  return(hr.kde.90)
}
