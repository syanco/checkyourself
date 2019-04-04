gettissue <- function(foodlist, window, time){
  relevant <- foodlist[((time-window+1):time)] #separate the relevant window
  #pull out the carbon values
  carbon <- lapply(relevant, FUN = function(values){carb <- values[1]})
  #pull out the nitrogen values
  nitrogen <- lapply(relevant, FUN = function(values){nit <- values[2]})
  #get the means
  carbmean <- mean(unlist(carbon))
  nitmean <- mean(unlist(nitrogen))
  #combine means into  vector for return
  isos <- c(carbmean, nitmean)
  return(isos)
}

#iterates gettissue across individuals
getiso <- function(window, time, food){
  if(window >= time)
    stop("Integration window must be less than time step of analysis")
  isos <- lapply(food, gettissue, window = window, time = time)
  return(isos)
}
