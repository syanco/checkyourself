gettissue <- function(foodvec, window, time){
  releventisos <- foodvec[(time-window:time)]
  tissueval <- mean(releventisos)
  return(tissueval)
}

getiso <- function(window, time, food){
  if(window >= time)
    stop("Integration window must be less than time step of analysis")
  isos <- apply(food, MARGIN = 2, gettissue, window = window, time = time)
  return(isos)
}
