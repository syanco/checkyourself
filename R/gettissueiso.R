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

#convenience fucntion to format "observed" isodata into a dataframe with columns
#for d13C and d15N and each row being an indivdiual observation
#default is to write the object to a .csv file but can also return dataframe as
#an object
formatiso <- function(iso, filename, writefile = T, returnobject = F){
  #error checks
  if(!is.character(filename))
    stop("Filename must be supplied as a character string")

  #format data
  dat <- data.frame(matrix(unlist(iso), ncol = 2, byrow = T))
  names(dat) <- c("d13C", "d15N")

  #write file routine
  if(writefile == T)
    write.csv(dat, filename)

  #return dataframe routine
  if(returnobject == T)
    return(dat)
}
