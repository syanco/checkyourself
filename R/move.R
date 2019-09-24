#overall location selection function
choose.loc <- function(r, pd.rate, steps, lambda, coef.d, coef.r, blank.rast, matsize, ...) {
  move.list <- list(c(0,0)) #establish list to store locations
  for (i in 2:steps) {
    #check last location, if at nest - allow hab-based movement
    if (move.list[[i-1]][1] == 0 & move.list[[i-1]][2] == 0)  {
      list.tmp <- c() #create temp list to store each location per iteration
      loc <- probsel(r, move.list = move.list, coef.d = coef.d, coef.r = coef.r,
                     blank.rast = blank.rast, matsize = matsize,
                     lambda = lambda)
      #pull the coords out so that we can store as x, y not row, col
      list.tmp[1] <- loc[2]
      list.tmp[2] <- loc[1]
      #update main movelist
      move.list[[i]] <- list.tmp
    } else { #if not at nest, return to nest probabilistically (set by pd.rate)
      if (runif(1) <= pd.rate) {
        list.tmp <- list()
        list.tmp[1] <- 0
        list.tmp[2] <- 0
        move.list[[i]] <- list.tmp
      } else {
        list.tmp <- list() #create temp list to store each location per iteration
        loc <- probsel(r, move.list = move.list, coef.d = coef.d, coef.r = coef.r,
                       blank.rast = blank.rast, matsize = matsize,
                       lambda = lambda)
        #pull the coords out so that we can store as x, y not row, col
        list.tmp[1] <- loc[2]
        list.tmp[2] <- loc[1]
        #update main movelist
        move.list[[i]] <- list.tmp
      }
    }
  }
  return(move.list)
}

#(based on https://scrogster.wordpress.com/2012/04/22/using-r-for-spatial-sampling-with-selection-probabilities-defined-in-a-raster/)
probsel <- function(r, lambda, move.list, coef.d, coef.r, blank.rast, matsize, ...){
  #get distance probabilities
  d <- makeDistProb(matsize = matsize, position = 2, lambda = lambda)
  #combine distance probs with habitat probs, weighted by mixing coefficients
  x <- (coef.d * d) + (coef.r * r)
  x <- x/sum(x) #convert back to true probabilities
  samp <- sample(matsize^2, size=1, prob=x)
  idmat <- matrix(1:matsize^2, nrow = matsize, byrow = F) #create an index matirx
  coords <- which(idmat == samp, arr.ind = T)
}

# #select the distance raster based on owl's position
# dist.prob.alt <- function (move.list, blank.rast, ...) {
#   pos <- cellFromXY(blank.rast, as.numeric(tail(move.list, 1)[[1]]))
#   d.rast <- raster(paste0("prob_rast/p", pos,".tif"))
#   return(d.rast)
# }

#' makeDistProb
#'
#'Function to make a matrix of size `matsize` containing probabilities that
#'decline exponentially by distance from a given cell (`position`). Strength of
#'exponential decline is given by `lambda`
#' @param matsize numeric, size of one side of a square matrix on which to
#' perform calculations
#' @param position numeric, position from which to generate declining
#' probability
#' @param lambda numeric, strength of exponential decline ($\lambda$) where:
#' $x = \lambda^{-\lambda x}$
#'
#' @return
#' @export
#'
#' @examples
makeDistProb <- function (matsize, position, lambda) {
  #create new empty matrix
  d.mat <- matrix(1:matsize^2, nrow = matsize, byrow = F) #create an index matirx
  x <- sapply(d.mat, matrixPythagoras, IDmat = d.mat, position2 = position)
  decreasebydist <- lambda*exp((-lambda*x))
  probmat <- decreasebydist/sum(decreasebydist)
}

#' matrixPythagoras
#'
#' @param IDmat ID matrix on which to calculate distance, must have cell values
#' matching 1-D indexing
#' @param position1 first positions in distance calc
#' @param position2 second position in dtance calc
#'
#' @return returns the pythogorean distance between two coordinates (supplied as
#'matirx indices)
#' @export
#'
#' @examples
matrixPythagoras <- function(position1, position2, IDmat){
  x <- as.numeric(sqrt(((which(d.mat == position1, arr.ind = T)[,"row"]-
                           which(d.mat == position2, arr.ind = T)[,"row"])^2) +
                         (which(d.mat == position1, arr.ind = T)[,"col"]-
                            which(d.mat == position2, arr.ind = T)[,"col"])^2))
}
