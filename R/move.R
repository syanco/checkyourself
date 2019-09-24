#' chooseLoc
#'
#'
#' A function to simulate a single individual moving with varying mixtures of
#' preference for habitat and/or short steps.
#'
#' @param hab.prob matrix, the habitat matrix with cell values representing the
#' habitat-only probabilities (uncorrected for the mixing proportions with
#' distance). Typically the second element of the list returned by
#' simlandscapes.
#' @param pd.rate numeric, the probability of returning to the central place to
#' "deliver prey" to nest/den/etc. on any given step (given the agent isnot
#' currently located at the central place.) Must be <1.
#' @param steps integer, the number of steps the agent should take
#' @param lambda numeric, the strength of the declining exponental function that controls
#' how "reluctant" the agent is to move large distances. Must be <1.
#' @param coef.d numeric, can be factor or proportion representing the relative
#' strength of the preference to move short steps compared to the strength of
#' the preference for habitat.
#' @param coef.r numeric, can be factor or proportion representing the relative
#' strength of the preference for habitat compared to the strength of
#' the preference for short steps.
#' @param matsize integer, the size in number of cells of one side of the square
#' matrix on which the simulation is being run.
#' @param ...
#'
#' @return Returns the list of positions for the simulated individual.
#' @export
#'
#' @examples
chooseLoc <- function(hab.prob, pd.rate, steps, lambda, coef.d, coef.r, blank.rast, matsize, ...) {
  move.list <- list(c()) #establish list to store locations
  for (i in 2:steps) {
    #check last location, if at nest - allow hab-based movement
    if (move.list[[i-1]][1] == 0 & move.list[[i-1]][2] == 0)  {
      list.tmp <- c() #create temp list to store each location per iteration
      loc <- probsel(hab.prob = hab.prob, move.list = move.list,
                     coef.d = coef.d, coef.r = coef.r, blank.rast = blank.rast,
                     matsize = matsize, lambda = lambda)
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
        loc <- probsel(hab.prob, move.list = move.list, coef.d = coef.d, coef.r = coef.r,
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


#' probsel
#'
#' Makes probabalistic selections of a cell based on the mixture of two
#' probability matrices. Not intended to be called directly - used by
#' `chooseloc`
#'
#' (based on https://scrogster.wordpress.com/2012/04/22/using-r-for-spatial-sampling-with-selection-probabilities-defined-in-a-raster/)
#'
#' @param hab.prob matrix, the habitat matrix with cell values representing the
#' habitat-only probabilities (uncorrected for the mixing proportions with
#' distance). Typically the second element of the list returned by
#' simlandscapes.
#' @param lambda numeric, the strength of the declining exponental function that controls
#' how "reluctant" the agent is to move large distances. Must be <1.
#' @param move.list list, the current `move.list` being maintained by
#' `chooseloc`. Used to find the current agent position.
#' @param coef.d numeric, the relative weight to give the distance preference.
#' @param coef.r numeric, the relative weight to give the habitat preference.
#' @param matsize integer, the size in number of cells of one side of the square
#' matrix on which the simulation is being run.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
probsel <- function(hab.prob, lambda, move.list, coef.d, coef.r, matsize, ...){
  #get distance probabilities
  d <- makeDistProb(matsize = matsize,
                    position = move.list[[length(move.list)]],
                    lambda = lambda)
  #combine distance probs with habitat probs, weighted by mixing coefficients
  x <- (coef.d * d) + (coef.r * hab.prob)
  x <- x/sum(x) #convert back to true probabilities
  samp <- sample(matsize^2, size=1, prob=x)
  idmat <- matrix(1:matsize^2, nrow = matsize, byrow = F) #create an index matirx
  coords <- which(idmat == samp, arr.ind = T)
  return(coords)
}

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
  decreasebydist <- lambda*exp(((-lambda)*x))
  probmat <- decreasebydist/sum(decreasebydist)
}

#' matrixPythagoras
#'
#' Uses the Pythagorean theorem to calculate the strightline distance between
#' two matrix cells.
#'
#' @param IDmat matrix, the ID matrix on which to calculate distance, must have
#' cell values matching 1-D indexing. See examples
#' @param position1 first positions in distance calc
#' @param position2 second position in disttance calc
#'
#' @return returns the pythogorean distance between two coordinates (supplied as
#'matirx indices)
#' @export
#'
#' @examples
#' #create an example ID matrix for `IDmat`. Example uses a 100x100 matrix
#' mat <- matrix(1:10000, nrow = 100)
#'
matrixPythagoras <- function(position1, position2, IDmat){
  x <- as.numeric(sqrt(((which(IDmat == position1, arr.ind = T)[,"row"]-
                           which(IDmat == position2, arr.ind = T)[,"row"])^2) +
                         (which(IDmat == position1, arr.ind = T)[,"col"]-
                            which(IDmat == position2, arr.ind = T)[,"col"])^2))
}


#' matrixCellFromXY
#'
#' Gets the cell ID within a matrix using the position format (`c(X, Y)`)
#' suppplied by `move.list`.
#'
#' @param pos vector, a position in the format `c(X, Y)` i.e. `c(col, row)`.
#' Intended to accept a single position in `move.list`
#' @param matsize integer, the size of the square matrix on which the simulation
#' is being run
#'
#' @return integer, the cell ID indicated by the X, Y coordinates.
#' @export
#'
#' @examples
matrixCellFromXY <- function(pos, matsize) {
  cell <- ((pos[1]-1)*matsize) + pos[2]
  return(cell)
}

