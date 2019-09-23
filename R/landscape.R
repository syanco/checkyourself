#' expand
#'
#' cellular automata-like function to "grow" habitat patches
#'
#' Function to expand a patch randomly within indicator array x (wherein
#' 1=unoccupied) by n.size cells, beginning at index start.
#' Algorithm based on and modified from:
#' https://gis.stackexchange.com/questions/60562/creating-randomly-shaped-clumps-of-
#' cells-in-a-raster-from-seeds-of-1-cell-pixel
#'
#' @param x A matrix whose cell values all = 1
#' @param n.size A number representing the target size for each cluster
#' @param start A number between 1 and the number of cells in matrx x indicating the algorithm starting position.
#'
#' @return A a grid of generated numbers of approximately n.size
#' @export
expand <- function(x, n.size, start) {
  if (x[start] != 1) return(NA)
  n.rows <- dim(x)[1]
  n.cols <- dim(x)[2]
  nbrhood <- matrix(c(-1,-1, -1,0, -1,1, 0,-1, 0,1, 1,-1, 1,0, 1,1), nrow=2)

  # Adjoin one more random cell and update `state`, which records
  # (1) the immediately available cells and (2) already occupied cells.
  grow <- function(state) {

    # Find all available neighbors that lie within the extent of `x` and are unoccupied.
    neighbors <- function(i) {
      n <- c((i-1)%%n.rows+1, floor((i-1)/n.rows+1)) + nbrhood
      n <- n[, n[1,] >= 1 & n[2,] >= 1 & n[1,] <= n.rows & n[2,] <= n.cols,
             drop=FALSE]             # Remain inside the extent of `x`.
      n <- n[1,] + (n[2,]-1)*n.rows  # Convert to *vector* indexes into `x`.
      n <- n[x[n]==1]                # Stick to valid cells in `x`.
      n <- setdiff(n, state$occupied)# Remove any occupied cells.
      return (n)
    }

    # Select one available cell uniformly at random. Return an updated state.
    j <- ceiling(stats::runif(1) * length(state$available))
    i <- state$available[j]
    return(list(index=i,
                available = union(state$available[-j], neighbors(i)),
                occupied = c(state$occupied, i)))
  }

  # Initialize the state. (If `start` is missing, choose a value at random.)
  if(missing(start)) {
    indexes <- 1:(n.rows * n.cols)
    indexes <- indexes[x[indexes]==1]
    start <- sample(indexes, 1)
  }
  if(length(start)==2) start <- start[1] + (start[2]-1)*n.rows
  state <- list(available=start, occupied=c())

  # Grow for as long as possible and as long as needed.
  i <- 1
  indices <- c(NA, n.size)
  while(length(state$available) > 0 && i <= n.size) {
    state <- grow(state)
    indices[i] <- state$index
    i <- i+1
  }

  # Return a grid of generation numbers from 1, 2, ... through n.size.
  indices <- indices[!is.na(indices)]
  y <- matrix(NA, n.rows, n.cols)
  y[indices] <- 2
  return(indices)
}

#' create.land
#'
#' Creates a simulated patchy habitat matrix based on user defined inputs.
#' Uses the expand function to grow each habitat patch.
#'
#' @param n.clusters A number. Target number of clusters to grow on the landscape
#' @param size.clusters A number. Target size (measured in number cells) of each patch.
#' @param x.mat A matrix whose cell values all = 1
#' @param count.max A number. Maximum number of iterations through the growing process (default is 200).
#' @param ... Additional arguments as necessary.
#'
#' @return A matrix with each cell that is a member of a patch identified by it's patch ID. Cells not inculded in a patch are listed as NA
#' @export
create.land <- function (n.clusters, size.clusters, x.mat, count.max = 200, ...) {
  n <- nrow(x.mat) * ncol(x.mat) #total number of cells in matrix
  cells.left <- 1:n #create 'cells.left' object
  cells.left[x.mat != 1] <- -1 # Indicates occupancy of cells
  i <- 0 #i counts clusters created and should start at 0 always
  indices <- c() #create empty vector for indices
  ids <- c() #create empty vector for ids
  while(i < n.clusters && length(cells.left) >= size.clusters && count.max > 0) {
    count.max <- count.max-1 #countdown against max number of loops
    xy <- sample(cells.left[cells.left > 0], 1) #randomly draw an unoccupied cell
    cluster <- expand(x.mat, size.clusters, xy) #run expand function to grow that cluster
    if (!is.na(cluster[1]) && length(cluster)==size.clusters) {
      i <- i+1 #add to cluster count
      ids <- c(ids, rep(i, size.clusters)) #add cluster to id list
      indices <- c(indices, cluster) #add cluster to indices list
      cells.left[indices] <- -1 #remove all cells in the cluster grown from the available list
    }
  }
  y <- matrix(NA, nrow(x.mat), ncol(x.mat)) #create blank matrix

  #Add each cluster id to the cells indicated by the
  #vector 'indices' and leaves the rest of the cells as 'NA'
  y[indices] <- ids
  return(y)
}

#' pref.strength
#'
#' Converts the values of the simulated patchy habitat matirx to strength of preference values supplied by the number or vector A.coef.
#' If A.coef is supplied as vector is vectorizes the process to create a matrix for each value of A.coef
#'
#' @param A.coef Numeric. Can be a number of vector.  Describes the strength(s) of preference for Habitat A over Habitat B.
#' @param mat Matrix.  Tne y matrix produced by create.land.  Converts the Habitat patches of this matrix
#' (which are labelled) by patch ID to the value(s) of A.coef.
#' @param ... Additional arguments as necessary.
#'
#' @return A matrix with values of the grown habitat converted to relative strengths of preference OR
#' an array of matrix values if A.coef is supplied as a vector
#' @export
pref.strength <- function (A.coef, mat, ...) {
  mat[which(mat != is.na(mat))] <- A.coef #convert all non-NA cells (where patches were grown) to the value of A.coef
  #set Habitat B to coef = 1 - these are all the spaces not filled by the patches "grown" above
  mat[which(is.na(mat))] <- 1
  return(mat)
}

#' convert.cell
#'
#' Converts the values of each cell in a matrix
#' to probabilities by normalizing each cell's value by the sum of the entire matrix.
#' @param mat Matrix to be converted (typically used to convert the output of pref.strength)
#' @param ... Additional arguments as necessary.
#'
#' @return A matrix with the values of each cell converted to relative probabilities,
#'  normalized byt the sum of the values of the matrix
#' @export
convert.cell <- function (mat,  ...) {
  converted <- mat/sum(mat)
  return(converted)
}

#' startmat
#'
#'Creates a square matirx of supplied size with all cell values being 1.  This
#'is the matrix anticipated by the habitat growth functions.
#'
#' @param matsize integer, the size of a side of the suare matrix to be created.
#'
#' @return a square matrix whose sides are of size `matsize` and all cell values
#'  are "1".
#'
#' @export
startmat <- function(matsize){
  x <- matrix(1, matsize, matsize)
  return(x)
}

#' simpatches
#'
#'Function to simulate a patchy habitat matrix by "growing" one habitat type
#'from a cellular automota-like loop.
#'
#' @param matsize integer, the size of a side of the suare matrix to be created.
#' @param count.max interger, the maximum number of iterations throught the
#' while-loop.  Default value is 200
#' @param n.clusters integer, the number of indpendent habitat patches of
#' "type A" to grow.
#' @param size.clusters integer, the target average size (in number of matrix
#' cells) for each patch.
#'
#' @return Returns a matrix of size `matsize` x `matsize` containing
#' approximately `n.clusters` patches of approximately `size.clusters` size.
#' @export
simpatches <- function(matsize, count.max = 200, n.clusters, size.clusters){
  x <- startmat(matsize = matsize) #run startmat to get starting matrix
  n <- matsize^2 #calc total number of cells
  cells.left <- 1:n #create 'cells.left' object
  cells.left[x!=1] <- -1 # Indicates occupancy of cells
  i <- 0 #i counts clusters created and should start at 0 always
  indices <- c() #create empty vector for indices
  ids <- c() #create empty vector for ids
  max <- count.max

  while(i < n.clusters &&
        length(cells.left) >= size.clusters &&
        max > 0) {
    max <- max-1 #countdown against max number of loops
    xy <- sample(cells.left[cells.left > 0], 1) #randomly draw an unoccupied cell
    #run expand function to grow that cluster
    cluster <- expand(x, size.clusters, xy) #
    if (!is.na(cluster[1]) && length(cluster)==size.clusters) {
      i <- i+1 #add to cluster count
      ids <- c(ids, rep(i, size.clusters)) #add cluster to id list
      indices <- c(indices, cluster) #add cluster to indices list
      #remove all cells in the cluster grown from the available list
      cells.left[indices] <- -1
    }
  }
  y <- matrix(NA, matsize, matsize) #create blank matrix of the same size as `x`.

  #Add the cluster ids to the matrix at locations in indices - this adds each
  #cluster id to the cells indicated by the
  #vector 'indices' and leaves the rest of the cells as 'NA'
  y[indices] <- ids

  return(y)
}

#' simlandscapes
#'
#' @param A.coef Numeric. Can be a number of vector.  Describes the strength(s)
#' of preference for Habitat A over Habitat B.
#' @param matsize integer, the size of a side of the suare matrix to be created.
#' @param count.max interger, the maximum number of iterations throught the
#' while-loop.  Default value is 200
#' @param n.clusters integer, the number of indpendent habitat patches of
#' "type A" to grow.
#' @param size.clusters integer, the target average size (in number of matrix
#' cells) for each patch.
#'
#' @return A list each element of which contains the "habitat matrix" and
#' "preference strength" matrix for the given value of A.coef. There is an
#' element in the list for each value of A.coef supplied.
#' @export
simlandscapes <- function(A.coef, matsize, count.max = 200, n.clusters,
                         size.clusters, writeraster = F, filename = "1"){
  y <- simpatches(matsize=matsize, n.clusters = n.clusters,
                  size.clusters = size.clusters)
  hab.mat <- sapply(A.coef, pref.strength, mat = y)
  p.mat <- apply(hab.mat, 2, convert.cell)
  IDmat <- matrix(1:matsize^2, nrow = matsize)
  # if(writeraster = T) {
  #   writeRaster(hab.mat, filename = paste0(filename,".tif"))
  # }
  return(list(hab.mat, p.mat, IDmat))
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
