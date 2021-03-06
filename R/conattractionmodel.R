#' neighbors
#'
#' Function called by settle.con that identifies all matrix cells within `radius` number of cells of any
#' already occupied cell. Flexible to the number of prevously "settled" cells.  This generates the set of
#' available cells from which subsequently settling agents can choose. Can only work on one occupied cell at a time.
#' @param i The set of currently occupied cells.
#' @param radius sets how many cells away from the occupied cells should be included in the "neighborhood"
#' (e.g., 1 means the 8 adjoining neighbors, 2 gives the 8 adjoining neighbors plus the 16 cells outside
#' those - 24 total, etc.).
#' @param mat A matrix of the same size as the simulated landscape but with the value of each cell corresponding
#' to it's cell ID.
#' @param ... Additional arguments as necessary.
#'
#' @return A numeric vector of the matrix cell IDs which are available for selection.
#' @export
neighbors <- function (i, radius, mat, ...) {
  id <- arrayInd(i, .dim = dim(mat)) #gets row and col for occupied cell
  nhood <- expand.grid(-radius:radius, -radius:radius) #generates index for all neighbors w/i specified size
  nhood <- nhood[-which(nhood[1]==0 & nhood[2]==0),] #remove the 0,0 entry so  animals cannot locate ontop of each other
  nhood <- data.matrix(nhood)
  idx <- t(apply(nhood, 1, function(id, nhood) {id + nhood}, id = id)) #applies the abstract 'nhood' to each occupied cell to get 'mat'-specific cell IDs
  #ensure 'idx' doesn't exceed size of landscape matrix
  idx <- idx[abs(idx[,1]) <= dim(mat) &
               abs(idx[,2]) <= dim(mat) &
               idx[,1] > 0 &
               idx[,2] > 0,]
  return(as.numeric(mat[idx]))
}

#' settle.con
#'
#' Agents settle landscape based on conspecific attraction
#' This function acts iteratively to allow agents to settle the landscape selecting randomly from cells defined as
#' available based on proximity to previously settled cells.  The first agent settles the landscape completely
#' randomly. Instanciation as a function allows the entire ABM to be run an arbitray number of times (see 'reps' below)
#'  in order to generate a sampling distribution of apparent selections generated by the model.
#'
#' @param radius Sets how many cells away from the occupied cells should be included in the "neighborhood".
#' Passed to `neighbors`.
#' @param ID.mat A matrix of the same size as the simulated landscape but with the value of each cell corresponding
#' to it's cell ID. Passed to `neighbors`.
#' @param hab.mat The matrix containing the "grown" habitats enumerated by the strength of preference
#' (e.g., the output from pref.strength, not convert.cell.) Note that which hab.mat to use (i.e., which
#' preference strength to use) is arbitrary for the null model, but the choice must match the `A.coef` supplied (below).
#' @param n.individ The number of individuals included in the ABM.
#' @param A.coef The habitat preference parameter that matches the hab.mat supplied.  This is the "key" that
#' allows the function to interpret each cell of the matrix as either "Habitat A" or "Habitat B".
#' @param ... Additional arguments as necessary.
#'
#' @return The proportion of simulated settled locations located in Habitat A as well as the
#' raw locations for each simulation run
#' @export
settle.con <- function(radius, ID.mat, hab.mat, n.individ, A.coef, ...) {
  avail <-ID.mat
  locs <- c(rep(NA, n.individ))
  locs[1] <- sample(avail, size = 1)
  for (l in 2:n.individ) {
    avail <- unlist(sapply(locs[which(locs != is.na(locs))], neighbors,
                           radius = radius, mat = ID.mat))
    locs[l] <- as.numeric(sample(avail, size = 1))
  }

  #calculate proportion habitat selected
  hab.sel <- hab.mat[locs]
  sel.p <- mean(hab.sel == A.coef)
  result <- list(locs, sel.p)
  return(result)
}


#' repCon
#'
#' @param reps Number of model runs to simulate (each run conisting of `n.individ` number of individuals settling)
#' @param radius Sets how many cells away from the occupied cells should be included in the "neighborhood".
#' Passed to `neighbors`. Passed to `settle.con`.
#' @param ID.mat A matrix of the same size as the simulated landscape but with the value of each cell corresponding
#' to it's cell ID. Passed to `settle.con`.
#' @param hab.mat The matrix containing the "grown" habitats enumerated by the strength of preference
#' (e.g., the output from pref.strength, not convert.cell.) Note that which hab.mat to use (i.e., which
#' preference strength to use) is arbitrary for the CA model, but the choice must match the `A.coef` supplied (below) andmust
#' be entered into this function as only the single matrix and not the vector of matirces used to iterate through the HP model.
#' @param n.individ The number of individuals included in the ABM.
#' @param A.coef The habitat preference parameter that matches the hab.mat supplied.  This is the "key" that
#' allows the function to interpret each cell of the matrix as either "Habitat A" or "Habitat B".
#' @param ... Additional arguments as necessary.
#'
#' @return Returns a list, each entry of which is the output of `reps` number of
#' iterations of `settle.hab` (if vecotrized, for a a single parameterization).
#' @export
repCon <- function(reps, radius, ID.mat, hab.mat, n.individ, A.coef, ...) {
  sampling.con <- list()
  for (i in 1:length(radius)) {
    sampling.con[[i]] <- replicate(n=reps, settle.con(radius = radius[i], ID.mat = ID.mat, hab.mat = hab.mat,
                                                      n.individ = n.individ, A.coef = A.coef))
  }
  names(sampling.con) <- as.character(radius)
  return(sampling.con)
}
