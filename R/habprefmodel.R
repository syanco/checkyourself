#' prob.sel
#'
#' Simulate an individual settling the landscape based on habitat preference. The strength of habitat preference
#' is supplied as part of the generation of the simulated landscape in the term A.coef included in the function
#' pref.strength. This function only chooses loctions within the matrix based on the supplied parameterization.
#' @param mat The array of probabilities of selection (output from convert.cell). Anticipates that it will be supplied
#' as part of muli-model vecotrization.
#' @param ... Additional arguments as necessary.
#'
#' @return Returns the cell ID for the chosen cell.
#' @export
prob.sel <- function(mat, ...){
  samp <- sample(length(mat), size=1, prob=(mat)) #choose a random cell in the matrix weighted by the probabilites
  return(samp)
}

#' settle.hab
#'
#' Repeates prob.sel for n.individ number of individuals. Completes a single model run.
#' @param n.individ The number of individuals the model should simulate.
#' Repeats prob.sel this number of times for a single model run.
#' @param p.mat The probability matrix passed to prob.sel. Typically the output from convert.cell
#' @param hab.mat The matrix containing the "grown" habitats enumerated by the strength of preference
#' (e.g., the output from pref.strength, not convert.cell).
#' @param A.coef The habitat preference parameter that matches the hab.mat supplied.
#' This is the "key" that allows the function to interpret each cell of the matrix as either "Habitat A" or "Habitat B".
#' @param ... Additional arguments as necessary.
#'
#' @return A list, the first element of which is the set of coordinates for the settled location of each individual and
#' the 2nd element of which is the proportion of those indviduals that settled within "Habitat A".  If the model is iterated,
#' returns a list of such lists.
#' @export
settle.hab <- function(n.individ, p.mat, hab.mat, A.coef, ...){
  loc <- replicate(n.individ, prob.sel(mat = p.mat)) #run probsel n.individ number of times

  #calculate proportion habitat selected
  hab.sel <- hab.mat[loc] #subset habitat matrix to settled locations
  sel.p <- sum(hab.sel[hab.sel == A.coef])/length(hab.sel) #calculate proportion within Habitat A
  result <- list(loc, sel.p)
  return(result)
}

#' repHab
#'
#' Generate sampling distributions from settle.hab vectorized across parameterizations.
#' Allows settle.hab to be iterated reps number of times across multiple parameterizations of A.coef.
#' @param p.mat The probability matrix passed to prob.sel. Typically the output from `convert.cell`.
#' Passed to `settle.hab`
#' @param hab.mat The matrix containing the "grown" habitats enumerated by the strength of preference
#' (e.g., the output from pref.strength, not convert.cell). Passed to `settle.hab`
#' @param reps Number of model runs to simulate (each run conisting of `n.individ` number of individuals settling)
#' @param n.individ The number of individuals to simulate for a single model run. Passed to `settle.hab`.
#' @param A.coef The habitat preference parameter vector .  This is the "key" that
#' allows the function to interpret each cell of the matrix as either "Habitat A" or "Habitat B". Passed to settle.hab.
#' @param ... Additional arguments as necessary.
#'
#' @return Returns a list, each entry of which is the output of `reps` number of
#' iterations of `settle.hab` (if vecotrized, for a a single parameterization).
#' @export
repHab <- function(p.mat, hab.mat, reps, n.individ, A.coef, ...) {
  sampling.hab <- list()
  for (i in 1:length(A.coef)) {
    sampling.hab[[i]] <- replicate(n=reps, expr = settle.hab(n.individ = n.individ, p.mat = p.mat[,i], hab.mat = hab.mat[,i],
                                        A.coef = A.coef[i]))
  }
  return(sampling.hab)
}
