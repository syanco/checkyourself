#' getrangeSDM
#'
#'Function to extract distribution range from simulated SDM output (i.e., the
#'list output by `repRand`, `repHab`, and/or `repCon`. Expects data of the SDM
#'simulation format and reports the minimum and maximum observed preference for
#'"Habitat A".
#'
#' @param dat Dataset from which to extract range
#' @param ...
#'
#' @return Returns a list containing the minimum and maximum values for the
#' input distribution.
#' @export
getrangeSDM <- function(dat, ...) {
  reduced <- dat[c(FALSE, TRUE)] #exlude raw location data
  min <- min(unlist(reduced))
  max <- max(unlist(reduced))
  return(list(min, max))
}

#' getlow
#'
#'Convenience function to extract low values from `getrangeSDM` output.
#'
#' @param range object output by getrangeSDM
#'
#' @return a single value for the low range estimate for a component of the
#' output of getrangeSDM
#' @export
getlow <- function(range){
  return(range[1])
}

#' gethigh
#'
#'Convenience function to extract high values from `getrangeSDM` output.
#'
#' @param range object output by getrangeSDM
#'
#' @returna single value for the high range estimate for a component of the
#' output of getrangeSDM
#' @export
gethigh <- function(range){
  return(range[2])
}

#' compilerangesSDM
#'
#' @param nulldata data object output by `repRand`.
#' @param habdata data object output by `repHab`.
#' @param condatat data object output by `repCon`.
#' @param A.coef numeric, set of habitat preference coefficients supplied to
#' `repHab`.
#' @param radius numeric, set of conspecific attraction radii supplied to
#' `repCon`.
#'
#' @return a dataframe containg the sampling distribution ranges for all
#' simulation models supplied.
#' @export
compilerangesSDM <- function(nulldata, habdata, condata, A.coef, radius){
  #get null range
  nullRange <- lapply(nulldata, getrangeSDM)

  #get hab pref range
  habRange <- lapply(habdata, getrangeSDM)

  #get conspecific attraction range
  conRange <- lapply(condata, getrangeSDM)

  rangedata <- data.frame("model" = factor(c("Null",
                                              rep("HP Model", length(A.coef)),
                                              rep("CA Model", length(radius)))),
                           "parameter" = c("--", A.coef, radius),
                           "low" = c(unlist(sapply(nullRange, FUN = getlow)),
                                     unlist(sapply(habRange, FUN = getlow)),
                                     unlist(sapply(conRange, FUN = getlow))),
                           "upp" = c(unlist(sapply(nullRange, FUN = gethigh)),
                                     unlist(sapply(habRange, FUN = gethigh)),
                                     unlist(sapply(conRange, FUN = gethigh))))

  return(rangedata)
}

#' singleheatmat
#'
#' Function to calculate matrix of overlap between one distribution and the
#' range of another. Builds the matrix upon which a heatmap can act.
#'
#' @param data1 The first data set to compare. Anticipates the output of one
#' of the `rep...` calls (e.g., `repRand`).
#' @param data2 The second data set to compare. Anticipates the output of one
#' of the `rep...` calls (e.g., `repRand`).
#' @param param1 The numeric parameter object originally fed to the simulation
#' model that generated data set 1.
#' @param param2 The numeric parameter object originally fed to the simulation
#' model that generated data set 2.
#'
#' @return
#' @export
singleheatmat <- function(data1, data2){
  blankmat <- matrix(NA, nrow = length(data1), ncol = length(data2))
  # rownames(blankmat) <- as.factor(param1)
  # colnames(blankmat) <- as.factor(param2)
  range2 <- lapply(data2, getrangeSDM)
  for (i in 1:length(data1)) {
    reduce.dat1 <- unlist(data1[[i]][c(FALSE, TRUE)]) #get just the hab prefs
    for (j in 1:length(range2)) {
      blankmat[i,j] <- sum(reduce.dat1 <= range2[[j]][2] &
                             reduce.dat1 >= range2[[j]][1])/
        length(reduce.dat1)
    }
  }
  return(blankmat)
}

allheatmat <- function(modellist = list(NULLmod, HABmod, CAmod)){
  mat <- outer(X = modellist, Y = modellist, FUN = Vectorize(singleheatmat))
  return(mat)
}
