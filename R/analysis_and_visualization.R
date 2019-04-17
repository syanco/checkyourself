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
compilerangesSDM <- function(nulldata, habdata, condatat, A.coef, radius){
  #get null range
  nullRange <- sapply(nulldata, get.range)

  #get hab pref range
  habRange <- sapply(habdata, get.range)

  #get conspecific attraction range
  conRange <- sapply(condata, get.range)

  rangedatat <- data.frame("model" = factor(c("Null",
                                              rep("HP Model", length(A.coef)),
                                              rep("CA Model", length(radius))),
                                            levels = sort(rbind(A.coef, radius))),
                           "parameter" = c(NA, A.coef, radius),
                           "low" = unlist(c(nullRange[1,], habRange[1,],
                                            conRange[1,])),
                           "upp" = unlist(c(nullRange[2,], habRange[2,],
                                            conRange[2,])))

  return(rangedata)
}

singleheat <- function(data1, data2, param1, param2){
  blankmat <- matrix(NA, nrow = length(data1), ncol = length(data2))
  rownames(blankmat) <- as.factor(param1)
  colnames(p.hab.over.con) <- as.factor(param2)
  for (i in 1:length(sampling.mult.hab)) {
    reduce.hab <- unlist(sampling.mult.hab[[i]][c(FALSE, TRUE)])
    for (j in 1:length(range.con[1,])) {
      p.hab.over.con[i,j] <- length(reduce.hab[reduce.hab < range.con[2,j] & reduce.hab > range.con[1,j]]) /length(reduce.hab)
    }
  }
}
