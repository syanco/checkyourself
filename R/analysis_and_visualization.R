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

  rangedatat <- data.frame("model" = factor(c("Null",
                                              rep("HP Model", length(A.coef)),
                                              rep("CA Model", length(radius)))),
                           "parameter" = c(NA, A.coef, radius),
                           "low" = c(unlist(sapply(nullRange, getlow)),
                                     unlist(sapply(habRange, getlow)),
                                     unlist(sapply(conRange, getlow))),
                           "upp" = c(unlist(sapply(nullRange, gethigh)),
                                     unlist(sapply(habRange, gethigh)),
                                     unlist(sapply(conRange, gethigh))))

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
