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
