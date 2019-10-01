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
#' @return Returns a matirx or vector (depending on number of paraeters) of
#'  proportion of `data1` contained in `data2`
#' @export
singleheatmat <- function(data1, data2){
  blankmat <- matrix(NA, nrow = length(data1), ncol = length(data2))
  rownames(blankmat) <- as.factor(names(data1))
  colnames(blankmat) <- as.factor(names(data2))
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

#' getcombos
#'
#'Creates 2 lists which, paired, provide each unidirectional pairwise
#'combination (outer product) of the two lists. Meant as a component function to
#'assist in creating heatmaps.
#'
#' @param X The first list of model objects. It is suggested that each list
#' element be named.
#' @param Y The second list.
#'
#' @return A two element list which contains the paired combinations of the
#' supplied model lists
#'
#' @export
getcombos <- function(X, Y){
  Y <- rep(Y, rep.int(length(X), length(Y)))
  if (length(X))
    X <- rep(X, times = ceiling(length(Y)/length(X)))
  names(X) <- paste0(names(X),"|",names(Y)) #set names of X to the combo
  return(list(X,Y))
}

#' vectorheat
#'
#'The vectorized wrapper function for `singleheatmat` that can walk through the
#'paired list of model combiantions produced by `getcombos` and produce matrices
#'of overlap data with which heat maps can be generated.
#'
#' @param mod1 A list of model objects, each element of which will be the first
#' data set to compare in the `singleheatmat` function. This function is
#' intended to recieve the first element of the list produced by `getcombos`.
#' @param mod2 The second list of model objects, each element of which will be
#' the second data set to compare in the `singleheatmat` function. This
#' function is intended to recieve the second element of the list produced by
#'  `getcombos`.
#'
#' @return Returns a list, each element of which is matrix of the sampling
#' distribution overlaps as produced by `singleheatmat`.
#' @export
#'
#' @examples modelcomb <- getcombos(X = list("NULL"=NULLmod, "HP"=HABmod,
#' "CA"=CAmod), Y = list("NULL"=NULLmod, "HP"=HABmod, "CA"=CAmod))
#'  heats <- vectorheat(modelcomb[[1]], modelcomb[[2]])
vectorheat <- function(mod1, mod2){
  heats <- mapply(FUN = singleheatmat, mod1, mod2, SIMPLIFY = F, USE.NAMES = T)
  return(heats)
}

#' meltheatmats
#'
#'Convenience function to `melt` all the heatmap matirces into a long-form
#'data.frame format that can be used by `ggplot2`. Retains list naming so that
#'variables are correctly identified in created data.frames.
#'
#' @param mats A list of overlap matrices as is produced by `vectorheat`.
#'
#' @return Rturns a list of dataframes in long form, suitable for heat map
#' plotting via `ggplot2`.
#'
#' @export
meltheatmats <- function(mats){
  meltedmats <- list()
  for(i in 1:length(mats)){
    meltedmats[[i]] <- reshape2::melt(data = mats[[i]],
                                      varnames = c(
                                        sub('\\|.*','',
                                            names(mats[i])),
                                        sub('.*\\|','',
                                            names(mats[i]))),
                                      value.name = "overlap", as.is = T)
    meltedmats[[i]][,1] <- factor(meltedmats[[i]][,1])
    meltedmats[[i]][,2] <- factor(meltedmats[[i]][,2])
  }
  return(meltedmats)
}

#' plotsingleheat
#'
#'Function to creat a `ggplot2` heatmap from a single  "melted" matrix of
#'sampling distribution overlap. Intended to be called by vectorized wrapper function, but could be called on a single element of the list produced by
#'`meltheatmats` to produce a single heatmap
#'
#' @param data A single long-form dataframe object containing the
#' parameterizations of the models to be compared and the calculated pairwise
#'  unidirectional overlap in distributions.
#'
#' @return
#' @export
#'
#' @examples
plotsingleheat <- function(data){
  if(names(data)[1] != names(data)[2]){
    if(length(levels(data[,1])) < length(levels(data[,2]))){
      ggplot(data, aes_string(x = names(data)[1],
                              y = names(data)[2])) +
        geom_tile(aes_string(fill = names(data)[3]), color="black", size = 0.25) +
        coord_equal() +
        geom_text(aes_string(label = names(data)[3]), size = 2) +
        scale_fill_gradient2(low="white", high="olivedrab", limits=c(0,1),
                             guide = F) +
        xlab(as.character(names(data)[1])) +
        ylab(as.character(names(data)[2])) +
        ggtitle(paste0("p(", as.character(names(data)[1]), "|",
                       as.character(names(data)[2]), ")")) +
        theme(#axis.text.x = element_text(color="white"),
              axis.title = element_text(size=8),
              axis.ticks = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0.2, l = 0, unit = "in"),
              plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
              panel.border = element_blank(),
              panel.background = element_blank())
    } else {
      ggplot(data, aes_string(x = names(data)[2],
                              y = names(data)[1])) +
        geom_tile(aes_string(fill = names(data)[3]), color="black", size = 0.25) +
        coord_equal() +
        geom_text(aes_string(label = names(data)[3]), size = 2) +
        scale_fill_gradient2(low="white", high="olivedrab", limits=c(0,1),
                             guide = F) +
        xlab(as.character(names(data)[2])) +
        ylab(as.character(names(data)[1])) +
        ggtitle(paste0("p(", as.character(names(data)[1]), "|",
                       as.character(names(data)[2]), ")")) +
        theme(#axis.text.x = element_text(color="white"),
              axis.title = element_text(size=8),
              axis.ticks = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0.2, l = 0, unit = "in"),
              plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
              panel.border = element_blank(),
              panel.background = element_blank())
    }
  }
}


#' plotheatvector
#'
#' @param meltlist A list ofdataframe objects of the format required by
#' `plotsingleheat`.
#'
#' @return Returns a list of `ggplot2` plot objects
#' @export
#'
#' @examples
plotheatvector <- function(meltlist){
  maps <- lapply(meltlist, FUN = plotsingleheat)
  return(maps)
}

#' getPropUsed
#'
#' Function to get the proportion of used locations (from movement model) within
#'  the habitat matching the value of `A.coef`
#'
#' @param movelist The list of locations visitied by agent in a single model
#' iteration.
#' @param hab.mat The habitat matrix on which the simulation was run, needs
#' version with 1's and `A.coef` (not probabilities of cell IDs)
#' @param A.coef `A.coef` used to simulate habitat - used as key to `hab.mat`
#' @param matsize The size of one side of the square matrix on which the
#' simualtion is run
#'
#' @return
#' @export
#'
#' @examples
getPropUsed <- function(movelist, hab.mat, A.coef, matsize){
  cells <- sapply(movelist, matrixCellFromXY, matsize = matsize)
  habsused <- hab.mat[cells]
  prop <- sum(habsused == A.coef)/length(habsused)
  return(prop)
}

#' createKDE
#'
#' Ceeates a kernel density estimate based on simulated movement data
#' @param mov movement data supplied in format produced by `chooseLoc`
#'
#' @return
#' @export
#'
#' @examples
createKDE <- function (mov) {
  sp <- sp::SpatialPoints(matrix(unlist(mov), ncol = 2, byrow = T))
  kde <- adehabitatHR::kernelUD(sp, h = "href")
  return(kde)
}

#' getHRUsed
#'
#' Extracts the habitat classification of habitat matrix cells contained within
#' the boundaries of a supplied polygon (intended to be the isopleth of a
#' KDE-based home range)
#' @param poly A polygon of the format `SpatialPolygons` or
#' `SpatialPolygonsDataFrame`. Typically the output of the `getverticesHR`
#' function in the `adehabitatHR` package.
#' @param hab.mat The habitat classification matrix (first element of the list
#' generated by `simlandscapes`)
#' @param matsize The size of one side of the square matrix on which the
#' simualtion is run.
#'
#' @return
#' @export
#'
#' @examples
getHRUsed <-function(poly, hab.mat, matsize) {
  #convert hab.mat to a raster object
  r <- raster::raster(matrix(hab.mat, ncol = matsize), xmn = 0, xmx = matsize,
                      ymn = 0, ymx = matsize)

  #extract the raster cells contained within the boundaries of poly
  used <- raster::extract(r, poly)

  #return the "used" cells
  return(used)
}
get
