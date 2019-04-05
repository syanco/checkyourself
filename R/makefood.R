#' Make food
#'
#' Simulates the 13C isotope values for an arbitray number of potential diet
#' endpoints. Each endpoint population is simulated as a normally distributed
#' variable
#'
#' @param num_sources integer, number of potential food sources to simulate
#' @param popsize vector whose length = 'num_sources', each entry of which is
#' the size of the diet source population to simulate
#' @param mu_carb vector whose length = 'num_sources', each entry of which is
#' the mean 13C value of the diet source population to simulate
#' @param sd_carb vector whose length = 'num_sources', each entry of which is
#' the standard deviation of the 13C of the diet source population to simulate
#'#' @param mu_nit vector whose length = 'num_sources', each entry of which is
#' the mean 15N value of the diet source population to simulate
#' @param sd_nit vector whose length = 'num_sources', each entry of which is
#' the standard deviation of the 15N of the diet source population to simulate
#'
#' @return A list, each element of which is a vector of isotope vales for each
#' simulated diet endpoint individual
#' @export
#'
#' @examples
makefood <- function(num_sources, popsize, mu_carb, sd_carb, mu_nit, sd_nit) {
  #error checks
  if(length(mu_carb) != num_sources)
    stop("You must provide exactly 1 d13C mean for each proposed source")
  if(length(sd_carb) != num_sources)
    stop("You must provide exactly 1 d13C sd for each proposed source")
  if(length(mu_nit) != num_sources)
    stop("You must provide exactly 1 d15N mean for each proposed source")
  if(length(sd_nit) != num_sources)
    stop("You must provide exactly 1 d15N sd for each proposed source")

  #initialize list to hold simulated food data
  food <- list()
  #loop through each source to generate carban and nitrogen values
  for(i in 1:num_sources){
    food[[i]] <- list(rnorm(n = popsize[i], mean = mu_carb[i], sd = sd_carb[i]),
                      rnorm(n = popsize[i], mean = mu_nit[i], sd = sd_nit[i]))
    }
  return(food)
  }

#' savesources
#'
#' @param food list object of simulated food sources
#' @param filename name of file to save
#' @param mu_carb vector whose length = 'num_sources', each entry of which is
#' the mean 13C value of the diet source population to simulate
#' @param sd_carb vector whose length = 'num_sources', each entry of which is
#' the standard deviation of the 13C of the diet source population to simulate
#'#' @param mu_nit vector whose length = 'num_sources', each entry of which is
#' the mean 15N value of the diet source population to simulate
#' @param sd_nit vector whose length = 'num_sources', each entry of which is
#' the standard deviation of the 15N of the diet source population to simulate
#' @param popsize vector whose length = 'num_sources', each entry of which is
#' the size of the diet source population to simulate
#' @param writefile boolean, whether to write the csv to file needed by MixSIAR
#' default is TRUE
#' @param returnobject boolean, whether to return an object within the current
#' environment, default is FALSE
#'
#' @return either writes a MixSIAR formatted csv to file or a dataframe object
#' @export
#'
#' @examples
savesources <- function(food, filename = NULL, mu_carb, sd_carb, mu_nit,
                        sd_nit, popsize, writefile = T, returnobject = F){


  #format data
  dat <- data.frame("Sources" = 1:length(food), "Meand15N" = mu_nit,
                    "SDd15N" = sd_nit, "Meand13C" = mu_carb,
                    "SDd13C" = sd_carb, "n" = popsize)

  #write file routine
  if(writefile == T)
    #error checks
    if(!is.character(filename))
      stop("Filename must be supplied as a character string")
    write.csv(dat, filename, row.names=FALSE)

  #return dataframe routine
  if(returnobject == T)
    return(dat)
}

#' savediscrimination
#'
#' function to simulate isotope discrimination data in the format needed by
#' MixSIAR
#'
#' @param food list object of simulated food sources
#' @param filename name of file to save
#' @param discr_carb vector of carbon discrimination values,
#' one per diet source.  Default is no discrimintation
#' @param discr_nit vector of nitrogen discrimination values,
#' one per diet source.  Default is no discrimintation
#' @param writefile boolean, whether to write the csv to file needed by MixSIAR
#' default is TRUE
#' @param returnobject boolean, whether to return an object within the current
#' environment, default is FALSE
#'
#' @return either writes a MixSIAR formatted csv to file or a dataframe object
#' @export
#'
#' @examples
savediscrimination <- function(food, filename = NULL, discr_carb = c(0,0,0),
                               discr_nit = c(0,0,0), writefile = T,
                               returnobject = F){


  #format data
  dat <- data.frame("Sources" = 1:length(food), "Meand15N" = discr_nit,
                    "SDd15N" = discr_nit, "Meand13C" = discr_carb,
                    "SDd13C" = discr_carb)

  #write file routine
  if(writefile == T)
    #error checks
    if(!is.character(filename))
      stop("Filename must be supplied as a character string")
    write.csv(dat, filename, row.names=FALSE)

  #return dataframe routine
  if(returnobject == T)
    return(dat)
}
