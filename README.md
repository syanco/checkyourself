# checkyourself
Repository for the SDM simulation model associated with Yanco et al. 2018 (in prep).

## Installation
If you don't already have `devtools` installed and loaded execute the following:

```{r}
install.packages("devtools")
library(devtools)
```
Install the `checkyourself` package:
```{r}
install_github("syanco/checkyourself")
```
## How to Use

**Citation:** Yanco, SW, AL McDevitt, LM Hartley, and MB Wunder. Check yourself: using agent-based models to test the logic of hypotheses. *In prep*.  

Questions about this vignette or the `chekyourself` package should be directed to Scott Yanco at: Scott.Yanco@ucdenver.edu  

### Introduction

This vignette demonstrates how the `checkyourself` package can be used to simulate a simple agent-based species distribution model (SDM) in a hypothesis vetting process (Yanco et al. *in prep*).  Running this vignette with the parameterizations provided will match the workflow described in Yanco et al. Note that the number of parameterizations and model iterations is large and may take a substantial amount of time to complete. Analysis of the data produced by these simulations is not included in the package or this vignette, but may be found as part of the supplementary information accompanying the manuscript. Note also that many of the processes documented here are stochastic in nature so results will not match exactly those described in the manuscript.  

This package contains code to simulate a patchy landscape consisting of two habitat types and to implement three SDM models on that landscape:  
1. **Null Model:** Individuals settle the landscape randomly with no influence of habitat or neighbors.  
1. **Habitat Preference (HP) Model:** Individuals settle the landscape preferring "Habitat A" over "Habitat B".
1. **Conspecific Attraction (CA) Model:** Individuals settle the landscape preferring to settle near already-settled locations.

Many of the functions included in this package are called by the top-level functions described in this vignette and are not intended to be called directly by the user in running the models as constructed.  Therefore, in this document we only describe how to run the models implemented in the Yanco et al. manuscript.  However, we encourage readers interested in modifying, expanding, or evaluating these models to view the code directly.  All functions are documented in the code such that the model logic should be clear.  

### Simulating a Patchy Landscape  

The package allows the user to simulate a landscape of arbitrary size containing two types of habitats.  The habitats are "grown" as patches and the user may set the target number and size of patches.  

The first step is to load the `checkyourself` package and create the starting matrix on which the habitat will be grown:  
```{r}
library(checkyourself)

#Create starting matrix of specified size, all cell values = '1'
mat.side <- 100 #set starting matrix size
x <- matrix(1, mat.side, mat.side) #create the matrix
```

Now we need to set the declarations for the landscape we are targeting as well as some of the initial objects needed by the habitat growing process.

```{r}
#Declarations
size.clusters <- 100 #define target size of each cluster to be grown (in # of cells)
n.clusters <- 50 #define approximate number of clusters of Habitat A to "grow"
count.max <- 200 #set maximum number of iterations throught he while-loop

#Required Initial Objects
n <- mat.side * mat.side #total number of cells in matrix
cells.left <- 1:n #create 'cells.left' object
cells.left[x!=1] <- -1 # Indicates occupancy of cells
i <- 0 #i counts clusters created and should start at 0 always
indices <- c() #create empty vector for indices
ids <- c() #create empty vector for ids
```

The landscape generation runs as a while loop that either meets its target number of patches of maxes out the number of iterations through the while loop. This loop calls the function `expand` which is part of the `checkyourself` package.  

```{r}
while(i < n.clusters && length(cells.left) >= size.clusters && count.max > 0) {
  count.max <- count.max-1 #countdown against max number of loops
  xy <- sample(cells.left[cells.left > 0], 1) #randomly draw an unoccupied cell
  cluster <- expand(x, size.clusters, xy) #run expand function to grow that cluster
  if (!is.na(cluster[1]) && length(cluster)==size.clusters) {
    i <- i+1 #add to cluster count
    ids <- c(ids, rep(i, size.clusters)) #add cluster to id list
    indices <- c(indices, cluster) #add cluster to indices list
    cells.left[indices] <- -1 #remove all cells in the cluster grown from the available list
  }
}

y <- matrix(NA, mat.side, mat.side) #create blank matrix of the same size as `x`.

#Add the cluster ids to the matrix at locations in indices - this adds each cluster id to the cells indicated by the
#vector 'indices' and leaves the rest of the cells as 'NA'
y[indices] <- ids
```

Now that we have "the "Habitat A" grow, we can create a version of this matrix for each habitat preference strength we intend to run as part of the HP Model. First we need to declare the vector of parameterizations we intend to run for the HP Model.  Note that this declaration here is the declaration of the parameterization for that model run - this step is incorporated in the landscape generation process. The coefficients in `A.coef` (should be the factor(s) of selection, e.g., a 3 in `A.coef` means habitat A is 3 times more likely to be selected than habitat B).

```{r}
#Set the relative strength of selection for Habitat A relative to Habitat B -set as a vector of values through which
#the model iterates
A.coef <- c(1.5, 2, 2.5, 3, 3.5)
```

Now we can convert the value of all cells where Habitat A was "grown" to the value of the selection coefficient(s) `A.coef` using the function `pref.strength`. We implement this within a call to `sapply` in order to vectorize the process across the multiple declared values of `A.coef`. 

```{r}
hab.mat <- sapply(A.coef, pref.strength, mat = y)
```

`hab.mat` is now an array containing the values of each of the five habitat matrices, one for each value of `A.coef` supplied. We must also create a version of each of these matrices that converts the factors of selection to relative probabilities normalized by the sum of each of the matrices.  This can be accomplished using the function `convert.cell` which we call inside `apply` in order to vectorize across each element of the `hab.mat` array which contains, in this case, five sets of matrix values. Each of the individual habitat matrices is stored in the second dimension of the `hab.mat` array, hence the 2 in the second argument of `apply`.

```{r}
p.mat <- apply(hab.mat, 2, convert.cell)
```

### Run Simulations

We are now ready to simulate the multiple parameterizations of the three models on this simulated landscape.  We first establish the declarations which parameterize the models and determine the number of individual agents included in each model as well as the number of iterations of each model parameterization to run.  We also create an ID matrix (`IDmat`) whose cell values mach the cell IDs which is used to get cell locations for the CA model. Not also that we have already declared the parameterization for the HP Model as part of the landscape generation process.  We print these parameterizations along with the others again for convenience.

```{r}
reps <- 1000 #number of iterations of the simulation
n.individ <- 100 #number of animals to simulate settling
radius <- c(1,3,5,10,25) #define settlment radius for conspecific attraction
print(paste("Number of individuals per model run: ", n.individ, sep = ""))
print(paste("Number of model iterations per parameterization: ", reps, sep = ""))
print(paste(c("HP Parameters: ", A.coef), sep = ""))
print(paste(c("CA Parameters: ", radius), sep = ""))
print("Null Model contains no free parameters")

#create ID matrix with cell values corresponding to cell ID
IDmat <- matrix(1:dim(y)[1]^2, nrow = dim(y)[1])
```

#### Null Model

The Null Model can be run with a single call to `rep.rand`.  We supply `rep.rand` with several of the parameterizations from above.  We also need to supply the function with one of the habitat matrices in `hab.mat` as well as the corresponding `A.coef`.  While the simulation itself does not utilize this habitat matrix, because part of the output reports the proportion of settled location within Habitat A one of the habitat matrices is require.  It does not matter which one since they are identical in their spatial geometry and the magnitude of the preference encoded is not considered. The only requirement is that the habitat matrix supplied must match the value of `A.coef` supplied, since that value serves as the "key" to which cells are part of Habitat A.  Remember that the habitat matrices are stored in the second dimension of the `hab.mat` array. Below we simply use the first habitat matrix and the first value of `A.coef`.

```{r eval=FALSE}
NULLmod <- repRand(reps=reps, n.individ = n.individ, hab.mat = hab.mat[,1], A.coef = A.coef[1])
```

#### Habitat Preference Model

The HP Model is run by a call to `rep.hab`. This function takes all the same arguments as `rep.rand` plus the probability matrices we calculated from `hab.mat`.  Note that for this model we must supply all the habitat matrices and values of `A.coef` so that the model can iterate through each parameterization.

```{r eval=FALSE}
HABmod <- repHab(p.mat=p.mat, hab.mat = hab.mat, reps = reps, n.individ = n.individ, A.coef = A.coef)
```

#### Conspecific Attraction Model

The CA Model is run with a call to `rep.con`. This model again requires one of the habitat matrices contained in `hab.mat` and the corresponding value from `A.coef` as was the case for the Null Model.  In addition to the `reps` and `n.individ` required by all the models, this function also required the vector of parameterizations for the CA Model contained in `radius`.  Finally, this matrix required that we supply a matrix with cell values corresponding to cell IDs in the argument `ID.mat`.  Here we give the function the `IDmat` we created in our initial declarations.

```{r eval=FALSE}
CAmod <- repCon(reps = reps, radius = radius, ID.mat = IDmat, hab.mat = hab.mat[,1], n.individ = n.individ, A.coef = A.coef[1])
```

#### Save Data

Once these models have finished running we can save the results and other necessary object to an .Rdata file.  For convenience, in addition to the model outputs, we save the parameterizations we used and the matrices on which the model ran (both `p.mat` and `hab.mat`).  This .Rdata file contains all the output needed (and object names) to run the analysis and visualization code included in Supplement 3 of the Yanco et al. manuscript. The code below will generate an .Rdata file in the working directory that is date stamped.

```{r eval=FALSE}
save(p.mat, hab.mat, A.coef, radius, NULLmod, HABmod, CAmod, file = paste0("sdm_sim_", Sys.Date(), ".Rdata"))
```

[![DOI](https://zenodo.org/badge/152293679.svg)](https://zenodo.org/badge/latestdoi/152293679)
