source("R/makefood.R")
source("R/eatfood.R")
source("R/gettissueiso.R")

#declarations
num_sources <- 3 #number of diet endpoints
popsize <- c(1000,1000,1000) #size of each pop to simulate
mu_carb <- c(-1, -10, -30) #vector of mean carbon values for each source
sd_carb <- c(1,1,1) #vector of carbon source sds
mu_nit <- c(10,6,11) #vector of carbon discrimination factors (set to 0)
sd_nit <-c(1,1,1) #vector of nitrogen discrimination factors (set to 0)
diet_prop <- c(.1, .25, .65) #"real" proportional diet composition
steps <- 100 #number of steps per model run
num_individ <- 100 #number of consumers to simulate

#make the food sources
food <- makefood(3, popsize, mu_carb, sd_carb, mu_nit, sd_nit)
#save diet source data to file for import to MixSIAR
savesources(food, filename = "simulated_sources.csv", mu_carb, sd_carb, mu_nit,
            sd_nit, popsize)

#create and save discrimination factor data
#(NOTE: model only accomodates TDR = 0)
savediscrimination(food, filename="simulated_discrimination.csv")

#simulate agents over time
specimens <- eatfood(num_individ, steps, num_sources, diet_prop, food)

#collect "observed" isotop values
obs_iso <- getiso(specimens, time = steps, window = 10)

#format to data frame ans save as csv for import to MixSIAR
formatiso(obs_iso, filename = "simulated_iso.csv")

#sample analysis
library(MixSIAR)
library(splancs)

#load consumer/mixture data
mix <- load_mix_data(filename="simulated_iso.csv",
                     iso_names=c("d13C","d15N"),
                     factors=NULL,
                     fac_random=NULL,
                     fac_nested=NULL,
                     cont_effects=NULL)

#load diet source data
source <- load_source_data(filename="simulated_sources.csv",
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix)
#load discrimination data
discr <- load_discr_data(filename="simulated_discrimination.csv", mix)

#make plot of the isoscape
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE,
          mix, source, discr)

#calculate convex hull of isoscape
calc_area(source=source, mix=mix, discr=discr)

#plot priors (using generalist/uninformative prior)
plot_prior(alpha.prior=1, source)

# Write the JAGS model file
model_filename <- "MixSIAR_sim_test_model_040419.txt"
resid_err <- FALSE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

#run the model (using a 'test' lengthed chain)
jags.1 <- run_model(run="very short", mix, source, discr, model_filename,
                    alpha.prior = 1, resid_err, process_err)

#if test works, run analysis length
#jags.1 <- run_model(run="normal", mix, source, discr, model_filename,
#                    alpha.prior = 1, resid_err, process_err)

#set JAGS output options
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

#view the output
output_JAGS(jags.1, mix, source, output_options)
