library(rstan)
library(loo) # development version w/ stacking
library(dplyr)
library(magrittr)
library(purrr)

source("R/results_processing.R")

# setup ------------------------------------------------------------------------
cmd_arg <- commandArgs(trailingOnly = T)
sex <- cmd_arg[1]
config_file <- cmd_arg[2]
config <- yaml::read_yaml(config_file)

start_year <- config$start_year
end_year <- config$end_year
n_forecast_years <- config$n_forecast_years
model_stem <- config$model

path <- file.path(config$out_dir, config$type, tolower(sex))

cutpoints <- 80:95

## read in log-likelihoods -----------------------------------------------------
lls <- map(cutpoints, get_log_lik, model_stem, path)
lls %<>% map(function(ll) ll[,,ll[1,1,] != 0])

## Calc looic and weights ------------------------------------------------------
options(loo.cores = 4)
loos <- imap(lls, function(ll,i) {
    print(paste("loo", i))
    return(loo(ll, r_eff=relative_eff(ll)))
  }
)
weights <- loo_model_weights(loos)

saveRDS(loos, paste0(path, "/all_loos.rds"))
saveRDS(weights, paste0(path, "/model_weights.rds"))

rm(lls)

## Do stacking -----------------------------------------------------------------
n_total_samples <- 5000
sampled_weights <- rmultinom(1, n_total_samples, weights)
log_rates <- list()
dispersion <- list()
for (i in 1:length(cutpoints)){
  n_samples <- sampled_weights[i]
  if (n_samples > 0){
    cp <- cutpoints[i]
    model <- readRDS(file.path(path, paste0(model_stem, "_", cp, ".rds" )))
    fore <-  get_forecast_log_rates_for_stacking(model$mort_fit, 
                                                 model$stan_data,
                                                 sex)
    disp <- as.matrix(model$mort_fit, "dispersion")
    # choose a random ordering..
    sample_index <- sample.int(dim(fore)[[2]], n_samples)
    log_rates <- c(log_rates, list(fore[,sample_index]))
    dispersion <- c(dispersion, list(disp[sample_index,]))
  }
  print(cutpoints[i])
  rm(fore)
  rm(model)
  gc()
}

log_rates_a <- do.call(cbind, log_rates)
dispersion <- do.call(c, dispersion)

saveRDS(log_rates_a, paste0(path,"/stacked_forecast.rds"))
saveRDS(dispersion, paste0(path, "/dispersion.rds"))
