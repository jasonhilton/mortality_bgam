library(tibble)
library(magrittr)
library(dplyr)
library(rstan)
library(loo)

source("R/data_construction.R")
source("R/utilities.R")

cmd_arg <- commandArgs(trailingOnly = T)

# if command line arguments are specified use these to determine config file 
# sex, and run_no (cutpoint). If not, use defaults in default_config.yaml
config_file <- ifelse(is.na(cmd_arg[1]), "default_config.yaml", cmd_arg[1])
config <- yaml::read_yaml(config_file)
sex <- ifelse(is.na(cmd_arg[2]), config$sex, cmd_arg[2])
run_no <- ifelse(is.na(cmd_arg[3]), config$run_no, as.numeric(cmd_arg[3]))

# for ease of reading
start_year <- config$start_year
end_year <- config$end_year
cutoff <- config$first_cut_point + run_no - 1
n_forecast_years <- config$n_forecast_years
model <- config$model

## Setup data --------------------------------------------------------
if (tolower(sex)=="two_sex"){
  expos_m <- get_data_mat(start_year,end_year,data_type="exposures", sex="Male")
  deaths_m <- get_data_mat(start_year,end_year,data_type="deaths", sex="Male")
  expos_f <- get_data_mat(start_year,end_year,data_type="exposures", sex="Female")
  deaths_f <- get_data_mat(start_year,end_year,data_type="deaths", sex="Female")
  params <- construct_params(cutoff = cutoff,
                             n_forecast_years=n_forecast_years,
                             period_forecast_covar = T
                             )
  stan_data <- get_2_sex_inputs(expos_m, deaths_m, expos_f, deaths_f, params)
  # check model name contains string "two_sex". If no match, grep returns 
  # a logical vector of length 0. 
  if (length(grep("two_sex", model)) < 1){
    stop("Specified model probably not suitable for two sex data\n",
         "Please use config file where the model argument refers to stan file",
         " containing a two-sex model.")
  }
} else {
  expos <- get_data_mat(start_year,end_year,data_type="exposures", 
                        sex=capwords(sex, strict=T))
  deaths <- get_data_mat(start_year,end_year,data_type="deaths", 
                         sex=capwords(sex, strict=T))
  
  params <- construct_params(cutoff = cutoff,
                             n_forecast_years=n_forecast_years)
  stan_data <- get_mortality_inputs(expos, deaths, params)
}


## Compile and run model ---------------------------------------------
model_path <- file.path("stan", paste0(model, ".stan"))

mort_fit <- stan(file=model_path, data=stan_data,
                   cores=config$cores,
                   chains=config$chains,
                   iter=config$iters,
                   thin=config$thin,
                   save_warmup=F)

## Save Results ------------------------------------------------------
out_path <- file.path(config$out_dir, config$type, tolower(sex))
dir.create(out_path, recursive=T)
saveRDS(list(mort_fit=mort_fit, stan_data=stan_data), 
        file=file.path(out_path, paste0(config$model, "_", cutoff, ".rds")))


