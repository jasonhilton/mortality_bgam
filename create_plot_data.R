library(rstan)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(tibble)
library(yaml)
library(loo)

#source("R/results_processing.R")


cmd_arg <- commandArgs(trailingOnly = T)

# if command line arguments are specified use these to determine results_dir
results_dir <- ifelse(is.na(cmd_arg[1]), "results", cmd_arg[1])


source("R/results_processing.R")
source("R/utilities.R")

## Components ------------------------------------------------------------------
results_m_path <- file.path(results_dir, "standard", "male")
results_f_path <- file.path(results_dir, "standard", "female")

## TODO double check cutpoints are the right ones
model_m <- readRDS(file.path(results_m_path, "mortality_ll_93.rds"))
model_f <- readRDS(file.path(results_f_path, "mortality_ll_91.rds"))

comp_m <- get_components(model_m$mort_fit, model_m$stan_data,
                         output_sex="single_sex")
comp_m_df <- get_all_components_df(components = comp_m,
                                   stan_data=model_m$stan_data,
                                   end_year = 2013)

comp_f <- get_components(model_f$mort_fit, model_f$stan_data,
                         output_sex="single_sex")
comp_f_df <- get_all_components_df(components = comp_f,
                                   stan_data=model_f$stan_data,
                                   end_year = 2013)

process_dir <- file.path(results_dir, "processed")
dir.create(process_dir)

saveRDS(comp_m_df,file=file.path(process_dir,"comp_m_df.rds"))
saveRDS(comp_f_df,file=file.path(process_dir,"comp_f_df.rds"))
rm(comp_m_df, comp_f_df)
gc()

## log rates for single-run plots (quantilised) --------------------------------

log_rate_m_a <- get_log_rates(mort_fit = model_m$mort_fit, stan_data=model_m$stan_data,
                            output_sex = "single_sex")
log_rate_m_df <- get_log_rates_df(log_rate_m_a, years = 1961:2063, quantilise = T)
saveRDS(log_rate_m_df, file.path(process_dir,"log_rate_m_93.rds"))
rm(log_rate_m_a, log_rate_m_df, model_m)

log_rate_f_a <- get_log_rates(mort_fit = model_f$mort_fit, stan_data=model_f$stan_data,
                              output_sex = "single_sex")
log_rate_f_df <- get_log_rates_df(log_rate_f_a, years = 1961:2063, quantilise = T)
saveRDS(log_rate_f_df, file.path(process_dir,"log_rate_f_91.rds"))
rm(log_rate_f_a, log_rate_f_df, model_f)
gc()

## log rates for stacked forecasts (quantilised) -------------------------------

# only two dimensions in stacked forecast.

log_rate_m <- readRDS(file.path(results_m_path,"stacked_forecast.rds"))
config <- yaml::read_yaml("config/single_sex.yaml")
n_years <- config$n_forecast_years
n_iters <- dim(log_rate_m)[2]
n_ages <- dim(log_rate_m)[1] / n_years
log_rate_m_a <- array(log_rate_m, c(n_ages, n_years, n_iters))
log_rate_m_df <- get_log_rates_df(log_rate_m_a,
                                  years = ((config$end_year + 1):
                                           (config$end_year + n_years)), 
                                  quantilise = T)
saveRDS(log_rate_m_df, file.path(process_dir,"log_rate_m_stacked.rds"))
rm(log_rate_m, log_rate_m_a, log_rate_m_df)
gc()

log_rate_f <- readRDS(file.path(results_f_path,"stacked_forecast.rds"))
log_rate_f_a <- array(log_rate_f, c(n_ages, n_years, n_iters))
log_rate_f_df <- get_log_rates_df(log_rate_f_a,  
                                  years = ((config$end_year + 1):
                                           (config$end_year + n_years)), 
                                  quantilise = T)
saveRDS(log_rate_f_df, file.path(process_dir,"log_rate_f_stacked.rds"))
rm(log_rate_f_a, log_rate_f, log_rate_f_df)

## log rates for two-sex stacked forecasts (quantilised) -----------------------
results_two_sex_path <- file.path(results_dir, "standard", "two_sex")
log_rate <- readRDS(file.path(results_two_sex_path, "stacked_forecast.rds"))

# first half are female
halfway <- dim(log_rate)[1]/2
end <- dim(log_rate)[1]

log_rate_2f <- log_rate[1:halfway,]
log_rate_2m <- log_rate[(halfway + 1):end,]
rm(log_rate)
gc()
config <- yaml::read_yaml("config/two_sex.yaml")
n_years <- config$n_forecast_years
n_iters <- dim(log_rate_2m)[2]
n_ages <- dim(log_rate_2m)[1] / n_years

log_rate_2f_a <- array(log_rate_2f, c(n_ages, n_years, n_iters))
log_rate_f_df <- get_log_rates_df(log_rate_2f_a,  
                                  years = ((config$end_year + 1):
                                             (config$end_year + n_years)), 
                                  quantilise = T)
saveRDS(log_rate_f_df, file.path(process_dir,"log_rate_2_f_stacked.rds"))
rm(log_rate_f_df, log_rate_2f_a, log_rate_2f)

log_rate_2m_a <- array(log_rate_2m, c(n_ages, n_years, n_iters))
log_rate_m_df <- get_log_rates_df(log_rate_2m_a,  
                                  years = ((config$end_year + 1):
                                             (config$end_year + n_years)), 
                                  quantilise = T)
saveRDS(log_rate_m_df, file.path(process_dir,"log_rate_2_m_stacked.rds"))
rm(log_rate_m_df, log_rate_2m_a, log_rate_2m)
gc()


## parameter comparison - fitting period ---------------------------------------

cutpoints <- 80:95

get_components_mean <- function(cp, results_dir, model_stem, output_sex, end_year){
  model <- readRDS(file.path(results_dir,paste0(model_stem, "_",cp,".rds")))
  comp <- get_components(model$mort_fit, model$stan_data,output_sex)
  comp_df <- get_all_components_df(comp, model$stan_data, end_year)
  comp_df %<>% group_by(Component,x) %>% summarise(fx=mean(fx)) %>% 
    mutate(Cutpoint=cp)
  return(comp_df)
}

comp_m <- map_df(cutpoints, get_components_mean, results_m_path, "mortality_ll", 
              "single_sex", 2013) %>%
  mutate(`Fitting Period`="1961-2013") %>%
  filter((Component!="Period" & Component!="Cohort")| x < 2014)

results_m_hb_path <- file.path(results_dir, "holdback", "male")
comp_m_hb <- map_df(cutpoints, get_components_mean, results_m_hb_path, "mortality_ll", 
              "single_sex", 2003) %>% 
  mutate(`Fitting Period`="1961-2003")

comparison_df <- rbind(comp_m, comp_m_hb)
saveRDS(comparison_df, file = file.path(process_dir,"comparison.rds"))


## holdback comparison joint vs single -----------------------------------------
# requires additional of neg_bin uncertainty.


## single-----------------------------------------------------------------------
emp_rate_df <- load_empirical_rates()
# male

log_rates_m_a <- readRDS(file=file.path(results_m_hb_path, "stacked_forecast.rds"))
disp_m <- readRDS(file=file.path(results_m_hb_path, "dispersion.rds"))
ass_m_df <- get_assessment_df(log_rates=log_rates_m_a,
                              disp = disp_m,top_age = 125,
                              emp_rate_df = emp_rate_df %>% filter(Sex=="Male"),
                              end_year = 2013)
ass_q_m <- ass_m_df %>% select(-Sex,-Sim) %>% nest(-Age, -Year) %>%
  mutate(quantiles=map(data, get_quantiles)) %>% select(-data) %>% unnest()
rm(ass_m_df)


results_f_hb_path <- file.path(results_dir, "holdback", "female")
# female
log_rates_f_a <- readRDS(file=file.path(results_f_hb_path, "stacked_forecast.rds"))
disp_f <- readRDS(file=file.path(results_f_hb_path, "dispersion.rds"))
ass_f_df <- get_assessment_df(log_rates=log_rates_f_a,
                              disp = disp_f,top_age = 125,
                              emp_rate_df = emp_rate_df %>% filter(Sex=="Female"),
                              end_year = 2013)
ass_q_f <- ass_f_df %>% select(-Sex,-Sim) %>% nest(-Age, -Year) %>%
  mutate(quantiles=map(data, get_quantiles)) %>% select(-data) %>% unnest()
rm(ass_f_df)


## joint------------------------------------------------------------------------
#
results_two_sex_hb_path <- file.path(results_dir, "holdback", "two_sex")

log_rate_a <- readRDS(file=file.path(results_two_sex_hb_path, "stacked_forecast.rds"))
# first half are female
halfway <- dim(log_rate_a)[1]/2
end <- dim(log_rate_a)[1]
log_rates_2f_a <- log_rate_a[1:halfway,]
log_rates_2m_a <- log_rate_a[(halfway + 1):end,]
disp <- readRDS(file=file.path(results_two_sex_hb_path, "dispersion.rds"))
ass_f_df <- get_assessment_df(log_rates=log_rates_2f_a,
                              disp = disp,top_age = 125,
                              emp_rate_df = emp_rate_df %>% filter(Sex=="Female"),
                              end_year=2013)
ass_2q_f <- ass_f_df %>% select(-Sex,-Sim) %>% nest(-Age, -Year) %>%
  mutate(quantiles=map(data, get_quantiles)) %>% select(-data) %>% unnest()
rm(ass_f_df)

ass_m_df <- get_assessment_df(log_rates=log_rates_2m_a,
                              disp = disp,top_age = 125,
                              emp_rate_df = emp_rate_df %>% filter(Sex=="Male"),
                              end_year=2013)
ass_2q_m <- ass_m_df %>% select(-Sex,-Sim) %>% nest(-Age, -Year) %>%
  mutate(quantiles=map(data, get_quantiles)) %>% select(-data) %>% unnest()
rm(ass_m_df)

ass_q <- rbind(ass_q_m %>% mutate(Sex="Male"),
               ass_q_f %>% mutate(Sex="Female"))

ass_q2 <- rbind(ass_2q_m %>% mutate(Sex="Male"),
                ass_2q_f %>% mutate(Sex="Female"))

ass_all <- rbind(ass_q %>% mutate(Model="Single"),
                 ass_q2 %>% mutate(Model="Joint"))

# replace a couple of sampled 0s with very low rate
ass_all$PP_Log_rate[is.infinite(ass_all$PP_Log_rate)] <- -10


saveRDS(ass_all, file=file.path(process_dir,"ass_df.rds"))
rm(ass_all)
## qx and e0 plots vs ons. -------------------------------------------------------------
log_rate <- readRDS(file=file.path(results_two_sex_path, 
                                   "stacked_forecast.rds"))

halfway <- dim(log_rate)[1]/2
end <- dim(log_rate)[1]

log_rate_2f <- log_rate[1:halfway,]
log_rate_2m <- log_rate[(halfway + 1):end,]
rm(log_rate)
gc()
config <- yaml::read_yaml("config/two_sex.yaml")
n_years <- config$n_forecast_years
n_iters <- dim(log_rate_2m)[2]
n_ages <- dim(log_rate_2m)[1] / n_years

log_rate_2f_a <- array(log_rate_2f, c(n_ages, n_years, n_iters))
log_rate_f_df <- get_log_rates_df(log_rate_2f_a,  
                                  years = ((config$end_year + 1):
                                           (config$end_year + n_years)),
                                  quantilise = F)

rm(log_rate_2f_a, log_rate_2f)

log_rate_2m_a <- array(log_rate_2m, c(n_ages, n_years, n_iters))
log_rate_m_df <- get_log_rates_df(log_rate_2m_a,
                                  years = ((config$end_year + 1):
                                           (config$end_year + n_years)),
                                  quantilise = F)
rm(log_rate_2m_a, log_rate_2m)

qx_f_df <- get_qx_df(log_rate_f_df)
rm(log_rate_f_df)
e0_f <- get_e0(qx_df = qx_f_df)
qx_f_df %<>% select(-Sim, -Log_rate, -ax) %>% nest(-Age,-Year) %>%
  mutate(quantiles=map(data, get_quantiles)) %>% select(-data) %>% unnest()

qx_m_df <- get_qx_df(log_rate_m_df)
rm(log_rate_m_df)
e0_m <- get_e0(qx_df = qx_m_df)
qx_m_df %<>% select(-Sim, -Log_rate, -ax) %>% nest(-Age,-Year) %>%
  mutate(quantiles=map(data, get_quantiles)) %>% select(-data) %>% unnest()


qx_df <- rbind(qx_f_df %>% mutate(Sex="Female"),
               qx_m_df %>% mutate(Sex="Male"))

rm(qx_m_df, qx_f_df)
saveRDS(qx_df, file=file.path(process_dir,"qx_df.rds"))
rm(qx_df)
gc()

e0_df <- rbind(e0_f %>% mutate(Sex="Female"),
               e0_m %>% mutate(Sex="Male"))
saveRDS(e0_df, file=file.path(process_dir, "e0_df.rds"))

## save out session info.
library(abind)
library(devtools)
sesh_out <- session_info()
saveRDS(sesh_out, file.path(process_dir, "iridis_sesh.rds"))

