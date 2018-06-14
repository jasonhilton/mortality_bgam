#' Get all components of gam model
#' 
#' @param mort_fit The stanfit object containing the fitted model
#' @param stan_data The data used to fit the model
#' @param output_sex Character string indicating which sex to extract if the 
#' model includes both sexes jointly fitted. Default
#'
#' @return A named list with each element a matrix of samples of a parameter
get_components <- function(mort_fit, stan_data,
                           output_sex=c("single_sex",
                                        "joint_male",
                                        "joint_female")){
  output_sex <- match.arg(output_sex)
  smooth_mu <- get_smooth_function(mort_fit, stan_data$age_basis, "mu",
                                   output_sex=output_sex)
  smooth_aa <- get_smooth_function(mort_fit, stan_data$age_basis, "aa",
                                   output_sex=output_sex)
  smooth_gam <- get_smooth_function(mort_fit, stan_data$cohort_basis, "gam",
                                    output_sex=output_sex)
  
  kk <- t(as.matrix(mort_fit, "kk"))
  if (output_sex!="single_sex"){
    start_index = ifelse(output_sex=="joint_female", 1, 2)
    indicies = seq(start_index, dim(kk)[1], 2)
    kk <- kk[indicies,]
  }
  
  components <- list(smooth_mu= smooth_mu,
                     smooth_aa= smooth_aa,
                     smooth_gam= smooth_gam,
                     kk = kk)
  
  return(components)
}

#' Get dataframe of all components of gam model
#' 
#' @param components named list containing matrices of samples of each parameter
#' @param stan_data The data used to fit the model
#' @param end_year 
#' @param diff_period Boolean indicating whether to difference period effects
#' @param diff_cohort Boolean indicating whether to difference cohort effects
#' 
#' @return A single data frame containing all components, with columns x, fx, 
#' Component and Sim
get_all_components_df <- function(components, stan_data, end_year, 
                                  diff_period=T, diff_cohort=T){
  
  # Age effect -----------------------------------------------------------------
  
  mu_df <- get_component_df(components$smooth_mu, covar_name = "x", 
                            covar_values = 1:(stan_data$cutoff-1),
                            component_name="fx") %>% 
    mutate(Component="Age")
  
  # Improvements ---------------------------------------------------------------
  
  aa_df <- get_component_df(components$smooth_aa, covar_name = "x", 
                            covar_values = 1:(stan_data$cutoff-1),
                            component_name="fx") %>% 
    mutate(Component="Improvement")
  
  # adjust for standardised time index
  adjust_factor <- sd(1:stan_data$n_years)
  aa_df %<>% mutate(fx = fx / adjust_factor)
  
  # Cohorts --------------------------------------------------------------------
  
  smooth_gam <-components$smooth_gam
  if (diff_cohort){
    smooth_gam <- rbind(smooth_gam[1,], apply(smooth_gam, 2, diff))
  }
  
  cohorts <- ((end_year - stan_data$N - stan_data$n_years + 
                 stan_data$first_cohort + 1): 
                (end_year + stan_data$n_forecast_years))
  gam_df <- get_component_df(smooth_gam, covar_name = "x", 
                             covar_values = cohorts,
                             component_name="fx") %>%
    mutate(Component="Cohort")
  
  # Periods --------------------------------------------------------------------
  
  if (diff_period){
    kk <- apply(components$kk,2,diff)
  } else {
    kk <- components$kk
  }
  
  years <- ((end_year - stan_data$n_years + 1 + diff_period):
              (end_year + stan_data$n_forecast_years))
  kk_df <- get_component_df(kk, covar_name = "x", 
                            covar_values = years,
                            component_name="fx") %>%
    mutate(Component="Period")
  
  component_df <- rbind(mu_df, aa_df, gam_df, kk_df)
  return(component_df)
}


#' Extract posterior samples of the log-rates for the main gam model
#' 
#' @param components Named list of gam component samples
#' @param stan_data The data used to fit the model
#'
#'
#' @return The array where each matrix slice contains one posterior sample of
#' of the age x time lexis surface of log-mortality rates
get_young_log_rates <- function(components, stan_data){
  mu_rep <- (rep(1, stan_data$n_years + stan_data$n_forecast_years) %o% 
             components$smooth_mu)
  kk_rep <- rep(1, stan_data$cutoff - 1) %o% components$kk
  time_std <- get_standardised_time_index(stan_data$n_years, 
                                          stan_data$n_forecast_years)
  imp_mat <- time_std %o% components$smooth_aa 
  cohort_effects <- get_cohorts(components$smooth_gam, stan_data)
  log_rates <- mu_rep + imp_mat + cohort_effects[,1:(stan_data$cutoff-1),]
  log_rates <- aperm(log_rates, c(2,1,3)) + kk_rep
  return(log_rates)
}


#' Get array of cohort effects
#' 
#' @param smooth_gam Matrix with each row a sample from the posterior of the 
#'  smooth cohort function
#' @param stan_data Data used to fit the stan model
#' 
#' @return array with each matrix slice containing the time x age lexis surface
#' of cohort effects
get_cohorts <- function(smooth_gam, stan_data){
  iters <- dim(smooth_gam)[2]
  pad <- matrix(0, 
                stan_data$top_age + 1 - stan_data$N + stan_data$first_cohort,
                iters)
  smooth_gam <- rbind(pad, smooth_gam)
  cohort_effects <- array(0, c(stan_data$n_years + stan_data$n_forecast_years,
                               stan_data$top_age + 1,
                               iters))
  total_years <- stan_data$n_years + stan_data$n_forecast_years
  for(x in 1:(stan_data$top_age + 1)){
    cohort_effects[,x,] <- smooth_gam[(1 + stan_data$top_age + 1 - x):
                                     (total_years + stan_data$top_age + 1 - x),]
  }
  return(cohort_effects)
}

#' Calculate time index transformed to have zero mean and unit standard 
#' deviation
#' 
#' @param n_years Number of years in the data
#' @param n_forecast_years  number years to forecast for
#' 
#' @return Vector standardised using mean and sd of the in-sample time index.
get_standardised_time_index <- function(n_years, n_forecast_years){
  time <- 1:(n_years)
  time_fore <- 1:(n_years + n_forecast_years)
  time_std <- (time_fore - mean(time))/sd(time)
  return(time_std)
}


#' Extract posterior samples of smooth gam component
#'
#' @param mort_fit A stanfit object
#' @param basis The relevant basis function
#' @param parameter A character string containing the name of parameter as 
#' defined in the stan model
#' 
#' @return A matrix containing with each column containing a sample of the smooth
#' gam component. 
get_smooth_function <- function(mort_fit, basis, parameter_name,
                                output_sex=c("single_sex",
                                             "joint_male",
                                             "joint_female")){
  
  parameter <- as.matrix(mort_fit, parameter_name)
  
  if (output_sex!="single_sex"){
    start_index = ifelse(output_sex=="joint_female", 1, 2)
    indicies = seq(start_index, dim(parameter)[2], 2)
    parameter <- parameter[,indicies]
  }
  
  if( dim(parameter)[2] != dim(basis)[2]){
    stop("Mismatch between basis function and coefficient vector length \n")
  }
  
  smooth <- basis %*% t(parameter)
  return(smooth)
}

#' Convert samples of a smooth from a matrix to a dataframe
#'
#' @param component_matrix Matrix of samples of a smooth function, such as
#' returned by \code{get_smooth_function}.
#' @param covar_name Character vector giving the name of the covariate
#' @param covar_values Numeric vector with values of the covariate
#' @param component_name Character string containing the name of the component in 
#' question
#'
#' @return A tibble with columns for the covariate, the sample number ("Sim")
#' the value of the smooth function
get_component_df <- function(component_matrix, covar_name, covar_values, 
                             component_name){
  component_df <- component_matrix %>% tibble::as_tibble() %>% 
    dplyr::mutate(!!covar_name := covar_values) %>% 
    tidyr::gather(Sim, !!component_name, -!!covar_name)
  return(component_df)
}


#' Extract posterior samples of the log-rates for the old age model
#' 
#' @param mort_fit The stanfit object containing the fitted model
#' @param stan_data The data used to fit the model
#'
#'
#' @return The array where each matrix slice contains one posterior sample of
#' of the age x time lexis surface of log-mortality rates
get_old_log_rates <- function(mort_fit, stan_data, components, 
                              output_sex=c("single_sex",
                                           "joint_male",
                                           "joint_female")){
  
  output_sex <- match.arg(output_sex)
  
  n_ages <- stan_data$top_age - stan_data$cutoff + 1
  ones_age <- rep(1, n_ages)
  kk_rep <- ones_age %o%  components$kk
  cohort_effects <- get_cohorts(components$smooth_gam, stan_data)
  
  old_base <- get_old_base(mort_fit, stan_data, output_sex)
  
  cohort_index <- (stan_data$cutoff + 1):(stan_data$top_age + 1)
  
  log_asymptote <- as.matrix(mort_fit, "log_asymptote")
  ones_time <- rep(1, stan_data$n_years  + stan_data$n_forecast_years)
  log_asymptote <- ones_age %o% ones_time %o% log_asymptote[,1] 
  log_old_rates <- (old_base + kk_rep + 
                      aperm(cohort_effects[,cohort_index,], c(2,1,3))
                      - log1p(exp(old_base - log_asymptote)))
  return(log_old_rates)
  
}

#' Get parametric logistic part of the old age model
#' 
#' @param mort_fit The stanfit object containing the fitted model
#' @param stan_data The data used to fit the model
#' 
#' @return Array where each matrix slice contains one sample of the parametric 
#' part of the old age model
get_old_base <- function(mort_fit, stan_data,
                         output_sex=c("single_sex",
                                      "joint_male",
                                      "joint_female")){
  
  output_sex <- match.arg(output_sex)
  
  time_std <- get_standardised_time_index(stan_data$n_years, 
                                          stan_data$n_forecast_years)
  n_ages <- stan_data$top_age + 1 - stan_data$cutoff
  ones_age <- rep(1, n_ages)
  age_std <- get_standardised_age_index(n_ages)
  
  
  
  if (output_sex != "single_sex"){
    index  <- ifelse(output_sex=="joint_female", 1, 2)
    a1_name <- paste0("a_old_1", index)
    
  } else {
    index <- 1 
    a1_name <- "a_old_1"
  }
  m_old_0 <- as.matrix(mort_fit, "m_old_0")[,index]
  m_old_1 <- as.matrix(mort_fit, "m_old_1")[,index]
  a_old_0 <- as.matrix(mort_fit, "a_old_0")[,index]
  
  a_old_1 <- as.matrix(mort_fit, a1_name)
  ones_time <- rep(1, length(time_std))
  old_base <- (ones_age %o% ones_time %o% m_old_0 + 
                 age_std %o% ones_time  %o% m_old_1 +
                 ones_age %o% time_std %o% a_old_0 +
                 age_std %o% time_std %o% a_old_1[,1])
  return(old_base)
}

get_standardised_age_index <- function(n_ages){
  ages <- 1:n_ages
  return((ages- mean(ages))/sd(ages))
}

#' Extract posterior samples of the log-rates for the infant model
#' 
#' @param mort_fit The stanfit object containing the fitted model
#' @param stan_data The data used to fit the model
#'
#'
#' @return The array where each matrix slice contains one posterior sample of
#' of the age x time lexis surface of log-mortality rates
get_infant_log_rates <- function(mort_fit, stan_data, components,
                                 output_sex=c("single_sex",
                                              "joint_male",
                                              "joint_female")){
  
  output_sex <- match.arg(output_sex)
  time_std <- get_standardised_time_index(stan_data$n_years, 
                                          stan_data$n_forecast_years)
  
  
  m0 <- as.matrix(mort_fit, "m0")
  a0 <- as.matrix(mort_fit, "a0")
  
  if (output_sex!="single_sex"){
    index <- ifelse(output_sex=="joint_female", 1, 2)
    m0 <- m0[,index, drop=F]
    a0 <- a0[,index, drop=F]
  }
  
  smooth_gam <- components$smooth_gam
  gam_index <-((dim(smooth_gam)[1] - length(time_std) + 1):(dim(smooth_gam)[1]))
  infant_log_rates <- (rep(1, stan_data$n_years + stan_data$n_forecast_years) 
                       %o% m0  + 
                       time_std %o% a0)
  infant_log_rates <- aperm(infant_log_rates, c(3,1,2))
  infant_log_rates[1,,] <- infant_log_rates[1,,] + smooth_gam[gam_index,]
  return(infant_log_rates)
}


#' Extract posterior samples of the log-rates
#' 
#' @param mort_fit The stanfit object containing the fitted model
#' @param stan_data The data used to fit the model
#'
#'
#' @return The array where each matrix slice contains one posterior sample of
#' of the age x time lexis surface of log-mortality rates
get_log_rates <- function(mort_fit, stan_data,
                          output_sex=c("single_sex",
                          "joint_male",
                          "joint_female")){
  
  output_sex <- match.arg(output_sex)
  
  components <- get_components(mort_fit, stan_data, output_sex)
  infant_log_rates <- get_infant_log_rates(mort_fit, stan_data, components,
                                           output_sex)
  young_log_rates <- get_young_log_rates(components, stan_data)
  old_log_rates <- get_old_log_rates(mort_fit, stan_data, components,
                                     output_sex)
  
  log_rates <- abind::abind(infant_log_rates, young_log_rates, old_log_rates, 
                     along=1)
  return(log_rates)
}

#' Convert and array of log rates to a matrix
#' 
#' @param log_rate_array An array of log rates of dimension (n_ages, n_years,
#' n_iters [or n_quantiles])
#' 
#' @return A matrix of log rates of dimension (n_ages.n_years, 7
#' n_iters [or quantiles])
#'
log_rates_array_to_matrix <- function(log_rate_array){
  n_ages <- dim(log_rate_array)[1]
  n_years <- dim(log_rate_array)[2]
  n_cols <- dim(log_rate_array)[3]
  log_rates <- matrix(log_rate_array, n_ages*n_years, n_cols)
  return(log_rates)
}


#' Covert from an array to tidy data format
#' 
#' @param log_rate_array An array of log rates of dimension (n_ages, n_years,
#' n_iters)
#' @param years Vector of years to which rates apply
#' @param quantilise Compute quantiles over iterations?
#' @param probs Probabilities for which to compute quantiles
#'
get_log_rates_df <- function(log_rate_array, years, quantilise=T,
                             probs=1:99/100){
  var_name <- ifelse(quantilise==T, "Quantile", "Sim")
  if (quantilise){
    log_rate_array %<>% apply(c(1,2), quantile, probs=probs) %>% 
      aperm(c(2,3,1))
  }
  n_ages <- dim(log_rate_array)[1]
  n_years <- dim(log_rate_array)[2]
  log_rates <- log_rates_array_to_matrix(log_rate_array)
  log_rates_df <- log_rates %>% as_tibble() %>%
    mutate(Age=rep(0:(n_ages - 1), n_years),
           Year=rep(years, rep(n_ages,n_years))) %>%
    gather(!!var_name, Log_rate, -Year, -Age)
  if (quantilise){
    log_rates_df %<>% mutate(Quantile = as.numeric(gsub("V", "", Quantile))/100)
  }
  return(log_rates_df)
}

#' Extract log like-lihood from a saved model with a particular name 
#' 
#' 
get_log_lik <- function(cp, model_stem, path, sex){
  model <- readRDS(file.path(path, paste0(model_stem,"_", cp, ".rds" )))
  ll <- extract_log_lik(model$mort_fit, merge_chains=F)
  return(ll)
}

#' Extract log rates from model object and covert to long format.
#' @param mort_fit stan_fit object containing the relavent objects
#' @param stan_data 
#' @param sex The sex of the data the model refers to
#' 
#' @return A matrix of forecast log rates, with columns referring to sampels
#' and rows to particular ages and years.
get_forecast_log_rates_for_stacking <- function(mort_fit, stan_data, sex){
  fore_index <- ((stan_data$n_years + 1):
                   (stan_data$n_years + stan_data$n_forecast_years))
  if (tolower(sex)=="two_sex"){
    log_rates_m_a <- get_log_rates(mort_fit, stan_data, "joint_male")
    log_rates_f_a <- get_log_rates(mort_fit, stan_data, "joint_female")
    log_rates_fore_m <- log_rates_array_to_matrix(log_rates_m_a[,fore_index,])
    log_rates_fore_f <- log_rates_array_to_matrix(log_rates_f_a[,fore_index,])
    log_rates_fore <- rbind(log_rates_fore_f, log_rates_fore_m)
  } else {
    log_rates_a <- get_log_rates(mort_fit, stan_data, "single_sex")
    log_rates_fore <- log_rates_array_to_matrix(log_rates_a[,fore_index,])
  }
  return(log_rates_fore)
}

#' Get data frame rates including negative binomial uncertainty
#' 
#' @param log_rates
#' @param disp
#' @param emp_rate_df
#' @param top_age
#'
#'@return A data frame containing columns for expected deaths, sampled deaths,
#' implied rates and residuals.
get_assessment_df<- function(log_rates, disp, emp_rate_df, top_age, end_year){
  #n_total_years <- emp_rate_df$Year %>% unique() %>% length()
  #start_year <- min(emp_rate_df$Year)
  #end_year <- max(emp_rate_df$Year)
  n_ages <- top_age + 1
  n_years <- dim(log_rates)[1] / n_ages
  log_rates_df <- log_rates %>% as_tibble() %>%
    mutate(Age= rep(0:top_age, n_years),
           Year= rep((end_year - n_years + 1):end_year,
                     rep(top_age + 1, n_years))) %>%
    gather(Sim, Log_rate, -Year, -Age)
  disp_df <- tibble(Sim=1:length(disp), Dispersion=disp)
  log_rates_df <- left_join(log_rates_df %>%
                              mutate(Sim=as.numeric(gsub("V", "", Sim))),
                            disp_df)
  rates_df <- mutate(log_rates_df, Rate= exp(Log_rate)) %>%
    filter(Age < 111)
  rates_df <- left_join(rates_df,
                        emp_rate_df %>%
                          mutate(Obs_rate = Rate) %>%
                          select(-Rate, -Log_rate))
  rates_df %<>% mutate(E_deaths = Rate* Exposure,
                       S_deaths = rnbinom(n = n(), mu=E_deaths,size=Dispersion),
                       PP_Log_rate = log(S_deaths/Exposure),
                       S_resids = Deaths - S_deaths)
  
  return(rates_df)
}

#' Calculate quantiles of columns of a matrix or dataframe
#' 
#' @param df Data frame for which quantiles are to be calculated.
#' @param probs Probabilities to be calculated, passed to \code{quantile}
get_quantiles <- function(df, probs=1:99/100){
  df <- df %>% apply(2,function(x) quantile(x, probs, na.rm=T))
  df <- mutate(as_tibble(df), quantiles=as.numeric(gsub("%","",rownames(df)))/100)
  return(df)
}


get_qx <- function(mx, a){
  mx / (1+(1-a)*mx)
}

get_qx_alt <- function(mx){
  1-exp(-mx)
}

#' A function converting from the rates to qx death probabilities
#' 
#' @param log_rates_df a dataframe containing Age and Log_rate columns
#' 
#' @return A data frame containg and additional qx column.
get_qx_df <- function(log_rates_df){
  top_age <- max(log_rates_df$Age)
  ax_df <- data.frame(Age=0:top_age, ax=c(0.1,0.3,rep(0.5,top_age - 1)))
  log_rates_df %<>% left_join(ax_df) %>% 
    mutate(qx = ifelse(Age < 2, 
                       get_qx(exp(Log_rate),ax),
                       get_qx_alt(exp(Log_rate))))
  return(log_rates_df)
}

#' Calculate life expectancy at birth from qxs
#' 
#' @param qx_df A dataframe  containing Age, Year, qx, and Sim (the index of 
#' samples)
#' 
#' @return A dataframe containing life expectancies at birth for each Year and 
#' sample
get_e0 <- function(qx_df){
  model_e0_df <- qx_df %>%
    group_by(Year, Sim) %>%
    mutate(lx = cumprod(1-qx), lx_1=lag(lx)) %>% 
    mutate(lx_1=ifelse(Age==0, 1,lx_1)) %>% 
    mutate(Lx = (lx*ax + lx_1*(1-ax))) %>% 
    mutate(ex=rev(cumsum(rev(Lx)))) %>% ungroup() %>% 
    filter(Age==0) %>% 
    select(Year, Sim, ex) %>% rename(e0=ex)
  return(model_e0_df)
}
