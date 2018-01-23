library(magrittr)
#' Construct data for stan mortality model
#' 
#' Function for creating stan input data for sampling from mortality model.
#'
#' @param expos The matrix of exposures, with ages = rows, and years = cols
#' @param deaths The matrix of deaths, with ages = rows, and years = cols
#' @param params A list of parameters controlling the construction of the data
#' inputs, including the gap between knots, the number of forecast years, and
#' the cutoff point between the young and old age model. Constructed by 
#' \code{construct_params}.
#' 
#' 
#' @seealso construct_params
#' 
#' @return A list to be passed as a data argument to stan(.)
#'
get_mortality_inputs <- function(expos, deaths, params = construct_params()){
  
  # TODO validation on expos and deaths matrix.
  
  # index vars
  N <- dim(deaths)[1]
  n_years <- dim(deaths)[2]
  
  # age gam basis data
  ages <- 1:(params$cutoff - 1)
  age_basis <- get_static_basis_function(ages, params$age_knot_gap)
  n_age_basis <- min(dim(age_basis))
  penalty_matrix_age <- make_penalty_matrix(n_age_basis)
  penalty_matrix_null <- make_null_penalty(penalty_matrix_age)
  
  # cohort gam basis data
  cohorts <- params$first_cohort:(N + n_years + params$n_forecast_years - 1)
  cohort_basis <- get_static_basis_function(cohorts, params$cohort_knot_gap)
  n_cohort_basis <- min(dim(cohort_basis))
  # note that constraints are relaxed for forecast years
  cohort_cons <- get_cohort_conditional_cov_matrix(cohort_basis, params)
  
  # use identity matrix as basis
  if (params$period_forecast_covar==T){
    period_cons <- get_period_conditional_cov_matrix(n_years, 
                                                     params$n_forecast_years)
  } else{
    period_cons <- get_period_conditional_cov_matrix(n_years)
  }
  # make output list
  stan_data <- list(N=N, n_years=n_years, top_age=params$top_age,
                    # spline basis function indexers
                    n_basis_age=n_age_basis, n_basis_cohort=n_cohort_basis,  
                    #basis functions
                    age_basis=age_basis, cohort_basis=cohort_basis,
                    k_inv_constraint = period_cons$inv_constraint,
                    Tau_k = period_cons$Tau,
                    gam_inv_constraint = cohort_cons$inv_constraint,
                    Tau_gam = cohort_cons$Tau,
                    expos = expos,                 # exposures
                    deaths = round(deaths),        # counts
                    penalty_matrix_age = penalty_matrix_age,
                    penalty_matrix_null=penalty_matrix_null,
                    n_forecast_years=params$n_forecast_years,
                    cutoff=params$cutoff,
                    first_cohort=params$first_cohort)
  return (stan_data)
}

#' Construct data for two sex stan mortality model
#' 
#' Function for creating stan input data for sampling from mortality model.
#'
#' @param expos_m The matrix of exposures for males, with ages = rows, and years
#'  = cols
#' @param deaths_m The matrix of deaths for males, with ages = rows, and years =
#'  cols
#' @param expos_f The matrix of exposures for females, with ages = rows, and 
#'  years = cols
#' @param deaths_f The matrix of deaths for females, with ages = rows, and years
#'  = cols
#' @param params A list of parameters controlling the construction of the data
#' inputs, including the gap between knots, the number of forecast years, and
#' the cutoff point between the young and old age model. Contructed by 
#' \code{construct_params}.
#' 
#' 
#' @seealso construct_params
#' 
#' @return A list to be passed as a data argument to stan(.)
#'
get_2_sex_inputs <- function(expos_m, deaths_m, expos_f, deaths_f,
                             params=construct_params(period_forecast_covar = T))
  {
  stan_data <- get_mortality_inputs(expos_f, deaths_f, params)
  stan_data$n_sexes <- 2
  expos <-abind::abind(list(expos_f, expos_m), along=3)
  expos <- aperm(expos, c(3,1,2))
  deaths <- abind::abind(list(round(deaths_f), round(deaths_m)),
                       along=3)
  deaths <- aperm(deaths, c(3,1,2))
  
  stan_data$deaths <- deaths
  stan_data$expos <- expos
  return(stan_data)
}


#' Construct parameters for stan mortality model
#' 
#' Create a list of parameters controlling the construction of the stan input
#' data, using the get_mortality_inputs function
#' 
#' @param cutoff The age at which the model transitions between the two separate
#' model
#' @param first_cohort The first cohort estimated using the gam model, 
#' considering the first cohort in the data as one. For earlier cohorts, all 
#' cohort effects are estimated to have the same value.
#' @param age_knot_gap The number of years of age between adjacent b-spline 
#' knots for the age smooth
#' @param cohort_knot_gap The number of cohorts between adjacent b-spline 
#' knots for the cohort smooth
#' @param n_forecast_years The number of years into the future to forecast
#' period and cohort effects.
#' @param top_age The highest age for which rates should be constructed.
#' @param zero_pad Number of zeros to add to the beginning of the constraint 
#' vectors
#' @param period_forecast_covariance Include future values in the construction 
#' of the period covariance matrix. Used for two-sex models.
#' 
#' @return A list of the parameters, including default values
#' 
#' @seealso get_mortality_inputs
#'
construct_params <- function(cutoff=90, first_cohort=6,
                             age_knot_gap=4, cohort_knot_gap=4,
                             n_forecast_years=0, top_age=125,
                             zero_pad = 0,
                             period_forecast_covar=F
){
  params <- as.list(environment())
  #  todo - some validation here.
  return(params)
}


#' Construct constraints for a basis matrix
#' 
#' Sometimes the output of a basis function might need to be constrained. This 
#' function constructs the constraint matrix needed to achieves this.
#' Combinations of constraints can be specificed. 
#' 
#' @param basis_matrix The basis matrix to be constrained
#' @param sum_zero Switch indicating whether sum-to-zero constraint should hold
#' @param zero_growth Switch indicating whether zero-growth constraint should hold
#' @param zero_quad Switch indicating whether zero-quadratic-growth constraint should hold
#' @param diff_sum_zero Switch indicating whether a sum-to-zero constraint 
#' should hold on the differenced scale
#' @param diff_zero_growth Switch indicating whether a zero-growth constraint 
#' should hold on the differenced scale
#' @param initial_zero Switch indicating whether the first value of the smooth 
#' function should be set equal to zero.
#' @param initial_value_equal Switch indicating whether the first and last values
#' should be equal
#' @param first_last_equal Switch indicating whether the first and last values
#' should be set to be equal.
#' @param last_zero Switch indicating whether the first and last values
#' should be set to be equal.
#' @param n_forecast how many years of forecasts are there, for which we do not
#' want constraints to hold
#' @param zero_pad number of zero to prepend to the constraint.
#' @return A constraint matrix for use in constructing a constrained basis
#' @export
#'
get_constraint <- function(basis_matrix, sum_zero=F, zero_growth=F, zero_quad=F,
                           diff_sum_zero=F, diff_zero_growth=F, 
                           initial_zero=F, initial_value_equal=F,
                           first_last_equal=F, last_zero=F,
                           n_forecast_years=0, zero_pad = 0){
  cons <- list()
  n_elems <- dim(basis_matrix)[1] - n_forecast_years - zero_pad
  n_basis <- dim(basis_matrix)[2]
  if (sum_zero){
    cons <- c(cons, list(c(rep(0,zero_pad),rep(1,n_elems), rep(0,n_forecast_years)) 
                         %*% basis_matrix))
  }
  if (zero_growth){
    cons <- c(cons ,list(c(rep(0,zero_pad), seq(- (n_elems - 1)/2,(n_elems - 1)/2),
                           rep(0,n_forecast_years))  %*% basis_matrix)  )
  }
  if (diff_zero_growth){
    cons <- c(cons, list(c(rep(0,zero_pad), 
                           (n_elems - 2) / 2, rep(-1, n_elems - 2),
                           (n_elems - 2) /2, 
                           rep(0, n_forecast_years)) %*% basis_matrix))
  }
  if (diff_sum_zero){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           -1, rep(0, n_elems - 2), 1, rep(0,n_forecast_years))
                         %*% basis_matrix ))
  }
  if (zero_quad){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           seq(-(n_elems - 1)/2, (n_elems - 1)/2)**2,
                           rep(0, n_forecast_years)) %*% basis_matrix))
  }
  if (initial_zero){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           1, rep(0, n_elems + n_forecast_years - 1)) 
                         %*% basis_matrix))
  }
  if (initial_value_equal){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           1,-1, rep(0, n_basis + n_forecast_years - 2))))
  }
  if (first_last_equal){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           1,rep(0, n_elems - 2),-1, rep(0,n_forecast_years)) 
                         %*% basis_matrix))
  }
  if (last_zero){
    cons <- c(cons, list(c(rep(0,zero_pad),
                           rep(0, n_elems-1), 1,
                           rep(0,n_forecast_years))
                         %*% basis_matrix))
  }
  
  C <- t(sapply(cons, function(x) x))
  return(C)
}


#' Convert a count dataframe to a matrix
#' 
#' The human mortality database uses long format for raw data. Convert to a 
#' matrix for ease of use.
#' 
#' @param count_df Dataframe with Year, Age and Count columns
#' @param start_year First year of data required
#' @param end_year Final year of data required
#' 
#' @return a matrix with ages in the rows and years across the columns.
#'
convert_to_matrix<-function(count_df, start_year, end_year){
  # expects count_df to have columns "Year", "Age" and "Count"
  
  # convert to Year by Age df, and filter to years required
  count_wide <- tidyr::spread(count_df, Age, Count)
  count_wide <- dplyr::filter(count_wide, Year>=start_year & Year<=end_year)
  
  row.names(count_wide) <- count_wide$Year
  count_wide <-  dplyr::select(count_wide, -Year)
  # for safety, check everything is in the right order
  count_wide <- count_wide[, order(as.numeric(names(count_wide)))]
  count_mat <- as.matrix(count_wide)
  return(count_mat)
}



#' Read in HMD data from disk and process
#' 
#' Read in a previously downloaded hmd dataframe
#' 
#' @param start_year First year of data required
#' @param end_year Final year of data required
#' @param sex Sex for which data is required
#' @param prefix string to prepend before the standard data path. This is useful
#' when processing a Rmd file that does not sit at the top level of the 
#' repository.
#' @param data_type Data required - exposures or deaths.
get_data_mat <- function(start_year, end_year,  sex=c("Male","Female"), 
                         prefix="./", data_type=c("exposures", "deaths")){
  data_type <- match.arg(data_type)
  sex <- match.arg(sex)
  data_hmd <- readRDS(file.path(prefix, "data", "hmd", 
                                 paste0(data_type, "_hmd.rds")))
  data <- data_hmd %>% dplyr::select(Year, Age, !!sex) %>% 
    dplyr::rename(Count:=!!sex)
  data_mat <- t(convert_to_matrix(data, start_year, end_year))  
  return(data_mat)
}

#' Get B-spline basis function with a fixed gap between knots
#' 
#' @param xx Covariate for which spline basis is to be defined
#' @param knot_gap The space between adjacent knots, using the same units as 
#' \code{xx}
#' 
#' @return A matrix of basis functions, with three knots outside the range of the 
#' data at each end.
get_static_basis_function <- function(xx, knot_gap){
  max_x <- max(xx)
  first_x <- min(xx)
  knot_locations <-  seq(first_x - (knot_gap*3), max_x + 1 + knot_gap*3, knot_gap)
  basis <- splines::splineDesign(knot_locations, xx, outer.ok=T)
  return(basis)
}

#' Return a penalty matrix for first differences.
#' 
#' @param n_basis the number of bases to be epanlised.
#' 
#' @return A matrix that can be used to penalise first differences.
make_penalty_matrix <- function(n_basis){
  penalty_matrix <- matrix(0, n_basis,n_basis)
  penalty_matrix <- penalty_matrix + diag(2, n_basis)
  diag(penalty_matrix[1:(n_basis-1), 2:n_basis]) <- rep(-1, n_basis - 1)
  diag(penalty_matrix[2:n_basis, 1:(n_basis - 1)]) <- rep(-1, n_basis - 1)
  penalty_matrix[1,1] <- 1
  penalty_matrix[n_basis,n_basis] <- 1
  return(penalty_matrix)
}

#' Find a penalty for the null space of an existing penalty matrix. See
#' \code{mgcv} package and the documentation for jagam.
#' 
#' @param penalty_matrix An existing penalty matrix, for example from 
#' \code{make_penalty_matrix}
#' 
#' @return A matrix that penalises the null space of \code{penalty_matrix}
make_null_penalty <- function(penalty_matrix){
  u0 <- eigen(penalty_matrix)$vector[,dim(penalty_matrix)[1]]
  penalty_matrix_null <- u0 %*% t(u0)
  return(penalty_matrix_null)
}

#' Construct the first-difference matrix operator
#' 
#' Returns square matrix $D$ with rows and columns $n$ that computes the first
#' difference of anything it is post-multiplied by.
#' 
#' @param n The number of rows in the difference matrix.
#' 
#' @return An $n$ x $n$ matrix with 1 on the diagonal and -1 on the first 
#' lower off-diagnoal
get_difference_matrix <- function(n){
  D <- matrix(0, n, n)
  D[2:n, 1:(n - 1)] <- - diag(n - 1)
  D <- D + diag(n)
  return(D)
}


#' Construct matrix to allow conditioning on constraints holding.
#' 
#' Constructs matrix that transforms from the initial differenced parameter 
#' space to one where some rows are zero if the constraints on the cumulative 
#' sums of the parameters (specifed in `constraint`) hold, allow conditioning on
#' the constraints holding true.
#' 
#' @param constraint  A matrix where each row is a constraint which must equal 
#' zero.
#' 
#' @returns A square matrix for which, for $k$ constraints, will contain $k$  
#' rows that will sum to zero if the desired constraints on the cumulative sums 
#' hold.
#'
get_constraint_transformation_matrix <- function(constraint, con_ind){
  n_basis <- dim(constraint)[2]
  ZZ <- diag(n_basis)
  ZZ[con_ind,] <- constraint
  return(ZZ)
}


#' Get the conditional covariance matrix for the cohort parameters,
#'
#' A distribution for cohort basis function coefficient innovation is needed for
#' model fitting, conditional on the identifiability constraints on the cohort 
#' effects. This function constructs the covariance matrix of this 
#' distribution for use in model fitting. This covariance matrix may be scaled 
#' by a variance parameter during the sampling process.
#'
#' @param cohort_basis The matrix of cohort basis functions.
#' @param params set of parameters constructed using construct_params
get_cohort_conditional_cov_matrix <- function(cohort_basis, params){
  n_cohort_basis <- min(dim(cohort_basis))
  # note that constraints are relaxed for forecast years
  constraint <- get_constraint(cohort_basis, sum_zero = T,
                               initial_zero = T, last_zero = T,
                               n_forecast_years = params$n_forecast_years,
                               zero_pad = params$zero_pad)  
  # swap sum to zero constraint and initial zero constraint.
  constraint <- constraint[c(2,1,3),]
  # get cumulative sum matrix
  S <- solve(get_difference_matrix(n_cohort_basis))
  CS <- constraint %*% S
  non_zero_basis <- which(constraint[2,] > 0)
  # How many pre-constraint parameters are there?
  i1 <- non_zero_basis[1] - 1
  # How many unconstrained (forecast) parameters are there? 
  i2 <- n_cohort_basis -  non_zero_basis[length(non_zero_basis)]
  # Which bases should be constrained 
  con_ind <- c(i1 + c(1,2), n_cohort_basis - i2)
  result <- get_conditional_cov_matrix(CS, con_ind)
  return(result)
}


#' Get the conditional covariance matrix for the period parameters,
#'
#' A distribution for the period innovations is needed for model fitting, 
#' conditional on the identifiability constraints. 
#' This function constructs the covariance matrix of this 
#' distribution for use in model fitting. This covariance matrix may be scaled 
#' by a variance parameter during the sampling process.
#'
#' @param n_years The number of years in the data.
#' @param n_forecast_years Number of forecast years to include in the covariance
#' matrix. Only useful in the two-sex case.
#' 
#' @return The raw covariance matrix for period effect, conditional on the 
#' constraints
get_period_conditional_cov_matrix <- function(n_years, n_forecast_years=0){

  constraint <- get_constraint(diag(n_years + n_forecast_years),
                               sum_zero = T, zero_growth=T,
                               n_forecast_years = n_forecast_years)
  # get cumulative sum matrix
  S <- solve(get_difference_matrix(n_years + n_forecast_years))
  CS <- constraint %*% S
  # Which bases should be constrained 
  con_ind <- c(1,2)
  result <- get_conditional_cov_matrix(CS, con_ind)
  return(result)
}

#' Construct a conditional covariance matrix
#' 
#' Given an unconditional covariance matrix, and the indexes to which the 
#' constraints apply, construct an unscaled conditional covariance matrix, based
#' on standard multivariate normal results. 
#' 
#' @param CS Unconditional covariance matrix
#' @param con_ind The indexes to which the constraints apply.
#'
get_conditional_cov_matrix <- function(CS, con_ind){
  ZZ <- get_constraint_transformation_matrix(CS, con_ind)
  Sigma <- ZZ %*% t(ZZ)
  Tau <- (Sigma[-con_ind,-con_ind] - Sigma[-con_ind,con_ind] %*%
            solve(Sigma[con_ind,con_ind]) %*%
            Sigma[con_ind,-con_ind])
  inv_constraint <- solve(ZZ)[,-con_ind]
  return(list(Tau=Tau, inv_constraint=inv_constraint))
}



