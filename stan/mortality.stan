data { 
  // indexes
  int N; // number of ages in data
  int n_years; // years of data
  int cutoff; // cutoff is the first age in the old age model
  int top_age; // the highest age for which to estimate rates (>= N)
  int n_basis_age; // number of basis functions of age
  int n_basis_cohort; // number of basis functions of age
  int n_forecast_years; // number of years to forecast
  int first_cohort; // first cohort in the model

  // empirical data
  int<lower=0> deaths[N, n_years];
  matrix<lower=0>[N, n_years] expos;

  // covariate bases, 
  matrix<lower=0>[cutoff - 1, n_basis_age] age_basis;
  matrix[(N + n_years - first_cohort + n_forecast_years),
         n_basis_cohort] cohort_basis;

  // period transformation and covariance matricies
  matrix[n_years, n_years - 2] k_inv_constraint;
  matrix[n_years - 2, n_years - 2] Tau_k; // unscaled covariance matrix for k

  // age penalty matrices
  matrix[n_basis_age, n_basis_age] penalty_matrix_age; 
  matrix[n_basis_age, n_basis_age] penalty_matrix_null;

  // cohort transformation and covariance matricies
  matrix[n_basis_cohort, n_basis_cohort - 3] gam_inv_constraint;
  matrix[n_basis_cohort - 3, n_basis_cohort - 3] Tau_gam; 
}

transformed data {
  matrix[N, n_years] log_expos;
  matrix[n_years - 2, n_years - 2] LK;
  matrix[n_basis_cohort - 3, n_basis_cohort - 3] Lgam;
  row_vector[n_years] time;
  row_vector[n_years] ones; // used in repeating matrix
  vector[n_years + n_forecast_years] zeroes;
  vector[top_age - cutoff + 1] old_age_index;

  // calculate log exposures in advance ----------------------------------------
  log_expos = log(expos);

  // cholesky decompose unscaled covariance mat for conditional period effects
  LK = cholesky_decompose(Tau_k);
  // cholesky decompose unscaled covariance mat for conditional cohort effects
  Lgam = cholesky_decompose(Tau_gam);

  // create standardised time index --------------------------------------------
  for (t in 1:(n_years)){
    time[t] = t; 
    ones[t] = 1; // used to create repeated matrices
  }
  time = (time - mean(time[1:n_years])) / sd(time[1:n_years]);

  // create standardised old age index -----------------------------------------
  for (x in 1:(top_age - cutoff + 1)){
    old_age_index[x] = x;
  }
  old_age_index = (old_age_index - mean(old_age_index)) / sd(old_age_index);

  // zeroes vector for use in prior as means.
  zeroes = rep_vector(0.0, n_years + n_forecast_years);
}

parameters {
  real log_dispersion;

  // static age-specific log-rate pattern mu -----------------------------------
  vector[n_basis_age] mu;
  real<lower=0> sigma_mu;
  real<lower=0> sigma_mu_null;

  // age-specific improvement rates  -------------------------------------------
  vector[n_basis_age] aa;
  real<lower=0> sigma_aa;
  real<lower=0> sigma_aa_null;

  // years-specific effects ----------------------------------------------------
  vector[n_years - 2] dkz;
  vector[n_forecast_years] delta_kk_fore;
  real<lower=0> sigma_dk;

  // cohort effects ------------------------------------------------------------
  vector[n_basis_cohort - 3] dgamz;
  real<lower=0> sigma_dgam;

  // infant parameters ---------------------------------------------------------
  real a0;
  real m0;

  // old age parameters --------------------------------------------------------
  real m_old_0;
  real<lower=0> m_old_1;
  real<upper=0> a_old_0;
  real<lower=-m_old_1 / max(time)> a_old_1;
  real log_asymptote;
}

transformed parameters {
  row_vector[n_years + n_forecast_years] kk;
  real asymptote;
  vector[n_basis_cohort] gam;
  real dispersion;

  // calculate cohort basis function coefficients ------------------------------
  {
  vector[n_basis_cohort] delta_gam;
  // equivalent to MVN distribution on delta_gam
  // as dgamz ~ normal(0, 1);
  delta_gam = gam_inv_constraint * sigma_dgam * (Lgam * dgamz);
  gam = cumulative_sum(delta_gam);
  }

  // convert asymptote and dispersion ------------------------------------------
  asymptote = exp(log_asymptote);
  dispersion = exp(log_dispersion);

  // calculate period effects
  {
  vector[n_years - 2] delta_kk;
  vector[n_years + n_forecast_years] delta_kk_full;
  // this VvV is equivalent to delta_kk ~ MVN(0, sigma.(LK'LK)), 
  // as dkz ~normal(0, 1);
  delta_kk = sigma_dk * (LK * dkz);
  delta_kk_full[1:n_years] = k_inv_constraint * delta_kk;
  delta_kk_full[(n_years + 1):(n_years + n_forecast_years)] = delta_kk_fore;
  kk = cumulative_sum(delta_kk_full * 0.1)';
  }
}
model {
  // holding matrices for intermediate results ---------------------------------
  vector[cutoff - 1] smooth_mu;
  matrix[cutoff - 1, n_years] smooth_imp;
  matrix[cutoff - 1, n_years] mu_rep;
  matrix[cutoff - 1, n_years] kk_rep;
  row_vector[top_age + n_years + n_forecast_years] smooth_gam;
  matrix[top_age + 1, n_years] cohort_effects;

  // final rate matrices
  matrix[(top_age - cutoff + 1), n_years] old_base_rate;
  matrix[(top_age - cutoff + 1), n_years] old_logit_rate;
  matrix[(top_age - cutoff + 1), n_years] old_log_rate;
  matrix[cutoff, n_years] young_log_rate;  

  // temporary rate variables
  real lambda;
  row_vector[n_years] lambda_v;

  // precision matrices for prior on basis function coefficients.
  matrix[n_basis_age, n_basis_age] full_penalty_mu;
  matrix[n_basis_age, n_basis_age] full_penalty_aa;
  matrix[n_basis_cohort, n_basis_cohort] full_penalty_gam;


  // priors --------------------------------------------------------------------
  
  // infants
  a0 ~ normal(0, 10);
  m0 ~ normal(0, 10);

  // old age
  m_old_0 ~ normal(0, 10);
  m_old_1 ~ normal(0, 10);
  a_old_0 ~ normal(0, 10);
  a_old_1 ~ normal(0, 10) T[- m_old_1 / max(time),]; // enforce positive gradient
  log_asymptote ~ normal(0, 1);

  // period priors 
  // Note: all sigma are half-normal, taking positive part only.
  sigma_dk ~ normal(0, 10);
  delta_kk_fore ~ normal(zeroes[1:n_forecast_years], sigma_dk);
  dkz ~ normal(0, 1);
  
  // cohort priors
  dgamz ~ normal(0,1);
  sigma_dgam ~ normal(0, 10);

  // basis function coefficient priors
  sigma_mu ~ normal(0, 10);
  sigma_mu_null ~ normal(0, 10);
  sigma_aa ~ normal(0, 10);
  sigma_aa_null ~ normal(0, 10);
  

  // implicit prior
  // log_dispersion ~ 1

  // penalties on rough vectors ------------------------------------------------
  
  full_penalty_mu = ((1 / pow(sigma_mu,2)) * penalty_matrix_age + 
                     (1 / pow(sigma_mu_null, 2)) * penalty_matrix_null); 

  full_penalty_aa = ((1 / pow(sigma_aa, 2)) * penalty_matrix_age + 
                     (1 / pow(sigma_aa_null, 2)) * penalty_matrix_null); 
  
  mu ~ multi_normal_prec(zeroes[1:n_basis_age], full_penalty_mu);
  aa ~ multi_normal_prec(zeroes[1:n_basis_age], full_penalty_aa);

  // year-specific effects -----------------------------------------------------
  kk_rep = rep_matrix(kk[1:n_years], cutoff - 1);

  // static age-specific log-rate pattern mu -----------------------------------
  smooth_mu = age_basis * mu; 
  mu_rep = rep_matrix(smooth_mu, n_years);

  // improvement rates ---------------------------------------------------------
  smooth_imp = age_basis * aa * time;


  // cohort --------------------------------------------------------------------
              // first cohort plus additional zeros above N
  smooth_gam[(first_cohort + ((top_age + 1) - N)):
               (top_age + n_years + n_forecast_years)] = (cohort_basis * gam)';
  
  for (i in 1:(first_cohort + (top_age - N))){
    smooth_gam[i] = 0;
  }  
  
  // fill in cohort matrix
  for (i in 1:(top_age + 1)){
    cohort_effects[i] = smooth_gam[(1 + (top_age + 1) - i):
                                   (n_years + (top_age + 1) - i)];
  }

  // data model ----------------------------------------------------------------

  // infants -------------------------------------------------------------------
  young_log_rate[1] = m0 + a0 * time + cohort_effects[1];
  
  // gam model -----------------------------------------------------------------
  young_log_rate[2:cutoff] = (mu_rep + smooth_imp + kk_rep + 
                              cohort_effects[2:cutoff]);

  // old age -------------------------------------------------------------------  
  old_base_rate = ((m_old_0 + m_old_1 * old_age_index) * ones + 
                   (a_old_0 + a_old_1  * old_age_index) * time);
  old_logit_rate = (old_base_rate  + kk_rep[1:(top_age - cutoff + 1)] + 
                    cohort_effects[(cutoff + 1):(top_age + 1)]);

  // log likelihoods -----------------------------------------------------------

  // main gam-based model ------------------------------------------------------
  for (i in 1:cutoff){
    lambda_v = young_log_rate[i] + log_expos[i];
    deaths[i] ~ neg_binomial_2_log(lambda_v, dispersion);
  }
  
  // old_age -------------------------------------------------------------------
  // calculate rates
  for (i in (cutoff + 1):(top_age + 1)){
    for (t in 1:n_years){
      old_log_rate[i - cutoff, t] = (old_logit_rate[i - cutoff, t] - 
                               log1p_exp(old_base_rate[i - cutoff, t] -
                                        log_asymptote));
    }
  }

  // old-age likelihood
  for (i in (cutoff + 1):N){
    for (t in 1:(n_years)){
      if (expos[i, t] != 0 && N - i + t >= first_cohort){
        deaths[i, t]  ~ neg_binomial_2_log(log_expos[i, t] +
                                           old_log_rate[i - cutoff, t], 
                                           dispersion);                    
      } 
    }
  }
}
