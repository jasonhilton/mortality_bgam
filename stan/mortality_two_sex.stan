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
  int n_sexes;

  // empirical data
  int<lower=0> deaths[n_sexes, N, n_years];
  matrix<lower=0>[N, n_years] expos[n_sexes];

  // covariate bases, 
  matrix<lower=0>[cutoff - 1, n_basis_age] age_basis;
  matrix[(N + n_years - first_cohort + n_forecast_years),
         n_basis_cohort] cohort_basis;

  // period transformation and covariance matricies
  matrix[n_years + n_forecast_years, 
         n_years + n_forecast_years - 2] k_inv_constraint;
  // unscaled covariance matrix for k
  matrix[n_years + n_forecast_years - 2,
         n_years + n_forecast_years - 2] Tau_k; 

  // age penalty matrices
  matrix[n_basis_age, n_basis_age] penalty_matrix_age; 
  matrix[n_basis_age, n_basis_age] penalty_matrix_null;

  // cohort transformation and covariance matricies
  matrix[n_basis_cohort, n_basis_cohort - 3] gam_inv_constraint;
  matrix[n_basis_cohort - 3, n_basis_cohort - 3] Tau_gam; 
}

transformed data {
  matrix[N, n_years] log_expos[n_sexes];
  matrix[n_years + n_forecast_years - 2, n_years + n_forecast_years - 2] LK;
  matrix[n_basis_cohort - 3, n_basis_cohort - 3] Lgam;
  row_vector[n_years] time;
  row_vector[n_years] ones; // used in repeating matrix
  vector[n_years + n_forecast_years] zeroes; // used for mean in MVN priors
  vector[top_age - cutoff + 1] old_age_index;
  int total_years;

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

  total_years = n_years + n_forecast_years;
}

parameters {
  real log_dispersion[n_sexes];

  // static age-specific log-rate pattern mu -----------------------------------
  vector[n_basis_age] mu[n_sexes];
  real<lower=0> sigma_mu[n_sexes];
  real<lower=0> sigma_mu_null[n_sexes];

  // age-specific improvement rates  -------------------------------------------
  vector[n_basis_age] aa[n_sexes];
  real<lower=0> sigma_aa[n_sexes];
  real<lower=0> sigma_aa_null[n_sexes];

  // years-specific effects ----------------------------------------------------
  vector[(n_years + n_forecast_years - 2) *n_sexes] dkz;
  real<lower=0,upper=1> rho_kk;
  vector<lower=0>[n_sexes] sigma_dk;

  // cohort effects ------------------------------------------------------------
  vector[n_basis_cohort - 3] dgamz[n_sexes];
  real<lower=0> sigma_dgam[n_sexes];

  // infant parameters ---------------------------------------------------------
  real a0[n_sexes];
  real m0[n_sexes];

  // old age parameters --------------------------------------------------------
  real m_old_0[n_sexes];
  real<lower=0> m_old_1[n_sexes];
  real<upper=0> a_old_0[n_sexes];
  real<lower=-m_old_1[1]/max(time)> a_old_11;
  real<lower=-m_old_1[2]/max(time)> a_old_12;
  real log_asymptote;
}

transformed parameters {
  row_vector[n_years + n_forecast_years] kk[n_sexes];
  real asymptote;
  vector[n_basis_cohort] gam[n_sexes];
  real dispersion[n_sexes];

  // calculate cohort basis function coefficients ------------------------------
  for (s in 1:n_sexes){

    vector[n_basis_cohort] delta_gam;
    // equivalent to MVN distribution on delta_gam
    // as dgamz ~ normal(0, 1);
    delta_gam = gam_inv_constraint * sigma_dgam[s] * (Lgam * dgamz[s]);
    gam[s] = cumulative_sum(delta_gam);
  }

  // convert asymptote and dispersion ------------------------------------------
  asymptote = exp(log_asymptote);
  dispersion = exp(log_dispersion);

  // calculate period effects
  {
  matrix[(total_years - 2) * 2, (total_years - 2) * 2] LK2;
  vector[(total_years - 2) * 2] delta_kk_raw;
  
  int i1;
  int i2;
  int i3;

  // index variables to help readability
  i1 = total_years - 2;
  i2 = total_years - 1;
  i3 = (total_years - 2) * 2;


  // cholesky covariance structure of conditional distributions
  LK2[1:i1, 1:i1] = LK * sigma_dk[1];
  LK2[i2:i3, 1:i1] = LK * rho_kk * sigma_dk[2];
  LK2[i2:i3, i2:i3] = LK * sigma_dk[2] * sqrt(1 - rho_kk^2);
  LK2[1:i1, i2:i3] = rep_matrix(rep_vector(0, i1), i1);
  // this VvV is equivalent to delta_kk ~ MVN(0, sigma.(LK'LK)), 
  // as dkz ~normal(0, 1);
  delta_kk_raw = LK2 * dkz;
  for (s in 1:n_sexes){
    int index[total_years - 2];
    for (i in 1:(total_years-2)){
      index[i] = (s - 1)*(total_years - 2) + i;
    }
    kk[s] =  cumulative_sum(k_inv_constraint * delta_kk_raw[index] * 0.1)';
  }
  }
}
model {
  // holding matrices for intermediate results ---------------------------------
  vector[cutoff - 1] smooth_mu[n_sexes];
  matrix[cutoff - 1, n_years] smooth_imp[n_sexes];
  matrix[cutoff - 1, n_years] mu_rep[n_sexes];
  matrix[cutoff - 1, n_years] kk_rep[n_sexes];
  row_vector[top_age + n_years + n_forecast_years] smooth_gam[n_sexes];
  matrix[top_age + 1, n_years] cohort_effects[n_sexes];

  // final rate matrices
  matrix[(top_age - cutoff + 1), n_years] old_base_rate[n_sexes];
  matrix[(top_age - cutoff + 1), n_years] old_logit_rate[n_sexes];
  matrix[(top_age - cutoff + 1), n_years] old_log_rate[n_sexes];
  matrix[cutoff, n_years] young_log_rate[n_sexes];  

    // precision matrices for prior on basis function coefficients.
  matrix[n_basis_age, n_basis_age] full_penalty_mu[n_sexes];
  matrix[n_basis_age, n_basis_age] full_penalty_aa[n_sexes];


  // priors --------------------------------------------------------------------
  
  // infants
  a0 ~ normal(0, 10);
  m0 ~ normal(0, 10);

  // old age
  m_old_0 ~ normal(0, 10);
  m_old_1 ~ normal(0, 10);
  a_old_0 ~ normal(0, 10);
  a_old_11 ~ normal(0, 10) T[- m_old_1[1] / max(time),]; // enforce positive gradient
  a_old_12 ~ normal(0, 10) T[- m_old_1[2] / max(time),]; // enforce positive gradient
  log_asymptote ~ normal(0, 1);

  // period priors -------------
  // Note: all sigma are half-normal, taking positive part only.
  sigma_dk ~ normal(0, 10);
  // this prior is used in combination with cholesky decomposition of 
  // precomputed covariance matrix to give the correct prior on kk
  dkz ~ normal(0, 1); 
  
  // cohort priors -------------
  // this prior is used in combination with cholesky decomposition of 
  // precomputed covariance matrix to give the correct prior on gam
  dgamz[1] ~ normal(0, 1);
  dgamz[2] ~ normal(0, 1);
  sigma_dgam ~ normal(0, 10);

  // basis function coefficient priors
  sigma_mu ~ normal(0, 10);
  sigma_mu_null ~ normal(0, 10);
  sigma_aa ~ normal(0, 10);
  sigma_aa_null ~ normal(0, 10);
  

  // implicit prior
  // log_dispersion ~ 1

  // penalties on rough vectors ------------------------------------------------
  for (s in 1:n_sexes){
    // two penalties - one on differences, one to give proper prior.
    full_penalty_mu[s] = ((1 / pow(sigma_mu[s], 2)) * penalty_matrix_age + 
                       (1 / pow(sigma_mu_null[s], 2)) * penalty_matrix_null); 
    // two penalties - one on differences, one to give proper prior.
    full_penalty_aa[s] = ((1 / pow(sigma_aa[s], 2)) * penalty_matrix_age + 
                         (1 / pow(sigma_aa_null[s], 2)) * penalty_matrix_null); 
      
    mu ~ multi_normal_prec(zeroes[1:n_basis_age], full_penalty_mu[s]);
    aa ~ multi_normal_prec(zeroes[1:n_basis_age], full_penalty_aa[s]);
  }

  // big loop over sexes--------------------------------------------------------
  for (s in 1:n_sexes){      
    // temporary rate variables
    real lambda;
    row_vector[n_years] lambda_v;
    real a_old_1;

    if (s==1){
      a_old_1 = a_old_11;
    } else {
      a_old_1 = a_old_12;
    }

    // year-specific effects -----------------------------------------------------
    kk_rep[s] = rep_matrix(kk[s][1:n_years], cutoff - 1);

    // static age-specific log-rate pattern mu -----------------------------------
    smooth_mu[s] = age_basis * mu[s]; 
    mu_rep[s] = rep_matrix(smooth_mu[s], n_years);

    // improvement rates ---------------------------------------------------------
    smooth_imp[s] = age_basis * aa[s] * time;


    // cohort --------------------------------------------------------------------
                
    {
    
    // The first few elements of the smooth cohort function are left blank.
    // These correspond to the top left hand corner of the age X period lexis 
    // surface, and do not enter the likelihood, but constructing the cohort 
    // effect matrix is made easier by adding zeroes here.
    int coh_i_1 = (first_cohort + ((top_age + 1) - N));
    int coh_i_2 = top_age + n_years + n_forecast_years;
    smooth_gam[s, coh_i_1:coh_i_2] = (cohort_basis * gam[s])';
    }
    for (i in 1:(first_cohort + (top_age - N))){
      smooth_gam[s][i] = 0;
    }  
    
    // fill in cohort matrix
    for (i in 1:(top_age + 1)){
      cohort_effects[s][i] = smooth_gam[s][(1 + (top_age + 1) - i):
                                     (n_years + (top_age + 1) - i)];
    }

    // data model ----------------------------------------------------------------

    // infants -------------------------------------------------------------------
    young_log_rate[s][1] = m0[s] + a0[s] * time + cohort_effects[s][1];
    
    // gam model -----------------------------------------------------------------
    young_log_rate[s, 2:cutoff] = (mu_rep[s] + smooth_imp[s] + kk_rep[s] + 
                                cohort_effects[s][2:cutoff]);

    // old age -------------------------------------------------------------------  
    old_base_rate[s] = ((m_old_0[s] + m_old_1[s] * old_age_index) * ones + 
                     (a_old_0[s] + a_old_1  * old_age_index) * time);

    old_logit_rate[s] = (old_base_rate[s]  + kk_rep[s][1:(top_age - cutoff + 1)] + 
                      cohort_effects[s][(cutoff + 1):(top_age + 1)]);

    // log likelihoods -----------------------------------------------------------

    // main gam-based model ------------------------------------------------------
    for (i in 1:cutoff){
      lambda_v = young_log_rate[s][i] + log_expos[s][i];
      deaths[s][i] ~ neg_binomial_2_log(lambda_v, dispersion[s]);
    }
    
    // old_age -------------------------------------------------------------------
    // calculate rates
    for (i in (cutoff + 1):(top_age + 1)){
      for (t in 1:n_years){
        old_log_rate[s][i - cutoff, t] = (old_logit_rate[s][i - cutoff, t] - 
                                 log1p_exp(old_base_rate[s][i - cutoff, t] -
                                          log_asymptote));
      }
    }

    // old-age likelihood
    for (i in (cutoff + 1):N){
      for (t in 1:(n_years)){
        // exclude early cohorts and points with zero exposure and
        if (expos[s][i, t] != 0 && N - i + t >= first_cohort){
          deaths[s][i, t]  ~ neg_binomial_2_log(log_expos[s][i, t] +
                                             old_log_rate[s][i - cutoff, t], 
                                             dispersion[s]);                    
        } 
      }
    }
  }
}
