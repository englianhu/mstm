data {
  int<lower = 1> n;
  int<lower = 1> T;
  int<lower = 1> L; // number of dimensions
  int<lower = 1> p; // number of columns in design matrix
  matrix[n * L, p] X[T];
  int<lower = 1> r; 
  matrix[n * L, r] S[T];
  matrix[r, r] M[T];
  matrix[n * L, T] Z;
}

parameters {
  vector[p] beta;
  matrix[r, T] etaR;
  cholesky_factor_corr[r] L_eta;
  vector<lower = 0>[r] sigma_0;
  vector<lower = 0>[r] sigma_eta;
  real<lower = 0> sigma_z;
}

transformed parameters {
  matrix[r, T] eta;
  matrix[n * L, T] mu;
  
  eta[, 1] = diag_pre_multiply(sigma_0, L_eta) *  etaR[, 1];
  for (t in 2:T) {
    eta[, t] = M[t] * eta[, t - 1] + 
    diag_pre_multiply(sigma_eta, L_eta) *  etaR[, t];
  }
  
  for (t in 1:T) {
    mu[, t] = X[t] * beta + S[t] * eta[, t];
  }
}

model {
  beta ~ normal(0, 1);
  to_vector(etaR) ~ normal(0, 1);
  L_eta ~ lkj_corr_cholesky(2);
  sigma_0 ~ normal(0, 1);
  sigma_eta ~ normal(0, 1);
  sigma_z ~ normal(0, 1);
  for (t in 1:T) Z[t] ~ normal(mu[t], sigma_z);
}

generated quantities {
  matrix[r, r] R_eta;
  
  R_eta = multiply_lower_tri_self_transpose(L_eta);
}
