data {
  int<lower=1> N;                       // number of observations
  vector[N] y;                          // continuous outcome
  int<lower=1> K;                       // number of individual predictors
  int<lower=1> Q;                       // number of cluster predictors
  int<lower=1> J;                       // number of clusters
  int<lower=1> M;                       // number of areal mid-points
  matrix[N, K] x;                       // individual-level predictors
  matrix[J, Q] z;                       // observed cluster predictors
  matrix[J, Q] zstar;                   // latent/true cluster predictors
  array[J] int<lower=1, upper=2> ur;    // urban-rural indicator
  matrix[J, M] D_JxM;
  matrix[J, M] Dstar;
  matrix[N, J] IMAT;                    // observation x cluster incidence matrix
  matrix[J, J] MJJ;
  matrix[M, Q] Sigma_diag;
  real phi;
  real<lower=0> sig_nu_shape;
  real<lower=0> sig_nu_scale;
  real<lower=0> sig_u_shape;
  real<lower=0> sig_u_scale;
  real<lower=0> sig_r_shape;
  real<lower=0> sig_r_scale;
  real beta_mu;
  real beta_sd;
  real zeta_mu;
  vector[J] zeta_sd; // this is for car prior
}
parameters {
  vector[K] beta;
  // J cluster-specific coefficients
  // for each of Q cluster-level covariates
  matrix[J, Q] zeta;
  real<lower=0> sigma_nu_hat;
  real<lower=0> sigma_u_hat;
  real<lower=0> sigma_r_hat;
  real<lower=0> sigma_y_hat;
  real<lower=0> sigma_zeta_hat;
}
transformed parameters {
  matrix[J, M] Dstar_phi;
  matrix[J, M] Psi;
  matrix[J, Q] Psi_eta;
  matrix[J, Q] zstar_hat;
  // Cluster-level contribution
  vector[J] eta_z;
  // Observation-level cluster contribution
  matrix[N, Q] zzstar;
  // Cluster-specific zeta mapped to observations
  matrix[N, Q] zeta_obs;
  vector[N] eta;

  for (j in 1:J) {
    for (m in 1:M) {
      if (Dstar[j,m] < phi) {
        Dstar_phi[j,m] = square(1 - square(Dstar[j,m]/phi));
      } else {
        Dstar_phi[j,m] = 0;
      }
    }
  }
  Psi = MJJ * Dstar_phi;
  Psi_eta = Psi * Sigma_diag;
  for (j in 1:J) {
    for (q in 1:Q) {
      zstar_hat[j,q] = z[j,q] + Psi_eta[j,q];
    }
  }
  // IMAT is N x J
  // zstar_hat is J x Q
  // result is N x Q
  zzstar = IMAT * zstar_hat;
  // zeta is J x Q
  // IMAT * zeta gives N x Q
  // Each observation receives the coefficients
  // belonging to its cluster.
  zeta_obs = IMAT * zeta;
  eta_z = rep_vector(0, J);
  for (q in 1:Q) {
    eta_z +=  zstar_hat[,q] .* zeta[,q];
  }
  eta = x * beta + rows_dot_product(zzstar, zeta_obs);
}
model {
  sigma_u_hat ~ inv_gamma(sig_u_shape,sig_u_scale);
  sigma_r_hat ~ inv_gamma(sig_r_shape,sig_r_scale);
  sigma_nu_hat ~ inv_gamma(sig_nu_shape,sig_nu_scale);
  sigma_y_hat ~ inv_gamma(sig_nu_shape,sig_nu_scale);
  sigma_zeta_hat ~ inv_gamma(sig_nu_shape,sig_nu_scale);
  beta ~ normal(beta_mu,beta_sd);
  for (q in 1:Q){
    for (j in 1:J){
      zeta[j,q] ~ normal(zeta_mu, zeta_sd[j]);
    }
  }
  for (j in 1:J) {
    for (m in 1:M) {
      if (ur[j] == 1) {
        Dstar[j,m] ~ normal(D_JxM[j,m], sigma_u_hat);
      } else {
        Dstar[j,m] ~ normal(D_JxM[j,m], sigma_r_hat);
      }
    }
  }
  for (j in 1:J) {
    for (q in 1:Q) {
      z[j,q] ~ normal(zstar[j,q], sigma_nu_hat);
    }
  }
  y ~ normal(eta, sigma_y_hat);
}
generated quantities {
  vector[J] z_effect;
  vector[N] zstar_hat_scaled;
  vector[N] y_mean;
  vector[N] y_pred;
  for (j in 1:J) {
    z_effect[j] = 0;
    for (q in 1:Q) {
      z_effect[j] += zstar_hat[j,q] * zeta[j,q];
    }
  }
  zstar_hat_scaled = IMAT * z_effect;
  y_mean = x * beta + zstar_hat_scaled;
  for (i in 1:N) {
    y_pred[i] = normal_rng(y_mean[i],sigma_y_hat);
  }
}
