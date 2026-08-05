data {
  int<lower=0> N;                  // number of observations
  vector[N] y;                      // continuous outcome
  int<lower=0> K;                  // number of individual predictors
  int<lower=0> Q;                  // number of cluster predictors
  int<lower=0> J;                  // number of clusters
  int<lower=0> M;                  // number of areal mid-points
  matrix[N,K] x;                   // individual-level predictors
  matrix[J,Q] z;                   // observed cluster predictors
  matrix[J,Q] zstar;               // latent true cluster predictors
  array[J] int<lower=1,upper=2> ur; // urban-rural indicator
  matrix[J,M] D_JxM;
  matrix[J,M] Dstar;
  matrix[N,J] IMAT;
  matrix[J,J] MJJ;
  matrix[M,Q] Sigma_diag;
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
  real zeta_sd;
}
parameters {
  vector[K] beta;
  vector[Q] zeta;
  real<lower=0> sigma_nu_hat;
  real<lower=0> sigma_u_hat;
  real<lower=0> sigma_r_hat;
  real<lower=0> sigma_y;          // residual SD
}
model {
  matrix[J,M] Dstar_phi;
  matrix[J,M] Psi;
  matrix[J,Q] Psi_eta;
  matrix[J,Q] zstar_hat;
  vector[N] eta;
  
  sigma_u_hat ~ inv_gamma(sig_u_shape, sig_u_scale);
  sigma_r_hat ~ inv_gamma(sig_r_shape, sig_r_scale);
  sigma_nu_hat ~ inv_gamma(sig_nu_shape, sig_nu_shape);
  sigma_y ~ inv_gamma(sig_nu_shape, sig_nu_scale);

  beta ~ normal(beta_mu, beta_sd);
  zeta ~ normal(zeta_mu, zeta_sd);

  for(m in 1:M){
    for(j in 1:J){
      if(ur[j] == 1){       // rural
        Dstar[j,m] ~ normal(D_JxM[j,m], sigma_r_hat);
      } else {              // urban
        Dstar[j,m] ~ normal(D_JxM[j,m], sigma_u_hat);
      }
      Dstar_phi[j,m] = square(1 - square(Dstar[j,m] / phi));
    }
  }
  
  Psi = MJJ * Dstar_phi;
  Psi_eta = Psi * Sigma_diag;

  for(q in 1:Q){
    for(j in 1:J){
      zstar_hat[j,q] =
        z[j,q] + Psi_eta[j,q];
      z[j,q] ~ normal(zstar[j,q], sigma_nu_hat);
    }
  }
  eta = x * beta + IMAT * (zstar_hat * zeta);
  y ~ normal(eta, sigma_y);
}
generated quantities {
  vector[N] y_mean;
  vector[N] y_pred;
  vector[J] z_effect;
  vector[N] zstar_hat_scaled;

  z_effect = zstar * zeta;
  zstar_hat_scaled =
    IMAT * z_effect;
  y_mean =
    x * beta + zstar_hat_scaled;
  for(i in 1:N){
    y_pred[i] =
      normal_rng(y_mean[i], sigma_y);
  }
}