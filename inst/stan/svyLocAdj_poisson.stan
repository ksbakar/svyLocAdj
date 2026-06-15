data {
  int<lower=0> N;
  array[N] int<lower=0> y;   // count outcome
  int<lower=0> K;
  int<lower=0> Q;
  int<lower=0> J;
  int<lower=0> M;
  matrix[N, K] x;
  matrix[J, Q] z;
  matrix[J, Q] zstar;
  array[J] int<lower=1, upper=2> ur;
  matrix[J, M] D_JxM;
  matrix[J, M] Dstar;
  matrix[N, J] IMAT;
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
  real zeta_sd;
}

parameters {
  vector[K] beta;
  vector[Q] zeta;
  real<lower=0> sigma_nu_hat;
  real<lower=0> sigma_u_hat;
  real<lower=0> sigma_r_hat;
}

transformed parameters {
  // incidence rate ratios (optional)
  vector[K] IRR_beta;
  vector[Q] IRR_zeta;
  IRR_beta = exp(beta);
  IRR_zeta = exp(zeta);
}
model {
  sigma_u_hat ~ inv_gamma(sig_u_shape, sig_u_scale);
  sigma_r_hat ~ inv_gamma(sig_r_shape, sig_r_scale);
  sigma_nu_hat ~ inv_gamma(sig_nu_shape, sig_nu_scale);
  matrix[J, M] Dstar_phi;
  for (m in 1:M) {
    for (j in 1:J) {
      if (ur[j] == 1) {
        Dstar[j, m] ~ normal(D_JxM[j, m], sigma_r_hat);
      } else {
        Dstar[j, m] ~ normal(D_JxM[j, m], sigma_u_hat);
      }
      Dstar_phi[j, m] =
        square(1 - square(Dstar[j, m] / phi));
    }
  }
  matrix[J, M] Psi;
  matrix[J, Q] Psi_eta;
  matrix[J, Q] zstar_hat;
  Psi = MJJ * Dstar_phi;
  Psi_eta = Psi * Sigma_diag;
  for (q in 1:Q) {
    for (j in 1:J) {
      zstar_hat[j, q] =
        z[j, q] + Psi_eta[j, q];
      z[j, q] ~ normal(zstar[j, q], sigma_nu_hat);
    }
  }
  matrix[N, Q] zzstar;
  zzstar = IMAT * zstar_hat;
  // Poisson likelihood
  for (i in 1:N) {
    real eta;
    eta = dot_product(x[i], beta) + dot_product(zzstar[i], zeta);
    y[i] ~ poisson_log(eta);
  }
  beta ~ normal(beta_mu, beta_sd);
  zeta ~ normal(zeta_mu, zeta_sd);
}
generated quantities {
  vector[N] lambda;
  array[N] int y_pred;
  vector[J] z_effect;
  vector[N] zstar_hat_scaled;
  z_effect = zstar * zeta;
  zstar_hat_scaled = IMAT * z_effect;
  for (i in 1:N) {
    real eta;
    eta =
      dot_product(x[i], beta) +
      zstar_hat_scaled[i];
    // expected count
    lambda[i] = exp(eta);
    // posterior predictive counts
    y_pred[i] = poisson_log_rng(eta);
  }
}
