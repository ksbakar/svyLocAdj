data {
  int<lower=0> N;                 // number of observations
  vector[N] y;                    // continuous response
  int<lower=0> K;                 // number of predictors
  int<lower=0> Q;                 // number of cluster-level predictors, it should be Q=1
  int<lower=0> J;                 // number of clusters
  int<lower=0> M;                 // number of 2nd stage areal mid-points
  matrix[N, K] x;                 // predictor matrix
  matrix[J, Q] z;                 // cluster/spatial predictors
  matrix[J, Q] zstar;
  array[J] int<lower=1, upper=2> ur;
  matrix[J, M] D_JxM;
  matrix[J, M] Dstar;
  matrix[J, J] D_JxJ; // for spatially varying model
  array[N] int<lower=1,upper=J> clst;     // cluster index for spatially varying model
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
  matrix[J] zeta;

  real<lower=0> sigma_zeta;
  real<lower=0> rho_zeta;

  real<lower=0> sigma_nu_hat;
  real<lower=0> sigma_u_hat;
  real<lower=0> sigma_r_hat;
  real<lower=0> sigma_y;          // residual SD
}

model {
  matrix[J, M] Dstar_phi;
  matrix[J, M] Psi;
  matrix[J, Q] Psi_eta;
  matrix[J] zstar_hat;
  matrix[N] zzstar;
  vector[N] eta;

  // Priors

  sigma_u_hat ~ inv_gamma(sig_u_shape, sig_u_scale);
  sigma_r_hat ~ inv_gamma(sig_r_shape, sig_r_scale);
  sigma_nu_hat ~ inv_gamma(sig_nu_shape, sig_nu_scale);
  sigma_y ~ inv_gamma(sig_nu_shape, sig_nu_scale);

  beta ~ normal(beta_mu, beta_sd);

  // spatial prior

  sigma_zeta ~ inv_gamma(sig_nu_shape, sig_nu_scale);
  matrix[J,J] Sigma_zeta;
  for (j1 in 1:J) {
    for (j2 in 1:J) {
      Sigma_zeta[j1,j2] =  square(sigma_zeta) * exp(-D_JxJ[j1,j2] / rho_zeta);
    }
  }
  zeta ~ multi_normal(rep_vector(zeta_mu, J),Sigma_zeta);

  // Spatial process

  for (m in 1:M) {
    for (j in 1:J) {
      if (ur[j] == 1) {
        Dstar[j, m] ~ normal(D_JxM[j, m], sigma_r_hat);
      } else {
        Dstar[j, m] ~ normal(D_JxM[j, m], sigma_u_hat);
      }
      Dstar_phi[j, m] = square(1 - square(Dstar[j, m] / phi));
    }
  }

  Psi = MJJ * Dstar_phi;
  Psi_eta = Psi * Sigma_diag;

  for (j in 1:J) {
      zstar_hat[j] = z[j, 1] + Psi_eta[j, 1];
      z[j, 1] ~ normal(zstar[j, 1], sigma_nu_hat);
  }
  zzstar = IMAT * zstar_hat; // N x Q = (NxJ)x(JxQ)
  //
  for (i in 1:N) {
    eta[i] = dot_product(x[i], beta) + zzstar[clst[i]] * zeta;
  }
  y ~ normal(eta, sigma_y);
}

generated quantities {
  vector[N] y_mean;
  vector[N] y_pred;
  vector[J] z_effect;
  vector[N] zstar_hat_scaled;
  for (j in 1:J){
    z_effect[j] = zstar[j,1] * zeta[j];
  }
  zstar_hat_scaled = IMAT * z_effect;
  for (i in 1:N) {
    y_mean[i] = dot_product(x[i], beta) + zstar_hat_scaled[i];
    y_pred[i] = normal_rng(y_mean[i], sigma_y);
  }
}
