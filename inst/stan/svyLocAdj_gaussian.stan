data {
  int<lower=0> N;                 // number of observations
  vector[N] y;                    // continuous response
  int<lower=0> K;                 // number of predictors
  int<lower=0> Q;                 // number of cluster-level predictors
  int<lower=0> J;                 // number of clusters
  int<lower=0> M;                 // number of 2nd stage areal mid-points
  matrix[N, K] x;                 // predictor matrix
  matrix[J, Q] z;                 // cluster/spatial predictors
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


transformed data {
  int Jr = 0;
  int Ju = 0;

  for (j in 1:J) {
    if (ur[j] == 1)
      Jr += 1;
    else
      Ju += 1;
  }

  array[Jr] int idx_r;
  array[Ju] int idx_u;

  {
    int r = 1;
    int u = 1;

    for (j in 1:J) {
      if (ur[j] == 1) {
        idx_r[r] = j;
        r += 1;
      } else {
        idx_u[u] = j;
        u += 1;
      }
    }
  }
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
  matrix[J, M] Dstar_phi;
  matrix[J, M] Psi;
  matrix[J, Q] Psi_eta;
  matrix[J, Q] zstar_hat;
  vector[N] eta;

  // Priors
  sigma_u_hat ~ inv_gamma(sig_u_shape, sig_u_scale);
  sigma_r_hat ~ inv_gamma(sig_r_shape, sig_r_scale);
  sigma_nu_hat ~ inv_gamma(sig_nu_shape, sig_nu_scale);
  sigma_y ~ inv_gamma(sig_nu_shape, sig_nu_scale);

  beta ~ normal(beta_mu, beta_sd);
  zeta ~ normal(zeta_mu, zeta_sd);

  // Spatial process
  Dstar[idx_r, ] ~ normal(D_JxM[idx_r, ], sigma_r_hat);
  Dstar[idx_u, ] ~ normal(D_JxM[idx_u, ], sigma_u_hat);
  Dstar_phi =
    square(
      rep_matrix(1.0, J, M)
      - square(Dstar / phi)
    );

  //Psi = MJJ * Dstar_phi;
  //Psi_eta = Psi * Sigma_diag;
  Psi_eta = MJJ * (Dstar_phi * Sigma_diag);

  z ~ normal(zstar, sigma_nu_hat);
  zstar_hat = z + Psi_eta;

  // Linear predictor
  eta = x * beta + IMAT * (zstar_hat * zeta);
  // Gaussian likelihood
  y ~ normal(eta, sigma_y);
}

generated quantities {
  vector[N] y_mean;
  vector[N] y_pred;
  vector[J] z_effect;
  vector[N] zstar_hat_scaled;
  z_effect = zstar * zeta;
  zstar_hat_scaled = IMAT * z_effect;
  y_mean = x * beta + zstar_hat_scaled;
  for (i in 1:N)
     y_pred[i] = normal_rng(y_mean[i], sigma_y);
}
