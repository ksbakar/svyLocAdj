  data {
    //
      int <lower = 0> N; // Defining the number of data points
      array[N] int<lower=0, upper=1> y; // A variable that describes [1] or [0]
      int <lower = 0> K;   // number of predictors
      int <lower = 0> Q;   // number of cluster level predictors
      int <lower = 0> J;   // number of clusters
      int <lower = 0> M;   // number of 2nd stage areal mid-points
      matrix[N, K] x;   // K predictor matrix
      matrix[J, Q] z;   // Q predictors at J cluster/spatial scale
      matrix[J, Q] zstar;   // Q predictors at J cluster/spatial scale
      array[J] int<lower=1, upper=2> ur; // urban-rural ID; 1=urban, 2=rural
      matrix[J, M] D_JxM; // distance matrix
      matrix[J, M] Dstar; // distance matrix
      matrix[N, J] IMAT; // matrix
      matrix[J, J] MJJ; // matrix
      matrix[M, Q] Sigma_diag; // matrix
      real phi; // decay
      real <lower = 0> sig_nu_shape; // prior
      real <lower = 0> sig_nu_scale; // prior: scale = 1/rate
      real <lower = 0> sig_u_shape; // prior
      real <lower = 0> sig_u_scale; // prior: scale = 1/rate
      real <lower = 0> sig_r_shape; // prior
      real <lower = 0> sig_r_scale; // prior: scale = 1/rate
      real beta_mu; // prior
      real beta_sd; // prior
      real zeta_mu; // prior
      real zeta_sd; // prior
      //
  }
  parameters {
    vector[K] beta;   // coefficients for predictors (including intercept, if applicable)
    vector[Q] zeta;   // coefficients for spatial predictors (including intercept, if applicable)
    real <lower = 0> sigma_nu_hat; // for zstar - true cluster level values
    real <lower = 0> sigma_u_hat; //  initiate - urban
    real <lower = 0> sigma_r_hat; //  initiate - rural
  }
  transformed parameters {
    vector[K] OR_beta;   // OR of coefficients for predictors (including intercept, if applicable)
    vector[Q] OR_zeta;   // OR of coefficients for spatial predictors (including intercept, if applicable)
    OR_beta = exp(beta);
    OR_zeta = exp(zeta);
  }
  model {
    //
    sigma_u_hat ~ inv_gamma(sig_u_shape, sig_u_scale);
    sigma_r_hat ~ inv_gamma(sig_r_shape, sig_r_scale);
    sigma_nu_hat ~ inv_gamma(sig_nu_shape, sig_nu_shape);
    //
    matrix[J, M] Dstar_phi; // matrix
    //
      for(m in 1:M){
        for(j in 1:J){
          if(ur[j] == 1){ // rural
            Dstar[j, m] ~ normal(D_JxM[j, m], sigma_r_hat);
          }
          else{ // urban
            Dstar[j, m] ~ normal(D_JxM[j, m], sigma_u_hat);
          }
          Dstar_phi[j, m] = (1-(Dstar[j, m]/phi)^2)^2;
        }
      }
    //
    matrix[J, M] Psi;
    matrix[J, Q] Psi_eta;
    matrix[J, Q] zstar_hat;
    Psi = MJJ * Dstar_phi; // (JxJ) x (JxM) = JxM
    Psi_eta = Psi * Sigma_diag; // (JxM) x (MxQ) = JxQ
    //
      for(q in 1:Q){
        for(j in 1:J){
          zstar_hat[j, q] = z[j, q] + Psi_eta[j, q]; // working
          z[j, q] ~ normal(zstar[j, q], sigma_nu_hat); // working
        }
      }
    //
    matrix[N, Q] zzstar;
    zzstar = IMAT * zstar_hat; // (NxJ) x (JxQ) = NxQ // IMAT => 0 and 1 design matrix
    //
    for(i in 1:N){
        y[i] ~ bernoulli_logit(x[i,] * beta + zzstar[i,] * zeta); // likelihood
    }
    for(k in 1:K){
        beta[k] ~ normal(beta_mu, beta_sd);
    }
    for(q in 1:Q){
      zeta[q] ~ normal(zeta_mu, zeta_sd);
    }
  }
  generated quantities {
    vector[N] y_prob;
    int<lower=0, upper=1> y_pred[N];
    //vector[N] log_lik;
    vector[J] z_effect;
    vector[N] zstar_hat_scaled;
    z_effect = zstar * zeta;          // Jxq * q => J
    zstar_hat_scaled = IMAT * z_effect;     // NxJ * J => N
    for (i in 1:N) {
      real eta = dot_product(x[i], beta) + zstar_hat_scaled[i];
      y_prob[i] = inv_logit(eta);
      y_pred[i] = bernoulli_rng(y_prob[i]);
    }
  }
