
################################################################################
## utility function
################################################################################

create_area_grid <- function(data_val, # data.frame(lon,lat,type)
                             coords, # names of the coords: lon, lat
                             grid_x = 3, grid_y = 4,
                             crs = 4326) {
  ##
  options(warn=-1)
  names(data_val)[3] <- "Area"
  data_val$id_J <- 1:nrow(data_val)
  points_sf <- st_as_sf(
    x = data_val,
    coords = coords,
    crs = crs
  )
  ## bounding box
  bbox <- st_bbox(points_sf)
  ## create grid
  grid <- st_make_grid(
    st_as_sfc(bbox),
    n = c(grid_x, grid_y),
    what = "polygons"
  )
  grid_sf <- st_sf(
    id_m = 1:length(grid),
    geometry = grid
  )
  ## assign points to grid areas
  points_sf <- st_join(
    points_sf,
    grid_sf,
    join = st_within
  )
  ## assign nearest area for points outside grid
  na_idx <- which(is.na(points_sf$id_m))
  if (length(na_idx) > 0) {
    centroids <- st_centroid(grid_sf)
    nearest_idx <- st_nearest_feature(
      points_sf[na_idx, ],
      centroids
    )
    points_sf$id_m[na_idx] <-
      centroids$id_m[nearest_idx]
  }
  ## grid midpoints
  grid_midpoints <- st_centroid(grid_sf)
  ## plot
  p <- ggplot() +
    geom_sf(
      data = grid_sf,
      fill = NA,
      color = "black",
      linetype = "dashed"
    ) +
    geom_sf(
      data = points_sf,
      aes(color = Area, shape = Area),
      size = 2,
      alpha = 1
    ) +
    geom_sf(
      data = grid_midpoints,
      shape = 7,
      color = "black",
      size = 2
    ) +
    labs(color = "", shape = "") +
    theme_classic() +
    theme(legend.position = "bottom")
  ##
  return(list(
    data = points_sf,
    grid = grid_sf,
    midpoints = grid_midpoints,
    plot = p
  ))
  ##
}

## prepare data

prepare_stan_data <- function(f,
                              data,
                              id_var,
                              coord_var,
                              cluster_var,
                              urban_rural_var,
                              grid_x = 5, grid_y = 5,
                              sigma_u_prior = c(20, 20), #shape, rate
                              sigma_r_prior = c(10, 25), # shape, rate
                              sigma_nu_prior = c(2, 1), #shape, rate
                              beta_prior = c(0, 2), # mean, sd
                              zeta_prior = c(0, 2) # mean, sd
                              )
{
  ## check for z variable - currently can take only one z
  if(length(cluster_var)>2){stop("\n only take one cluster level variable \n")}
  ## create layer coordinates m & cluster coordinate points J
  layer_data <- create_area_grid(data_val = unique(data[,c(coord_var,urban_rural_var)]),
                   coords = coord_var,
                   grid_x = grid_x, grid_y = grid_y)
  coord_m <- st_coordinates(layer_data$midpoints)
  coord_J <- unique(data[,c(coord_var)])
  d <- geodist::geodist(coord_J, coord_m, measure = "geodesic") / 1000 # in KM
  ## function to create W matrix
  fnc_w <- function(cluster_J, layer_m) {
    W_Jxm <- table(cluster_J, layer_m)
    J <- nrow(W_Jxm)
    m <- ncol(W_Jxm)
    W <- matrix(0, nrow = J, ncol = J)
    for (i in 1:m) {
      ck <- W_Jxm
      dimnames(ck)[[1]] <- 1:dim(ck)[1]
      number <- as.numeric(names(ck[, i][ck[, i] == 1]))
      W[number, number] <- 1
    }
    return(W)
  }
  W <- fnc_w(cluster_J = layer_data$data$id_J,
             layer_m = layer_data$data$id_m)
  ## prepare cluster-level data
  z <- unique(data[, c(cluster_var,coord_var)]); z <- z[,1]
  ## IMAT matrix - initial
  Ind <- data.frame(clst = data[, id_var])
  for(i in layer_data$data$id_J){
    Ind[[paste0("clst", i)]] <- as.numeric(Ind$clst == i)
  }
  IMAT <- as.matrix(Ind[, -1]) # N x J
  ## urban-rural index
  ur <- unique(data[,c(urban_rural_var,id_var)])
  ur <- as.numeric(as.factor(ur[,1]))
  ## MW matrix - inital
  z <- as.matrix(z)
  MW <- (1 - z[,1] %*% solve(t(z[,1]) %*% z[,1]) %*% t(z[,1])) %*%
    W %*%
    (1 - z[,1] %*% solve(t(z[,1]) %*% z[,1]) %*% t(z[,1]))
  ## prepare QQ and eigen decomposition - initial
  x <- model.matrix(f, data)
  phi <- max(d)
  QQ <- t(MW %*% (1 - (d/phi)^2)^2) %*% (diag(W)*1 - W) %*% (MW %*% (1 - (d/phi)^2)^2)
  eigen_QQ <- eigen(QQ)
  ## initial zstar
  J <- nrow(d)
  M <- ncol(d)
  M_JxJ <- eigen(MW)$vectors
  D_JxM <- d
  Dstar <- d
  Sigma_diag <- abs(diag(eigen_QQ$vectors))
  sigma_r <- 1 / rgamma(100, shape = sigma_r_prior[1], rate = sigma_r_prior[2])
  sigma_u <- 1 / rgamma(100, shape = sigma_u_prior[1], rate = sigma_u_prior[2])
  for(m in 1:M){
    for(j in 1:J){
      if(ur[j] == 1){ # rural
        sigma_r_val <- sd(rnorm(length(sigma_r), mean = D_JxM[j, m], sd = sigma_r))
        Dstar[j, m] <- rnorm(1, D_JxM[j, m], sigma_r_val)
      } else { # urban
        sigma_u_val <- sd(rnorm(length(sigma_u), mean = D_JxM[j, m], sd = sigma_u))
        Dstar[j, m] <- rnorm(1, D_JxM[j, m], sigma_u_val)
      }
    }
  }
  Psi <- M_JxJ %*% (1 - (Dstar / phi)^2)^2
  Psi_eta <- Psi %*% Sigma_diag
  zstar <- z[,1] + Psi_eta
  ## scale z - initial
  z_scale <- scale(cbind(z[,1], zstar))
  #X_block <- matrix(0, N, J)
  #cluster_factor <- as.numeric(factor(df$DHSCLUST))
  #X_block[cbind(1:N, cluster_factor)] <- z
  ## data list
  data_list <- list(
    N = nrow(x),
    y = model.frame(f, data)[,1],
    K = ncol(x),
    Q = 1,
    J = J,
    M = M,
    x = x,
    z = as.matrix(z_scale[,1]),
    zstar = as.matrix(z_scale[,2]),
    ur = ur,
    D_JxM = D_JxM,
    Dstar = Dstar,
    IMAT = IMAT,
    MJJ = M_JxJ,
    Sigma_diag = as.matrix(Sigma_diag),
    phi = phi,
    sig_nu_shape = sigma_nu_prior[1],
    sig_nu_scale = 1/sigma_nu_prior[2],
    sig_u_shape = sigma_u_prior[1],
    sig_u_scale = 1/sigma_u_prior[2],
    sig_r_shape = sigma_r_prior[1],
    sig_r_scale = 1/sigma_r_prior[2],
    beta_mu = beta_prior[1],
    beta_sd = beta_prior[2],
    zeta_mu = zeta_prior[1],
    zeta_sd = zeta_prior[2]
  )
  return(list(data_list_stan = data_list,
         data = data.frame(clst = data[, id_var],
                           clst_var = data[, cluster_var],
                           model.frame(f, data)),
         sp_data = layer_data))
}

## summary function for binary outcome

evolve_summary_binary <- function(data, out, cluster_var,
                           digits = 2, sim = 10000,
                           cluster_summary = FALSE) {
  ##
  K <- data$data_list_stan$K
  Q <- data$data_list_stan$Q
  ## covariance matrix
  opt_cov <- solve(-out$hessian[1:(K + Q), 1:(K + Q)])
  fixed_simulations <- list()
  variability_simulations <- list()
  post.beta.mn <- out$par[1:K]
  post.beta.sd <- sqrt(diag(opt_cov[1:K, 1:K]))
  beta_summary <- matrix(NA, nrow = K, ncol = 3)
  for (i in 1:K) {
    sims <- exp(rnorm(
      sim,
      mean = post.beta.mn[i],
      sd = post.beta.sd[i]
    ))
    fixed_simulations[[colnames(data$data_list_stan$x)[i]]] <- sims
    beta_summary[i, ] <- quantile(
      sims,
      prob = c(0.5, 0.025, 0.975),
      na.rm = TRUE
    )
  }
  rownames(beta_summary) <- colnames(data$data_list_stan$x)
  colnames(beta_summary) <- c(
    "OR_median",
    "OR_lower_95",
    "OR_upper_95"
  )
  post.zeta.mn <- out$par[(K + 1):(K + Q)]
  post.zeta.sd <- sqrt(diag(opt_cov)[(K + 1):(K + Q)])
  zeta_summary <- matrix(NA, nrow = Q, ncol = 3)
  for (i in 1:Q) {
    sims <- exp(rnorm(
      sim,
      mean = post.zeta.mn[i],
      sd = post.zeta.sd[i]
    ))
    fixed_simulations[[cluster_var[i]]] <- sims
    zeta_summary[i, ] <- quantile(
      sims,
      prob = c(0.5, 0.025, 0.975),
      na.rm = TRUE
    )
  }
  rownames(zeta_summary) <- cluster_var[1:Q]
  colnames(zeta_summary) <- c(
    "OR_median",
    "OR_lower_95",
    "OR_upper_95"
  )
  fixed_summary <- rbind(
    beta_summary,
    zeta_summary
  )
  var_names <- c(
    "sigma_r_hat",
    "sigma_u_hat",
    "sigma_nu_hat"
  )
  cov_para <- solve(
    -out$hessian[var_names, var_names]
  )
  sig_means <- out$par[var_names]
  sig_sds <- sqrt(diag(cov_para))
  variability_summary <- matrix(
    NA,
    nrow = length(var_names),
    ncol = 3
  )
  for (i in seq_along(var_names)) {
    sims <- rnorm(
      sim,
      mean = sig_means[i],
      sd = sig_sds[i]
    )
    sims[sims < 0] <- 0
    variability_simulations[[var_names[i]]] <- sims
    variability_summary[i, ] <- quantile(
      sims,
      prob = c(0.5, 0.025, 0.975),
      na.rm = TRUE
    )
  }
  rownames(variability_summary) <- var_names
  colnames(variability_summary) <- c(
    "median",
    "lower_95",
    "upper_95"
  )
  ##
  df <- data$data
  df$pred_prob <- out$par[grep("y_prob", names(out$par))]
  ## cluster-level summary by clst
  cluster_summary <- df %>%
    group_by(clst) %>%
    dplyr::summarise(
      prob_mean = mean(pred_prob, na.rm = TRUE),
      prob_sd = sd(pred_prob, na.rm = TRUE),
      prob_median = median(pred_prob, na.rm = TRUE),
      prob_low  = quantile(pred_prob, 0.025, na.rm = TRUE),
      prob_up   = quantile(pred_prob, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  ## log-lilekihood, approximate waic
  log_lik_vec <- dbinom(df[,1+Q+1], size = 1,
                        prob = df$pred_prob, log = TRUE)
  lppd <- sum(log_lik_vec); waic_approx <- -2 * (lppd)
  ## extract zstart values to store
  zstar_hat <- out$par[grep("zstar_hat_scaled", names(out$par))]
  zstar_hat <- matrix(zstar_hat, nrow = nrow(df))
  zstar_hat <- data.frame(clst = df$clst, zstar_hat)
  names(zstar_hat)[-1] = c(cluster_var)
  zstar_hat <- tibble(zstar_hat)
  ##
  return(list(
    fixed_parameters = round(fixed_summary, digits),
    variability_parameters = round(variability_summary, digits),
    simulations = tibble(data.frame(fixed_simulations,variability_simulations)),
    cluster_summary =  cluster_summary,
    zstar_scaled = zstar_hat,
    lppd = lppd, waic_approx = waic_approx
  ))
  ##
}

## run function for binary

run_fnc_binary <- function(data,
                      cluster_var = c("Precipitation_2020"),
                      digits=4,
                      draws = 10){
  ##
  stan_file <- system.file(
    "stan",
    "svyLocAdj_logistic.stan",
    package = "svyLocAdj"
  )
  model <- rstan::stan_model(stan_file)
  out <- rstan::optimizing(model, data = data$data_list_stan,
                           hessian = TRUE, draws = draws,
                           importance_resampling = TRUE)
  para <- evolve_summary_binary(data = data, out = out,
                         cluster_var = cluster_var,
                         digits = digits)
  para
  ##
}

## for gaussian: continuous outcome

evolve_summary_gaussian <- function(data, out, cluster_var,
                                  digits = 2, sim = 10000,
                                  cluster_summary = FALSE) {
  ##
  K <- data$data_list_stan$K
  Q <- data$data_list_stan$Q
  ## covariance matrix
  opt_cov <- solve(-out$hessian[1:(K + Q), 1:(K + Q)])
  fixed_simulations <- list()
  variability_simulations <- list()
  post.beta.mn <- out$par[1:K]
  post.beta.sd <- sqrt(diag(opt_cov[1:K, 1:K]))
  beta_summary <- matrix(NA, nrow = K, ncol = 3)
  for (i in 1:K) {
    sims <- rnorm(
      sim,
      mean = post.beta.mn[i],
      sd = post.beta.sd[i]
    )
    fixed_simulations[[colnames(data$data_list_stan$x)[i]]] <- sims
    beta_summary[i, ] <- quantile(
      sims,
      prob = c(0.5, 0.025, 0.975),
      na.rm = TRUE
    )
  }
  rownames(beta_summary) <- colnames(data$data_list_stan$x)
  colnames(beta_summary) <- c(
    "median",
    "lower_95",
    "upper_95"
  )
  post.zeta.mn <- out$par[(K + 1):(K + Q)]
  post.zeta.sd <- sqrt(diag(opt_cov)[(K + 1):(K + Q)])
  zeta_summary <- matrix(NA, nrow = Q, ncol = 3)
  for (i in 1:Q) {
    sims <- rnorm(
      sim,
      mean = post.zeta.mn[i],
      sd = post.zeta.sd[i]
    )
    fixed_simulations[[cluster_var[i]]] <- sims
    zeta_summary[i, ] <- quantile(
      sims,
      prob = c(0.5, 0.025, 0.975),
      na.rm = TRUE
    )
  }
  rownames(zeta_summary) <- cluster_var[1:Q]
  colnames(zeta_summary) <- c(
    "median",
    "lower_95",
    "upper_95"
  )
  fixed_summary <- rbind(
    beta_summary,
    zeta_summary
  )
  var_names <- c(
    "sigma_r_hat",
    "sigma_u_hat",
    "sigma_nu_hat",
    "sigma_y"
  )
  cov_para <- solve(
    -out$hessian[var_names, var_names]
  )
  sig_means <- out$par[var_names]
  sig_sds <- sqrt(diag(cov_para))
  variability_summary <- matrix(
    NA,
    nrow = length(var_names),
    ncol = 3
  )
  for (i in seq_along(var_names)) {
    sims <- rnorm(
      sim,
      mean = sig_means[i],
      sd = sig_sds[i]
    )
    sims[sims < 0] <- 0
    variability_simulations[[var_names[i]]] <- sims
    variability_summary[i, ] <- quantile(
      sims,
      prob = c(0.5, 0.025, 0.975),
      na.rm = TRUE
    )
  }
  rownames(variability_summary) <- var_names
  colnames(variability_summary) <- c(
    "median",
    "lower_95",
    "upper_95"
  )
  ##
  df <- data$data
  df$pred <- out$par[grep("y_pred", names(out$par))]
  ## cluster-level summary by clst
  cluster_summary <- df %>%
    group_by(clst) %>%
    dplyr::summarise(
      pred_mean = mean(pred, na.rm = TRUE),
      pred_low  = quantile(pred, 0.025, na.rm = TRUE),
      pred_up   = quantile(pred, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  ##
  log_lik_vec <- dnorm(df[,1+Q+1], mean = df$pred,
                       sd = variability_summary[4,1],log = TRUE)
  lppd <- sum(log_lik_vec); waic_approx <- -2 * (lppd)
  ## extract zstart values to store
  zstar_hat <- out$par[grep("zstar_hat_scaled", names(out$par))]
  zstar_hat <- matrix(zstar_hat, nrow = nrow(df))
  zstar_hat <- data.frame(clst = df$clst, zstar_hat)
  names(zstar_hat)[-1] = c(cluster_var)
  zstar_hat <- tibble(zstar_hat)
  ##
  return(list(
    fixed_parameters = round(fixed_summary, digits),
    variability_parameters = round(variability_summary, digits),
    #simulations = as_draws_df(data.frame(fixed_simulations,variability_simulations)),
    simulations = tibble(data.frame(fixed_simulations,variability_simulations)),
    cluster_summary =  cluster_summary,
    zstar_scaled = zstar_hat,
    lppd = lppd, waic_approx = waic_approx
  ))
  ##
}

## run function for cont. outcome

run_fnc_gaussian <- function(data,
                           cluster_var = c("Precipitation_2020"),
                           digits=4,
                           draws = 10){
  ##
  stan_file <- system.file(
    "stan",
    "svyLocAdj_gaussian.stan",
    package = "svyLocAdj"
  )
  model <- rstan::stan_model(stan_file)
  out <- rstan::optimizing(model, data = data$data_list_stan,
                           hessian = TRUE, draws = draws,
                           importance_resampling = TRUE)
  para <- evolve_summary_gaussian(data = data, out = out,
                                cluster_var = cluster_var,
                                digits = digits)
  para
  ##
}

## summary function for poisson: count outcome

evolve_summary_poisson <- function(data, out, cluster_var,
                                   digits = 2,
                                   sim = 10000,
                                   cluster_summary = FALSE) {
  ##
  K <- data$data_list_stan$K
  Q <- data$data_list_stan$Q
  opt_cov <- solve(-out$hessian[1:(K + Q), 1:(K + Q)])
  fixed_simulations <- list()
  variability_simulations <- list()
  post.beta.mn <- out$par[1:K]
  post.beta.sd <- sqrt(diag(opt_cov[1:K, 1:K]))
  beta_summary <- matrix(NA, nrow = K, ncol = 3)
  for (i in 1:K) {
    sims <- exp(
      rnorm(
        sim,
        mean = post.beta.mn[i],
        sd = post.beta.sd[i]
      )
    )
    fixed_simulations[[colnames(data$data_list_stan$x)[i]]] <- sims
    beta_summary[i, ] <- quantile(
      sims,
      prob = c(0.5, 0.025, 0.975),
      na.rm = TRUE
    )
  }
  rownames(beta_summary) <- colnames(data$data_list_stan$x)
  colnames(beta_summary) <- c(
    "IRR_median",
    "IRR_lower_95",
    "IRR_upper_95"
  )
  post.zeta.mn <- out$par[(K + 1):(K + Q)]
  post.zeta.sd <- sqrt(
    diag(opt_cov)[(K + 1):(K + Q)]
  )
  zeta_summary <- matrix(NA, nrow = Q, ncol = 3)
  for (i in 1:Q) {
    sims <- exp(
      rnorm(
        sim,
        mean = post.zeta.mn[i],
        sd = post.zeta.sd[i]
      )
    )
    fixed_simulations[[cluster_var[i]]] <- sims
    zeta_summary[i, ] <- quantile(
      sims,
      prob = c(0.5, 0.025, 0.975),
      na.rm = TRUE
    )
  }
  rownames(zeta_summary) <- cluster_var[1:Q]
  colnames(zeta_summary) <- c(
    "IRR_median",
    "IRR_lower_95",
    "IRR_upper_95"
  )
  fixed_summary <- rbind(
    beta_summary,
    zeta_summary
  )
  var_names <- c(
    "sigma_r_hat",
    "sigma_u_hat",
    "sigma_nu_hat"
  )
  cov_para <- solve(
    -out$hessian[var_names, var_names]
  )
  sig_means <- out$par[var_names]
  sig_sds <- sqrt(diag(cov_para))
  variability_summary <- matrix(
    NA,
    nrow = length(var_names),
    ncol = 3
  )
  for (i in seq_along(var_names)) {
    sims <- rnorm(
      sim,
      mean = sig_means[i],
      sd = sig_sds[i]
    )
    sims[sims < 0] <- 0
    variability_simulations[[var_names[i]]] <- sims
    variability_summary[i, ] <- quantile(
      sims,
      prob = c(0.5, 0.025, 0.975),
      na.rm = TRUE
    )
  }
  rownames(variability_summary) <- var_names
  colnames(variability_summary) <- c(
    "median",
    "lower_95",
    "upper_95"
  )
  df <- data$data
  df$lambda <- out$par[
    grep("lambda", names(out$par))
  ]
  cluster_summary <- df %>%
    group_by(clst) %>%
    dplyr::summarise(
      lambda_mean   = mean(lambda, na.rm = TRUE),
      lambda_sd     = sd(lambda, na.rm = TRUE),
      lambda_median = median(lambda, na.rm = TRUE),
      lambda_low    = quantile(lambda, 0.025, na.rm = TRUE),
      lambda_up     = quantile(lambda, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  ## assumes outcome variable is column (Q + 2)
  y_obs <- df[, 1 + Q + 1]
  log_lik_vec <- dpois(
    y_obs,
    lambda = df$lambda,
    log = TRUE
  )
  lppd <- sum(log_lik_vec)
  waic_approx <- -2 * lppd
  zstar_hat <- out$par[grep("zstar_hat_scaled", names(out$par))]
  zstar_hat <- matrix(zstar_hat,nrow = nrow(df))
  zstar_hat <- data.frame(
    clst = df$clst,
    zstar_hat
  )
  names(zstar_hat)[-1] <- cluster_var
  zstar_hat <- tibble(zstar_hat)
  return(
    list(fixed_parameters = round(fixed_summary, digits),
         variability_parameters = round(variability_summary, digits),
         simulations = tibble(data.frame(fixed_simulations,variability_simulations)),
         cluster_summary = cluster_summary,
         zstar_scaled = zstar_hat,
         lppd =  lppd, waic_approx = waic_approx
    )
  )
}

## run function for poisson model

run_fnc_poisson <- function(data,
                           cluster_var = c("Precipitation_2020"),
                           digits=4,
                           draws = 10){
  ##
  stan_file <- system.file(
    "stan",
    "svyLocAdj_poisson.stan",
    package = "svyLocAdj"
  )
  model <- rstan::stan_model(stan_file)
  out <- rstan::optimizing(model, data = data$data_list_stan,
                           hessian = TRUE, draws = draws,
                           importance_resampling = TRUE)
  para <- evolve_summary_poisson(data = data, out = out,
                                cluster_var = cluster_var,
                                digits = digits)
  para
  ##
}

################################################################################
################################################################################





