
################################################################################
## run function
################################################################################

svymodel <- function(formula = as.formula(ch_stunt ~ age_of_child),
                   data = y,
                   family = "binomial", #gaussian, #poisson
                   id_var = c("DHSCLUST"),
                   coord_var = c("LONGNUM","LATNUM"),
                   cluster_var = c("Precipitation_2020"),
                   urban_rural_var = c("URBAN_RURA"),
                   grid_x = 5, grid_y = 5,
                   sigma_u_prior = c(20, 20), #shape, rate
                   sigma_r_prior = c(10, 25), # shape, rate
                   sigma_nu_prior = c(2, 1), #shape, rate
                   beta_prior = c(0, 2), # mean, sd
                   zeta_prior = c(0, 2), # mean, sd
                   digits=4, draws = 10){
  ##
  if(family%in%"binomial"){
    cat(paste0("Bayesian model for binary outcome ...\n"))
    cat("Running Maximum A Posteriori (MAP) estimate ...\n")
    cat("...\n")
    dat <- prepare_stan_data(f = formula,  data = data,  id_var = id_var,
           coord_var = coord_var, cluster_var = cluster_var,
           urban_rural_var = urban_rural_var, grid_x = grid_x,
           grid_y = grid_y, sigma_u_prior = sigma_u_prior,
           sigma_r_prior = sigma_r_prior, sigma_nu_prior = sigma_nu_prior,
           beta_prior = beta_prior, zeta_prior = zeta_prior)
    out <- run_fnc_binary(data = dat, cluster_var = cluster_var,
                          digits = digits, draws = draws)
    out$family <- family
    res <- NULL
    res$results = out
    res$data = dat$data
    res$sp_data = dat$sp_data
    class(res) <- "svyLocAdj"
  }
  if(family%in%"gaussian"){
    cat(paste0("Bayesian model for continuous outcome ...\n"))
    cat("Running Maximum A Posteriori (MAP) estimate ...\n")
    cat("...\n")
    dat <- prepare_stan_data(f = formula,  data = data,  id_var = id_var,
                             coord_var = coord_var, cluster_var = cluster_var,
                             urban_rural_var = urban_rural_var, grid_x = grid_x,
                             grid_y = grid_y, sigma_u_prior = sigma_u_prior,
                             sigma_r_prior = sigma_r_prior, sigma_nu_prior = sigma_nu_prior,
                             beta_prior = beta_prior, zeta_prior = zeta_prior)
    out <- run_fnc_gaussian(data = dat, cluster_var = cluster_var,
                          digits = digits, draws = draws)
    out$family <- family
    res <- NULL
    res$results = out
    res$data = dat$data
    res$sp_data = dat$sp_data
    class(res) <- "svyLocAdj"
  }
  if(family%in%"poisson"){
    cat(paste0("Bayesian model for count data ...\n"))
    cat("Running Maximum A Posteriori (MAP) estimate ...\n")
    cat("...\n")
    dat <- prepare_stan_data(f = formula,  data = data,  id_var = id_var,
                             coord_var = coord_var, cluster_var = cluster_var,
                             urban_rural_var = urban_rural_var, grid_x = grid_x,
                             grid_y = grid_y, sigma_u_prior = sigma_u_prior,
                             sigma_r_prior = sigma_r_prior, sigma_nu_prior = sigma_nu_prior,
                             beta_prior = beta_prior, zeta_prior = zeta_prior)
    out <- run_fnc_poisson(data = dat, cluster_var = cluster_var,
                            digits = digits, draws = draws)
    out$family <- family
    res <- NULL
    res$results = out
    res$data = dat$data
    res$sp_data = dat$sp_data
    class(res) <- "svyLocAdj"
  }
  #if(family!="binomial" | family!="gaussian"){
  #  stop("family can take 'binomial' and 'gaussian' options")
  #}
  cat("Finished.\n")
  res
  ##
}
