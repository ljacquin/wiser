# function to simulate phenotypes and compute distance between simulated and
# observed values
simulate_and_compute_squared_l2_norm <- function(y, x_mat, z_mat, k_mat, beta_hat,
                                                 prior_sigma2_u, prior_sigma2_e) {
  # sample random values for variance components for prior ranges
  sigma2_u <- runif(1, prior_sigma2_u[1], prior_sigma2_u[2])
  sigma2_e <- runif(1, prior_sigma2_e[1], prior_sigma2_e[2])
  
  # simulate phenotypes
  y_sim <- simulate_y(x_mat, z_mat, beta_hat, sigma2_u, sigma2_e, k_mat)
  
  # compute distances
  dist_y_y_sim <- squared_l2_norm(y, y_sim)
  
  return(c(sigma2_u, sigma2_e, dist_y_y_sim))
}
