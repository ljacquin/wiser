# function to simulate phenotype data
simulate_y <- function(x_mat, z_mat, beta_hat, sigma2_u, sigma2_e, k_mat) {
  # get incidence matrices dimensions
  n <- nrow(x_mat)
  q <- ncol(z_mat)
  
  # simulate u and eps
  u <- mvrnorm(1, mu = rep(0, q), Sigma = sigma2_u * k_mat)
  eps <- rnorm(n, mean = 0, sd = sqrt(sigma2_e))
  
  # compute simulated y
  y_sim <- x_mat %*% beta_hat + z_mat %*% u + eps
  return(y_sim)
}
