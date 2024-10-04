# abc function to compute variance components
abc_variance_component_estimation <- function(y, x_mat, z_mat, k_mat, beta_hat,
                                              prior_sigma2_u, prior_sigma2_e,
                                              n_sim_abc, seed_abc,
                                              quantile_threshold_abc) {
  df_results <- future.apply::future_lapply(
    1:n_sim_abc,
    future.seed = T,
    function(sim_num) {
      set.seed(sim_num * seed_abc)
      simulate_and_compute_squared_l2_norm(
        y, x_mat, z_mat, k_mat, beta_hat, prior_sigma2_u, prior_sigma2_e
      )
    },
    future.packages = c("MASS")
  )
  df_results <- do.call(rbind, df_results)

  # assign colnames to df_results
  df_results <- as.data.frame(df_results)
  colnames(df_results) <- c("sigma2_u_hat", "sigma2_e_hat", "distance")

  # extract df_results
  vect_distances <- as.numeric(df_results[, "distance"])

  # get rejection threshold based on define quantile_threshold_abc
  reject_thresh <- quantile(vect_distances, quantile_threshold_abc)

  # get accepted variance components parameters for rejection threshold
  accepted_params <- as.data.frame(
    df_results[vect_distances <= reject_thresh, ]
  )

  # compute the average of the accepted parameters
  sigma2_u_hat_mean <- mean(accepted_params[, "sigma2_u_hat"])
  sigma2_e_hat_mean <- mean(accepted_params[, "sigma2_e_hat"])

  return(list(
    "complete_results" = df_results,
    "sigma2_u_hat_mean" = sigma2_u_hat_mean,
    "sigma2_e_hat_mean" = sigma2_e_hat_mean,
    "accepted_params" = accepted_params,
    "rejection_threshold" = reject_thresh
  ))
}
