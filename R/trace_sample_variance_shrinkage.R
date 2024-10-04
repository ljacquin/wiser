# function which makes a covariance matrix positive definite by using
# a convex shrinkage estimator which adds a positive diagonal matrix based on
# the sum of sample variances (i.e. trace of sample covariance matrix)
trace_sample_variance_shrinkage <- function(cov_mat_, alpha_) {
  n <- nrow(cov_mat_)
  cov_mat_ <- (1 - alpha_) * cov_mat_ + alpha_ * trace_mat(cov_mat_) * diag(n)
  return(cov_mat_)
}
