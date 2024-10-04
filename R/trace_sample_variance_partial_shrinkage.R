# function which applies partial shrinkage to the diagonal elements
# of a covariance matrix by adding a positive diagonal matrix based on
# the sum of sample variances (i.e. trace of sample covariance matrix).
# This does not modify the off-diagonal elements.
trace_sample_variance_partial_shrinkage <- function(cov_mat_, alpha_) {
  n <- nrow(cov_mat_)

  # create a new matrix where only the diagonal is modified
  cov_mat_reg <- cov_mat_
  diag(cov_mat_reg) <- (1 - alpha_) * diag(cov_mat_) + alpha_ * trace_mat(cov_mat_)

  return(cov_mat_reg)
}
