# function which applies partial shrinkage to the diagonal elements
# of a covariance matrix by adding a positive diagonal matrix based on
# the Frobenius norm. This does not modify the off-diagonal elements.
frobenius_norm_partial_shrinkage <- function(cov_mat_, alpha_) {
  n <- nrow(cov_mat_)
  
  # create a new matrix where only the diagonal is modified
  cov_mat_reg <- cov_mat_
  diag(cov_mat_reg) <- (1 - alpha_) * diag(cov_mat_) + alpha_ * frobenius_norm(cov_mat_)
  
  return(cov_mat_reg)
}