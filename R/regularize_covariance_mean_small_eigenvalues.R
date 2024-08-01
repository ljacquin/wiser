# function which makes a covariance matrix positive definite by adding a
# positive value delta to the diagonal, based on the mean of the percent_eig_%
# strictly positive smallest eigenvalues
regularize_covariance_mean_small_eigenvalues <- function(
    cov_mat_, k_mat,
    sigma2_u, percent_eig_, non_zero_precision_eig_) {
  # compute the eigen values from k_mat
  eig_val_ <- sigma2_u * (mixOmics::pca(k_mat, ncomp = ncol(k_mat))$sdev^2)
  
  # get the percent_eig_% strictly positive smallest ones
  thresh_ <- ceiling(percent_eig_ * length(eig_val_))
  small_eig_val_ <- sort(eig_val_, decreasing = F)[1:thresh_]
  mean_small_eigen_ <- mean(
    small_eig_val_[small_eig_val_ > non_zero_precision_eig_]
  )
  # compute the regularized covariance matrix
  n <- nrow(cov_mat_)
  cov_mat_ <- cov_mat_ + mean_small_eigen_ * diag(n)
  return(cov_mat_)
}