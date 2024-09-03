\name{compute_whitening_matrix_for_sig_mat_}
\alias{compute_whitening_matrix_for_sig_mat_}
\title{Compute whitening matrix based on selected whitening method}
\description{
  This function computes the whitening matrix, for a given covariance matrix, based on the selected whitening method. It includes options for regularizing the covariance matrix and applying different whitening methods, such as ZCA-correlation, PCA-correlation, or Cholesky decomposition.
}
\usage{
  compute_whitening_matrix_for_sig_mat_(
    whitening_method,
    regularization_method,
    parallelized_cholesky,
    sig_mat_,
    k_mat,
    sigma2_u,
    percent_eig_,
    non_zero_precision_eig_,
    alpha_frob_
  )
}
\arguments{
  \item{whitening_method}{A character string specifying the whitening method to use. Options are \code{"ZCA-cor"}, \code{"PCA-cor"}, or \code{"Cholesky"}.}
  \item{regularization_method}{A character string specifying the method for regularizing the covariance matrix. Options are \code{"mean_small_eigenvalues"}, \code{"mean_eigenvalues"}, or \code{"frobenius_norm"}.}
  \item{parallelized_cholesky}{A logical value indicating whether to use parallelized computation for Cholesky decomposition. Applies only if \code{whitening_method} is \code{"Cholesky"}.}
  \item{sig_mat_}{A numeric matrix representing the covariance matrix that will be whitened.}
  \item{k_mat}{A numeric matrix used in regularization when \code{regularization_method} is \code{"mean_small_eigenvalues"}.}
  \item{sigma2_u}{A numeric value representing the variance of the random effects. Used in regularization when \code{regularization_method} is \code{"mean_small_eigenvalues"}.}
  \item{percent_eig_}{A numeric value specifying the percentage of eigenvalues to consider in regularization. Used in \code{mean_small_eigenvalues} regularization.}
  \item{non_zero_precision_eig_}{A numeric value specifying the number of non-zero precision eigenvalues. Used in \code{mean_small_eigenvalues} regularization.}
  \item{alpha_frob_}{A numeric value specifying the Frobenius norm regularization parameter. Used in \code{frobenius_norm} regularization.}
}
\value{
  \item{sig_mat_}{ The numeric covariance matrix used to compute the whitening matrix \code{w_mat}.}
  \item{w_mat}{ A numeric matrix representing the whitening matrix, based on the selected whitening method, associated to \code{sig_mat_}.}
}
\details{
  The function performs the following steps:
    \itemize{
      \item Regularizes the covariance matrix \code{sig_mat_} based on the specified \code{regularization_method}.
      \item Computes the whitening matrix based on the selected \code{whitening_method}.
      \item For \code{"ZCA-cor"} and \code{"PCA-cor"}, it uses the \code{whiteningMatrix} function to compute the whitening matrix.
      \item For \code{"Cholesky"}, it computes the whitening matrix using Cholesky decomposition, optionally in parallel.
    }
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  
  # define example parameters
  whitening_method <- "ZCA-cor"
  regularization_method <- "frobenius_norm"
  parallelized_cholesky <- FALSE
  sig_mat_ <- matrix(rnorm(100), 10, 10)
  k_mat <- matrix(rnorm(100), 10, 10)
  sigma2_u <- 1.0
  percent_eig_ <- 0.95
  non_zero_precision_eig_ <- 5
  alpha_frob_ <- 0.01
  
  # compute the whitening matrix
  whitening_matrix <- compute_whitening_matrix_for_sig_mat_(
    whitening_method = whitening_method,
    regularization_method = regularization_method,
    parallelized_cholesky = parallelized_cholesky,
    sig_mat_ = sig_mat_,
    k_mat = k_mat,
    sigma2_u = sigma2_u,
    percent_eig_ = percent_eig_,
    non_zero_precision_eig_ = non_zero_precision_eig_,
    alpha_frob_ = alpha_frob_
  )
  
  # print the whitening matrix
  print(whitening_matrix)
}
