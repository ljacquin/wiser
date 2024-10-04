\name{compute_whitening_matrix_for_sig_mat_}
\alias{compute_whitening_matrix_for_sig_mat_}
\title{Compute whitening matrix based on selected whitening method}
\description{
  This function computes the whitening matrix for a given covariance matrix, based on the selected whitening method. The function also includes options for regularizing the covariance matrix using different regularization or shrinkage methods. The available whitening methods are ZCA-correlation, PCA-correlation, and Cholesky decomposition.
}
\usage{
  compute_whitening_matrix_for_sig_mat_(
    whitening_method,
    regularization_method,
    parallelized_cholesky,
    sig_mat_, 
    alpha_
  )
}
\arguments{
  \item{whitening_method}{A character string specifying the whitening method to use. Options are \code{"ZCA-cor"}, \code{"PCA-cor"}, or \code{"Cholesky"}.}
  \item{regularization_method}{A character string specifying the regularization or shrinkage method used to ensure the positive definiteness of the covariance matrix. Options include \code{"frobenius_norm_regularization"}, \code{"frobenius_norm_shrinkage"}, \code{"frobenius_norm_partial_shrinkage"}, \code{"trace_sample_variance_regularization"}, \code{"trace_sample_variance_shrinkage"}, and \code{"trace_sample_variance_partial_shrinkage"}.}
  \item{parallelized_cholesky}{A logical value indicating whether to use parallelized computation for Cholesky decomposition. This option is applicable only if \code{whitening_method} is \code{"Cholesky"}.}
  \item{sig_mat_}{A numeric matrix representing the covariance matrix that will be used to compute the whitening matrix.}
  \item{alpha_}{A numeric value specifying the regularization parameter associated with the selected regularization or shrinkage method.}
}
\value{
  A list containing:
  \item{sig_mat_}{The regularized covariance matrix, modified based on the selected regularization or shrinkage method.}
  \item{w_mat}{A numeric matrix representing the whitening matrix, computed based on the selected whitening method and the regularized \code{sig_mat_}.}
}
\details{
  This function performs the following steps:
  \itemize{
    \item Regularizes the covariance matrix \code{sig_mat_} based on the specified \code{regularization_method}. This ensures the matrix is positive definite, which is required for whitening.
    \item Computes the whitening matrix \code{w_mat} using the selected \code{whitening_method}:
    \itemize{
      \item For \code{"ZCA-cor"} and \code{"PCA-cor"}, it uses the \code{whiteningMatrix} function to compute the whitening matrix.
      \item For \code{"Cholesky"}, it computes the whitening matrix using Cholesky decomposition. If \code{parallelized_cholesky} is \code{TRUE}, the decomposition is performed in parallel.
    }
  }
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage of compute_whitening_matrix_for_sig_mat_
  
  # define example parameters
  whitening_method <- "ZCA-cor"
  regularization_method <- "frobenius_norm_shrinkage"
  parallelized_cholesky <- FALSE
  sig_mat_ <- matrix(rnorm(100), 10, 10)
  alpha_ <- 0.01
  
  # compute the whitening matrix
  whitening_matrix <- compute_whitening_matrix_for_sig_mat_(
    whitening_method = whitening_method,
    regularization_method = regularization_method,
    parallelized_cholesky = parallelized_cholesky,
    sig_mat_ = sig_mat_,
    alpha_ = alpha_
  )
  
  # print the whitening matrix
  print(whitening_matrix)
}
