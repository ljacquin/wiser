\name{regularize_covariance_mean_small_eigenvalues}
\alias{regularize_covariance_mean_small_eigenvalues}
\title{Regularize a covariance matrix using the mean of a percentage of its smallest eigenvalues to ensure positive definiteness }
\description{
  Regularizes a covariance matrix to ensure it is positive definite by adding a positive value \code{delta} to the diagonal. The value of \code{delta} is computed as the mean of a percentage of its smallest strictly positive eigenvalues, scaled by \code{sigma2_u}.
}
\usage{
  regularize_covariance_mean_small_eigenvalues(
    cov_mat_, k_mat,
    sigma2_u, percent_eig_, non_zero_precision_eig_
  )
}
\arguments{
  \item{cov_mat_}{A covariance matrix between individuals to be regularized.}
  \item{k_mat}{A Gram matrix, or genomic covariance matrix, used to compute the eigenvalues for regularization.}
  \item{sigma2_u}{A scalar value representing the variance component associated with \code{k_mat}.}
  \item{percent_eig_}{A numeric value between 0 and 1 representing the percentage of the smallest strictly positive eigenvalues to be considered, e.g. 0.05 represents the first 5\% of the smallest strictly positive eigenvalues}
  \item{non_zero_precision_eig_}{A small positive value used to filter out near-zero eigenvalues.}
}

\value{
  A numeric matrix representing the regularized covariance matrix, which is positive definite.
}

\details{
  To ensure that the covariance matrix \code{cov_mat_} is positive definite, a regularization term \code{delta} is added to the diagonal elements of the matrix. The value of \code{delta} is computed as follows:

  1. Compute the eigenvalues \eqn{\lambda_i(k\_mat)} of the matrix \code{k_mat} scaled by \code{sigma2_u}:

  \deqn{\lambda_i = \sigma^2_u \times \lambda_i(k\_mat) }

  where \eqn{\lambda_i} are the eigenvalues of \code{cov_mat_}

  2. Compute \eqn{\delta} as the mean of the \eqn{N_{s}} smallest strictly positive eigenvalues, associated to the first \code{percent_eig_}\% of them:

  \deqn{\text{Mean of the } N_s \text{ smallest strictly positive eigenvalues} = \frac{1}{N} \sum_{i=1}^{N} \lambda_i}

  3. Compute the regularized covariance matrix \eqn{\Sigma_{\text{reg}}} by adding this mean value to the diagonal of \eqn{\Sigma}:

  \deqn{\Sigma_{\text{reg}} = \Sigma + \delta \mathbf{I}}

  where \eqn{\delta} is the mean of the smallest eigenvalues and \eqn{\mathbf{I}} is the identity matrix of the same dimension as \eqn{\Sigma}. This ensures that the resulting matrix \eqn{\Sigma_{\text{reg}}} is positive definite.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # Example usage
  # simulate a covariance matrix
  set.seed(123)
  n <- 5
  cov_mat_ <- matrix(rnorm(n * n), n, n)
  cov_mat_ <- cov_mat_ %*% t(cov_mat_)  # make it symmetric positive definite

  # simulate a gram matrix
  k_mat <- matrix(rnorm(n * n), n, n)
  k_mat <- k_mat %*% t(k_mat_)  # make it symmetric positive definite

  # set parameters for the regularization
  sigma2_u <- 1
  percent_eig_ <- 0.5
  non_zero_precision_eig_ <- 1e-6

  # apply the regularization function
  regularized_cov_mat <- regularize_covariance_mean_small_eigenvalues(
    cov_mat_, k_mat, sigma2_u, percent_eig_, non_zero_precision_eig_
  )
  print(regularized_cov_mat)
}
