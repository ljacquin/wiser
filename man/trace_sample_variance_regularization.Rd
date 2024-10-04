\name{trace_sample_variance_regularization}
\alias{trace_sample_variance_regularization}
\title{Regularization of covariance matrix using trace and identity matrix}
\description{
  This function regularizes the covariance matrix by adding a scaled identity matrix, where the scaling factor is determined by the trace of the covariance matrix. The regularization intensity is controlled by \code{alpha_}.
  
  The formula for this regularization is:
    \deqn{\Sigma_{\text{reg}} = \Sigma + \alpha\text{tr}(\Sigma)I_n}
  where \eqn{\Sigma} is the original covariance matrix, \eqn{\text{tr}(\Sigma)} is the trace (sum of diagonal elements) of the matrix, and \eqn{I_n} is the identity matrix of size \eqn{n}.
}
\usage{
  trace_sample_variance_regularization(cov_mat_, alpha_)
}
\arguments{
  \item{cov_mat_}{A square numeric matrix representing the sample covariance matrix.}
  \item{alpha_}{A numeric value that controls the regularization intensity. Larger values result in more regularization.}
}
\value{
  A regularized covariance matrix, where the diagonal elements are scaled by the trace of the original matrix.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  sample_cov <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  regularized_cov <- trace_sample_variance_regularization(sample_cov, alpha_ = 0.1)
  print(regularized_cov)
}
