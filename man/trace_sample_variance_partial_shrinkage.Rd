\name{trace_sample_variance_partial_shrinkage}
\alias{trace_sample_variance_partial_shrinkage}
\title{Partial shrinkage of diagonal elements of a covariance matrix using the trace}
\description{
  This function applies partial shrinkage to the diagonal elements of a covariance matrix by adjusting them based on the trace (sum of diagonal elements) of the matrix. The off-diagonal elements remain unchanged. The shrinkage intensity is controlled by \code{alpha_}.
  
  The shrinkage is performed according to the formula:
    \deqn{\Sigma_{\text{shrink},ii} = (1 - \alpha) \Sigma_{ii} + \alpha \text{tr}(\Sigma)}
  where \eqn{\text{tr}(\Sigma)} is the trace of the covariance matrix \eqn{\Sigma}. The off-diagonal elements of \eqn{\Sigma} remain unchanged.
}
\usage{
  trace_sample_variance_partial_shrinkage(cov_mat_, alpha_)
}
\arguments{
  \item{cov_mat_}{A square numeric matrix representing the sample covariance matrix.}
  \item{alpha_}{A numeric value between 0 and 1 that controls the shrinkage intensity for the diagonal elements. A value of 0 corresponds to no shrinkage, and a value of 1 fully replaces the diagonal elements with the scaled trace of the matrix.}
}
\value{
  A covariance matrix with its diagonal elements partially shrunk towards the trace of the matrix, while the off-diagonal elements remain unchanged.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # Example usage of trace_sample_variance_partial_shrinkage
  
  sample_cov <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  
  # Apply partial shrinkage with alpha = 0.3
  partially_shrunk_cov <- trace_sample_variance_partial_shrinkage(sample_cov, alpha_ = 0.3)
  
  # Print the partially shrunken covariance matrix
  print(partially_shrunk_cov)
}
