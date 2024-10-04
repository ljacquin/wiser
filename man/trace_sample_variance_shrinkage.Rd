\name{trace_sample_variance_shrinkage}
\alias{trace_sample_variance_shrinkage}
\title{Shrinkage of covariance matrix towards a scaled identity matrix using trace}
\description{
  This function shrinks (regularizes) the covariance matrix towards a scaled identity matrix, where the scaling factor is determined by the trace of the covariance matrix. The shrinkage intensity is governed by \code{alpha_}.
  
  The formula for this shrinkage is:
  \deqn{\Sigma_{\text{shrink}} = (1 - \alpha) \Sigma + \alpha\text{tr}(\Sigma)I_n}
  where \eqn{\Sigma} is the original covariance matrix, \eqn{\text{tr}(\Sigma)} is the trace (sum of diagonal elements) of the matrix, and \eqn{I_n} is the identity matrix of size \eqn{n}.
}
\usage{
  trace_sample_variance_shrinkage(cov_mat_, alpha_)
}
\arguments{
  \item{cov_mat_}{A square numeric matrix representing the sample covariance matrix.}
  \item{alpha_}{A numeric value between 0 and 1 that controls the shrinkage intensity. A value of 0 corresponds to no shrinkage, while a value of 1 shrinks fully towards a scaled identity matrix.}
}
\value{
  A shrunken covariance matrix, with off-diagonal elements reduced and diagonal elements scaled by the trace of the original matrix.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  sample_cov <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  shrunk_cov <- trace_sample_variance_shrinkage(sample_cov, alpha_ = 0.5)
  print(shrunk_cov)
}
