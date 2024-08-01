\name{regularize_covariance_mean_eigenvalues}
\alias{regularize_covariance_mean_eigenvalues}
\title{Regularize a covariance matrix using the mean of its eigenvalues to ensure positive definiteness }
\description{
  Regularizes a covariance matrix to ensure it is positive definite by adding a positive value \code{delta} to the diagonal. The value of \code{delta} is computed as the mean of the eigenvalues of the covariance matrix.
}
\usage{
  regularize_covariance_mean_eigenvalues(cov_mat_)
}
\arguments{
  \item{cov_mat_}{A numeric matrix representing the covariance matrix to be regularized.}
}
\value{
  A numeric matrix representing the regularized covariance matrix, which is positive definite.
}
\details{
  To ensure that the covariance matrix \code{cov_mat_} is positive definite, a regularization term \code{delta} is added to the diagonal elements of the matrix. The value of \code{delta} is computed as the mean of the eigenvalues of \code{cov_mat_}:
    
    \deqn{\delta = \frac{\text{tr}(\Sigma)}{n} = \frac{\text{tr}(PDP')}{n} = \frac{\text{tr}(P'PD)}{n} = \frac{\text{tr}(D)}{n}}
  
  where  
  
  - \eqn{D} and \eqn{P} are the diagonal eigenvalue matrix and orthogonal change-of-basis matrix respectively,
  
  - \eqn{\text{tr}(\Sigma)} is the trace of the covariance matrix \eqn{\Sigma} which is equal to the sum of eigenvalues (i.e. \eqn{\text{tr}(\Sigma) = \text{tr}(D)}),
  
  - \eqn{n} is the dimension of the matrix \eqn{\Sigma}.
  
The regularized matrix \eqn{\Sigma_{\text{reg}}} is then computed by adding \eqn{\delta} to the diagonal elements of \eqn{\Sigma}:
    
    \deqn{\Sigma_{\text{reg}} = \Sigma + \delta \mathbf{I}}
  
Here, \eqn{\mathbf{I}} denotes the identity matrix of the same dimensions as \eqn{\Sigma}. This regularization step ensures that the resulting matrix \eqn{\Sigma_{\text{reg}}} is positive definite.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # Example usage
  cov_mat_ <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  regularized_cov_mat <- regularize_covariance_mean_eigenvalues(cov_mat_)
  print(regularized_cov_mat)
}
