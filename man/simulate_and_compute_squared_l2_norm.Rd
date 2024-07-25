\name{simulate_and_compute_squared_l2_norm}
\alias{simulate_and_compute_squared_l2_norm}
\title{Simulate phenotypes and compute distance between simulated and observed values}
\description{
  Simulates phenotypic values based on random draws for variance components and computes the squared L2 norm distance between the simulated and observed phenotypic values.
}
\usage{
  simulate_and_compute_squared_l2_norm(y, x_mat, z_mat, k_mat, beta_hat,
                                       prior_sigma2_u, prior_sigma2_e)
}
\arguments{
  \item{y}{A numeric vector of observed phenotypic measurements.}
  \item{x_mat}{A numeric matrix of incidence variables for fixed effects.}
  \item{z_mat}{A numeric matrix of incidence variables for random effects.}
  \item{k_mat}{A numeric matrix representing the genomic covariance matrix.}
  \item{beta_hat}{A numeric vector of estimated fixed effects.}
  \item{prior_sigma2_u}{A numeric vector specifying prior values for \eqn{\sigma^2_u}.}
  \item{prior_sigma2_e}{A numeric vector specifying prior values for \eqn{\sigma^2_e}.}
}
\value{
  A numeric vector containing:
    \item{sigma2_u}{The sampled variance component for random effects.}
  \item{sigma2_e}{The sampled variance component for residuals.}
  \item{dist_y_y_sim}{The squared L2 norm distance between observed and simulated phenotypic values.}
}
\examples{
  # Example usage
  # result <- simulate_and_compute_squared_l2_norm(y, x_mat, z_mat, k_mat, beta_hat,
  #                                                prior_sigma2_u = c(1e-2, 1),
  #                                                prior_sigma2_e = c(1e-2, 1))
}
