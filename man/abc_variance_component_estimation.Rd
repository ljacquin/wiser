\name{abc_variance_component_estimation}
\alias{abc_variance_component_estimation}
\title{Estimate variance components using approximate bayesian computation (ABC)}
\description{
  Estimates variance components by simulating phenotypes and distances using approximate bayesian computation (ABC) and computing the mean of accepted parameters based on a rejection threshold.
}
\usage{
  abc_variance_component_estimation(y, x_mat, z_mat, k_mat, beta_hat,
                                    prior_sigma2_u, prior_sigma2_e,
                                    n_sim_abc, seed_abc,
                                    quantile_threshold_abc)
}
\arguments{
  \item{y}{A numeric vector of observed phenotypic measurements.}
  \item{x_mat}{A numeric matrix of incidence variables for fixed effects.}
  \item{z_mat}{A numeric matrix of incidence variables for random effects.}
  \item{k_mat}{A numeric matrix representing the genomic covariance matrix.}
  \item{beta_hat}{A numeric vector of estimated fixed effects.}
  \item{prior_sigma2_u}{A numeric vector specifying prior values for \eqn{\sigma^2_u}.}
  \item{prior_sigma2_e}{A numeric vector specifying prior values for \eqn{\sigma^2_e}.}
  \item{n_sim_abc}{An integer specifying the number of simulations for ABC.}
  \item{seed_abc}{An integer to set the seed for random number generation.}
  \item{quantile_threshold_abc}{A numeric value for the quantile threshold to determine the rejection region in ABC.}
}
\value{
  A list containing:
    \item{complete_results}{A data frame with columns for variance component estimates and distances.}
  \item{sigma2_u_hat_mean}{The mean of accepted \eqn{\sigma^2_u} estimates.}
  \item{sigma2_e_hat_mean}{The mean of accepted \eqn{\sigma^2_e} estimates.}
  \item{accepted_params}{A data frame of accepted variance component parameters.}
  \item{rejection_threshold}{The rejection threshold based on the quantile threshold.}
}
\examples{
  # Example usage
  # result <- abc_variance_component_estimation(y, x_mat, z_mat, k_mat, beta_hat,
  #                                              prior_sigma2_u, prior_sigma2_e,
  #                                              n_sim_abc = 100, seed_abc = 123,
  #                                              quantile_threshold_abc = 0.05)
}
