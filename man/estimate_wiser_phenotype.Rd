\name{estimate_wiser_phenotype}
\alias{estimate_wiser_phenotype}
\title{Estimate phenotypes approximating genetic values using whitening, successive OLS and ABC}
\description{
  Estimates phenotypes approximating genetic values using whitening, successive ordinary least squares (OLS) and approximate bayesian computation (ABC) for variance component estimation.
}
\usage{
  estimate_wiser_phenotype(
    geno_df, raw_pheno_df, trait_,
    fixed_effects_vars = c("Envir", "Country", "Year", "Row", "Position", "Management"),
    random_effects_vars = "Genotype",
    init_sigma2_u = 1, init_sigma2_e = 1,
    n_sim_abc = 100, seed_abc = 123, quantile_threshold_abc = 0.05,
    nb_iter_abc = 1, kernel_type = "linear",
    rate_decay_kernel = 0.1, whitening_method = "Cholesky"
  )
}
\arguments{
  \item{geno_df}{A data frame of genomic data with genotypes as rows and SNPs as columns.}
  \item{raw_pheno_df}{A data frame of raw phenotypic data.}
  \item{trait_}{A character string specifying the trait of interest.}
  \item{fixed_effects_vars}{A character vector of fixed effects variables.}
  \item{random_effects_vars}{A character vector of random effects variables.}
  \item{init_sigma2_u}{A numeric value for the initial variance component of random effects.}
  \item{init_sigma2_e}{A numeric value for the initial variance component of residuals.}
  \item{n_sim_abc}{An integer specifying the number of simulations for ABC.}
  \item{seed_abc}{An integer to set the seed for random number generation.}
  \item{quantile_threshold_abc}{A numeric value for the quantile threshold in ABC.}
  \item{nb_iter_abc}{An integer specifying the number of iterations for ABC.}
  \item{kernel_type}{A character string specifying the kernel type ("linear", "gaussian", or "identity").}
  \item{rate_decay_kernel}{A numeric value for the rate decay parameter in the Gaussian kernel.}
  \item{whitening_method}{A character string specifying the whitening method ("Cholesky" or "ZCA-cor").}
}
\value{
  A list containing:
    \item{var_comp_abc_obj}{The results of the ABC variance component estimation.}
  \item{beta_hat}{The estimated fixed effects (adjusted for the genetic covariance structure between individuals in the experimental design).}
  \item{v_hat}{The estimated phenotypic values.}
}
\examples{
  # Example usage
  # result <- estimate_wiser_phenotype(
  #   geno_df, raw_pheno_df, trait_ = "Trait1",
  #   fixed_effects_vars = c("Envir", "Country"),
  #   random_effects_vars = "Genotype", init_sigma2_u = 1,
  #   init_sigma2_e = 1, n_sim_abc = 100, seed_abc = 123,
  #   quantile_threshold_abc = 0.05, nb_iter_abc = 1,
  #   kernel_type = "linear", rate_decay_kernel = 0.1,
  #   whitening_method = "Cholesky"
  # )
}
