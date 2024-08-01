\name{estimate_wiser_phenotype}
\alias{estimate_wiser_phenotype}
\title{Estimate phenotypes approximating genetic values using whitening, successive OLS, and ABC}
\description{
  Estimates phenotypes approximating genetic values using whitening, successive ordinary least squares (OLS), and approximate Bayesian computation (ABC) for variance component estimation. The function performs multiple iterations of ABC to refine the estimates of variance components and phenotypic values.
}
\usage{
  estimate_wiser_phenotype(
    geno_df, raw_pheno_df, trait_,
    fixed_effects_vars = c("Envir", "Country", "Year", "Row", "Position", "Management"),
    compute_row_and_position_as_factors = T,
    random_effects_vars = "Genotype",
    init_sigma2_u = 1, init_sigma2_e = 1,
    n_sim_abc = 100, seed_abc = 123, quantile_threshold_abc = 0.05,
    nb_iter_abc = 1, kernel_type = "linear",
    rate_decay_kernel = 0.1, whitening_method = "Cholesky",
    regularization_method = "frobenius_norm",
    alpha_frob_ = 0.01, percent_eig_ = 0.05, non_zero_precision_eig_ = 1e-5
  )
}
\arguments{
  \item{geno_df}{A data frame of genomic data with genotypes as rows and SNPs as columns.}
  \item{raw_pheno_df}{A data frame of raw phenotypic data.}
  \item{trait_}{A character string specifying the trait of interest.}
  \item{fixed_effects_vars}{A character vector of fixed effects variables. Defaults to \code{c("Envir", "Country", "Year", "Row", "Position", "Management")}.}
  \item{compute_row_and_position_as_factors}{A logical value indicating whether to treat rows and positions as factors. Defaults to \code{FALSE}.}
  \item{random_effects_vars}{A character vector of random effects variables. Defaults to \code{"Genotype"}.}
  \item{init_sigma2_u}{A numeric value for the initial variance component of random effects. Defaults to \code{1}.}
  \item{init_sigma2_e}{A numeric value for the initial variance component of residuals. Defaults to \code{1}.}
  \item{n_sim_abc}{An integer specifying the number of simulations for ABC. Defaults to \code{100}.}
  \item{seed_abc}{An integer to set the seed for random number generation. Defaults to \code{123}.}
  \item{quantile_threshold_abc}{A numeric value for the quantile threshold in ABC. Defaults to \code{0.05}.}
  \item{nb_iter_abc}{An integer specifying the number of iterations for ABC. Defaults to \code{1}.}
  \item{kernel_type}{A character string specifying the kernel type. Can be \code{"linear"}, \code{"gaussian"}, or \code{"identity"}. Defaults to \code{"linear"}.}
  \item{rate_decay_kernel}{A numeric value for the rate decay parameter in the Gaussian kernel. Defaults to \code{0.1}.}
  \item{whitening_method}{A character string specifying the whitening method. Can be \code{"Cholesky"} or \code{"ZCA-cor"}. Defaults to \code{"Cholesky"}.}
  \item{regularization_method}{A character string specifying the regularization method for the whitening. Can be \code{"frobenius_norm"} or other methods. Defaults to \code{"frobenius_norm"}.}
  \item{alpha_frob_}{A numeric value for the regularization parameter in Frobenius norm. Defaults to \code{0.01}.}
  \item{percent_eig_}{A numeric value for the percentage of eigenvalues retained in regularization. Defaults to \code{0.05}.}
  \item{non_zero_precision_eig_}{A numeric value for the precision of non-zero eigenvalues in regularization. Defaults to \code{1e-5}.}
}
\value{
  A list containing:
  \item{var_comp_abc_obj}{The results of the ABC variance component estimation.}
  \item{beta_hat}{The estimated fixed effects, adjusted for the genetic covariance structure between individuals in the experimental design.}
  \item{v_hat}{The estimated phenotypic values.}
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # Example usage
  # -- load wiser library and data subsets from the refpop dataset
  # load library
  library(wiser)

  # load raw individual phenotypic measurements (for subset of 80 genotypes from refpop data)
  data("refpop_raw_indiv_pheno_data_subset")
  head(refpop_raw_indiv_pheno_data_subset)

  # load SNP marker data (for subset of 80 genotypes from refpop data)
  data("refpop_geno_data_subset")
  head(refpop_geno_data_subset)[,1:10]

  # -- estimate wiser phenotype for a specified trait,i.e. Flowering_begin here
  # define trait
  trait_ <- "Flowering_begin"

  # apply wiser estimation function
  pheno_obj <- estimate_wiser_phenotype(
    geno_df = refpop_geno_data_subset,
    raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
    trait_ = trait_,
    fixed_effects_vars = c("Envir", "Country", "Year", "Row", "Position", "Management"),
    compute_row_and_position_as_factors = T,
    random_effects_vars = "Genotype",
    init_sigma2_u = 1, init_sigma2_e = 1,
    n_sim_abc = 100, seed_abc = 123, quantile_threshold_abc = 0.05,
    nb_iter_abc = 1, kernel_type = "linear",
    rate_decay_kernel = 0.1, whitening_method = "Cholesky",
    regularization_method = "frobenius_norm",
    alpha_frob_ = 0.01, percent_eig_ = 0.05, non_zero_precision_eig_ = 1e-5
  )

  # -- plot the density for the estimated phenotypes
  dev.new()
  plot(density(pheno_obj$v_hat), main = paste0(trait_, " v_hat"))

  # -- get the estimated fixed effects, from the whitening process and OLS and variance components from ABC

  # estimated fixed effects (from whitening and OLS)
  print(pheno_obj$beta_hat)

  # estimated variance components (from ABC)
  print(pheno_obj$var_comp_abc_obj)
}
