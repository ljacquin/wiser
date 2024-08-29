\name{estimate_wiser_phenotype}
\alias{estimate_wiser_phenotype}
\title{Estimate phenotypes approximating genetic values using whitening, successive OLS, and ABC}
\description{
  Estimates phenotypes that approximate genetic values by applying a combination of whitening, successive ordinary least squares (OLS), and approximate Bayesian computation (ABC) for variance component estimation. The function iteratively refines variance component and phenotypic value estimates through multiple rounds of ABC.
}
\usage{
  estimate_wiser_phenotype(
    omic_df, raw_pheno_df, trait_,
    fixed_effects_vars = c(
      "Envir", "Country", "Year",
      "Row", "Position", "Management"
    ),
    fixed_effects_vars_computed_as_factor = c(
      "Envir", "Country", "Year",
      "Row", "Position", "Management"
    ),
    site_var = "Country",
    fixed_effects_vars_computed_as_factor_by_site = c("Row", "Position"),
    random_effects_vars = "Genotype",
    init_sigma2_u = 1,
    init_sigma2_e = 1,
    n_sim_abc = 100,
    seed_abc = 123,
    quantile_threshold_abc = 0.05,
    nb_iter_abc = 1,
    kernel_type = "linear",
    whitening_method = "ZCA-cor",
    regularization_method = "frobenius_norm",
    alpha_frob_ = 0.01,
    percent_eig_ = 0.05,
    non_zero_precision_eig_ = 1e-5,
    parallelized_cholesky = TRUE,
    reduce_raw_dataset_size_ = TRUE,
    nrow_lim_raw_dataset_zca_cor = 10e3,
    nrow_lim_raw_dataset_pca_cor = 10e3,
    nrow_lim_raw_dataset_chol = 40e3
  )
}
\arguments{
  \item{omic_df}{A data frame containing the omic data (e.g. genotypes as rows and SNPs, metabolites, or wavelength reflectances as columns).}
  \item{raw_pheno_df}{A data frame containing the raw phenotypic data.}
  \item{trait_}{A character string specifying the trait of interest.}
  \item{fixed_effects_vars}{A character vector specifying the fixed effects variables. Defaults to \code{c("Envir", "Country", "Year", "Row", "Position", "Management")}.}
  \item{fixed_effects_vars_computed_as_factor}{A character vector specifying which fixed effects variables should be treated as factors. Defaults to \code{c("Envir", "Country", "Year", "Row", "Position", "Management")}. Note: declared fixed effects variables which are not specified to be computed as factors will be treated as quantitative variables.}
  \item{site_var}{A character string specifying the site variable. Defaults to \code{"Country"}.}
  \item{fixed_effects_vars_computed_as_factor_by_site}{A character vector specifying which fixed effects variables should be treated as factors by site. Defaults to \code{c("Row", "Position")}.}
  \item{random_effects_vars}{A character vector specifying the random effects variables. Defaults to \code{"Genotype"}.}
  \item{init_sigma2_u}{A numeric value for the initial variance component of random effects. Defaults to \code{1}.}
  \item{init_sigma2_e}{A numeric value for the initial variance component of residuals. Defaults to \code{1}.}
  \item{n_sim_abc}{An integer specifying the number of simulations for ABC. Defaults to \code{100}.}
  \item{seed_abc}{An integer to set the seed for random number generation. Defaults to \code{123}.}
  \item{quantile_threshold_abc}{A numeric value for the quantile threshold in ABC. Defaults to \code{0.05}.}
  \item{nb_iter_abc}{An integer specifying the number of iterations for ABC. Defaults to \code{1}.}
  \item{kernel_type}{A character string specifying the kernel type. Can be \code{"linear"} or \code{"identity"}. Defaults to \code{"linear"}.}
  \item{whitening_method}{A character string specifying the whitening method. Can be \code{"ZCA-cor"}, \code{"PCA-cor"}, or \code{"Cholesky"}. Defaults to \code{"ZCA-cor"}.}
  \item{regularization_method}{A character string specifying the regularization method for the whitening. Can be \code{"mean_small_eigenvalues"}, \code{"mean_eigenvalues"} or \code{"frobenius_norm"}. Defaults to \code{"frobenius_norm"}.}
  \item{alpha_frob_}{A numeric value for the regularization parameter in the Frobenius norm. Defaults to \code{0.01}.}
  \item{percent_eig_}{A numeric value between 0 and 1 representing the percentage of the smallest strictly positive eigenvalues to be considered. Defaults to \code{0.05}.}
  \item{non_zero_precision_eig_}{A small positive value used to filter out near-zero eigenvalues. Defaults to \code{1e-5}.}
  \item{parallelized_cholesky}{A logical value indicating whether to use parallelized Cholesky decomposition. Defaults to \code{TRUE}.}
  \item{reduce_raw_dataset_size_}{A logical value indicating whether to reduce the size of the raw dataset. Defaults to \code{TRUE}.}
  \item{nrow_lim_raw_dataset_zca_cor}{An integer specifying the maximum number of rows in the raw dataset when using the ZCA-cor whitening method. Defaults to \code{10e3}.}
  \item{nrow_lim_raw_dataset_pca_cor}{An integer specifying the maximum number of rows in the raw dataset when using the PCA-cor whitening method. Defaults to \code{10e3}.}
  \item{nrow_lim_raw_dataset_chol}{An integer specifying the maximum number of rows in the raw dataset when using the Cholesky whitening method. Defaults to \code{40e3}.}
}
\value{
  A list containing:
  \item{wiser_omic_data}{The omic data used in the estimation process.}
  \item{sig_mat_u}{The computed variance-covariance matrix for the individual random effects (i.e. \eqn{\Sigma_u}).}
  \item{wiser_fixed_effect_estimates}{A data frame of the estimated fixed effects.}
  \item{wiser_abc_variance_component_estimates}{The results of the ABC variance component estimation.}
  \item{wiser_phenotypes}{A data frame containing the estimated phenotypic values (i.e. \code{v_hat}), for the genotypes, from the second OLS regression.}
  \item{wiser_z_mat}{The design matrix for the random effects.}
  \item{wiser_x_mat}{The design matrix for the fixed effects.}
  \item{wiser_x_mat_tilde}{The transformed design matrix after whitening.}
  \item{wiser_xi_hat}{The estimated residuals from the first OLS regression.}
  \item{wiser_y_hat}{The predicted values from the first OLS regression.}
  \item{wiser_y}{The raw observed phenotype used by wiser.}
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage

  # -- load wiser library and data subsets from the refpop dataset
  # load library
  library(wiser)

  # load raw individual phenotypic measurements (for subset of 80 genotypes from refpop data)
  data("refpop_raw_indiv_pheno_data_subset")
  head(refpop_raw_indiv_pheno_data_subset)

  # load SNP marker data (for subset of 80 genotypes from refpop data)
  data("refpop_geno_data_subset")
  head(refpop_geno_data_subset)[,1:10]

  # -- estimate wiser phenotype for a specified trait, e.g., Flowering_begin
  # define trait
  trait_ <- "Flowering_begin"

  # apply wiser estimation function
  wiser_obj <- estimate_wiser_phenotype(
    omic_df = refpop_geno_data_subset,
    raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
    trait_ = trait_,
    fixed_effects_vars = c("Envir", "Country", "Year","Row", "Position", "Management"),
    fixed_effects_vars_computed_as_factor = c("Envir", "Country", "Year","Row", "Position", "Management"),
    site_var = "Country",
    fixed_effects_vars_computed_as_factor_by_site = c("Row", "Position"),
    random_effects_vars = "Genotype",
    kernel_type = "linear",
    whitening_method = "ZCA-cor"
  )

  # -- plot the density for the estimated phenotypes
  dev.new()
  plot(density(wiser_obj$wiser_phenotypes$v_hat), main = paste0(trait_, " v_hat"))

  # -- get the estimated fixed effects from the whitening process and OLS, and variance components from ABC
  # estimated fixed effects (from whitening and OLS)
  print(wiser_obj$wiser_fixed_effect_estimates)

  # estimated variance components (from ABC)
  print(wiser_obj$wiser_abc_variance_component_estimates)
}
