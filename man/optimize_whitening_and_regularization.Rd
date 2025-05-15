\name{optimize_whitening_and_regularization}
\alias{optimize_whitening_and_regularization}
\title{Optimization of whitening and regularization for WISER phenotype estimation}
\description{
  This function optimizes the selection of whitening methods and regularization parameters for phenotype estimation using the WISER method. It performs a grid search over specified combinations of whitening methods, regularization parameters, and prediction methods, using k-fold cross-validation to evaluate and identify the best-performing combination of whitening methods and regularization parameters. It allows parallel execution for faster computation.
}
\usage{
  optimize_whitening_and_regularization(
    omic_df, raw_pheno_df, trait_,
    fixed_effect_vars = c("Envir", "Row", "Position"),
    fixed_effect_vars_computed_as_factor = c("Envir", "Row", "Position"),
    envir_var = "Envir",
    fixed_effect_vars_computed_as_factor_by_envir = c("Row", "Position"),
    random_effect_vars = "Genotype",
    prediction_method = c("rf", "svr", "gblup", "rkhs", "lasso"),
    whitening_method_grid = c("ZCA-cor", "PCA-cor", "Cholesky"),
    regularization_method_ = "frobenius_norm_regularization",
    alpha_grid = c(0.01, 0.1),
    reduce_raw_dataset_size_ = TRUE,
    nrow_approx_lim_raw_dataset_ = 5e3,
    parallelized_cholesky = TRUE,
    k_folds_ = 5,
    nb_cores_ = 12
  )
}
\arguments{
  \item{omic_df}{A data frame containing the omic data (e.g., genotypes as rows and SNPs, metabolites, or wavelength reflectances/absorbances as columns).}
  \item{raw_pheno_df}{A data frame containing raw phenotypic data with observations across different environments and conditions.}
  \item{trait_}{A character string specifying the target trait to be analyzed and optimized.}
  \item{fixed_effect_vars}{A character vector listing the fixed effects variables to be included in the model. Default is \code{c("Envir", "Row", "Position")}.}
  \item{fixed_effect_vars_computed_as_factor}{A character vector listing the fixed effects variables to be treated as factors during modeling. Default is \code{c("Envir", "Row", "Position")}.}
  \item{envir_var}{A character string specifying the variable representing the environment (e.g., location) to be included in the model. Default is \code{"Envir"}.}
  \item{fixed_effect_vars_computed_as_factor_by_envir}{A character vector listing the fixed effects variables to be treated as factors within each environment. Default is \code{c("Row", "Position")}.}
  \item{random_effect_vars}{A character string specifying the random effects variable. Default is \code{"Genotype"}.}
  \item{prediction_method}{A character vector specifying the prediction methods to be evaluated. Default is \code{c("rf", "svr", "gblup", "rkhs", "lasso")}.}
  \item{whitening_method_grid}{A character vector specifying the whitening methods to be evaluated. Default is \code{c("ZCA-cor", "PCA-cor", "Cholesky")}.}
  \item{regularization_method_}{A character string specifying the regularization or shrinkage method. Default is \code{"frobenius_norm_regularization"}.}
  \item{alpha_grid}{A numeric vector specifying a grid of regularization parameter values for the specified regularization or shrinkage method. Default is \code{c(0.01, 0.1)}.}
  \item{reduce_raw_dataset_size_}{A logical value indicating whether to reduce the size of the raw phenotype dataset for computational efficiency. Default is \code{TRUE}.}
  \item{nrow_approx_lim_raw_dataset_}{An integer specifying the approximate maximum number of rows to retain in the dataset after reduction. Default is \code{5e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} used is computed as \code{K = nrow(raw_pheno_df) / nrow_approx_lim_raw_dataset_}.}
  \item{parallelized_cholesky}{A logical value indicating whether to perform Cholesky decomposition in parallel during whitening. Default is \code{TRUE}.}
  \item{k_folds_}{An integer specifying the number of folds to use in k-fold cross-validation. Default is \code{5}.}
  \item{nb_cores_}{An integer specifying the number of cores to use for parallel computation. Default is \code{12}.}
}
\value{
  A list containing:
  \item{opt_results}{A data frame summarizing, for each prediction method, the optimal combination of whitening method and regularization parameter based on k-fold cross-validation.}
  \item{opt_alpha_}{The regularization parameter associated with the most frequently selected combination of whitening method and regularization parameter across prediction methods.}
  \item{opt_whitening_method}{The whitening method associated with the most frequently selected combination of whitening method and regularization parameter across prediction methods.}
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # define inputs
  trait <- "Flowering_begin"
  whitening_method_grid <- c("ZCA-cor", "PCA-cor", "Cholesky")
  alpha_grid <- c(0.01, 0.1)
  prediction_method <- c("rf", "svr", "gblup", "rkhs", "lasso")

  # optimize whitening and regularization parameter
  opt_results <- optimize_whitening_and_regularization(
    omic_df = refpop_geno_data_subset,
    raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
    trait_ = trait,
    whitening_method_grid = whitening_method_grid,
    alpha_grid = alpha_grid,
    prediction_method = prediction_method,
    k_folds_ = 5,
    nb_cores_ = 12
  )

  # print the optimal results
  print(opt_results$opt_results)
  print(paste("Optimal alpha:", opt_results$opt_alpha_))
  print(paste("Optimal whitening method:", opt_results$opt_whitening_method))
}
