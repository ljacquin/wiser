\name{optimize_whitening_and_regularization}
\alias{optimize_whitening_and_regularization}
\title{Optimization of whitening methods and regularization parameters for WISER phenotype estimation}
\description{
  This function optimizes the selection of whitening methods and regularization parameters for phenotype estimation using the WISER method. It performs a grid search over specified combinations of whitening methods, regularization parameters, and prediction methods, using k-fold cross-validation to evaluate and identify the best-performing combination of whitening methods and regularization parameters.
}
\usage{
  optimize_whitening_and_regularization(
    geno_df, raw_pheno_df, trait_,
    fixed_effects_vars = c("Envir", "Country", "Year", "Row", "Position", "Management"),
    random_effects_vars = "Genotype",
    prediction_method = c("rf", "svr", "gblup", "rkhs", "lasso"),
    whitening_method_grid = c("ZCA-cor", "PCA-cor", "Cholesky"),
    regularization_method_ = "frobenius_norm",
    alpha_frob_grid = c(0.01, 0.1),
    reduce_raw_dataset_size_ = TRUE,
    nrow_lim_raw_dataset_ = 5e3,
    parallelized_cholesky_ = TRUE,
    k_folds_ = 5)
}
\arguments{
  \item{geno_df}{A data frame containing genomic data, with genotypes as rows and SNPs as columns.}
  \item{raw_pheno_df}{A data frame containing raw phenotypic data with observations across different environments and conditions.}
  \item{trait_}{A character string specifying the target trait to be analyzed and optimized.}
  \item{fixed_effects_vars}{A character vector listing the fixed effects variables to be included in the model. Default is \code{c("Envir", "Country", "Year", "Row", "Position", "Management")}.}
  \item{random_effects_vars}{A character string specifying the random effects variable. Default is \code{"Genotype"}.}
  \item{prediction_method}{A character vector specifying the prediction methods to be evaluated. Default is \code{c("rf", "svr", "gblup", "rkhs", "lasso")}.}
  \item{whitening_method_grid}{A character vector specifying the whitening methods to be evaluated. Default is \code{c("ZCA-cor", "PCA-cor", "Cholesky")}.}
  \item{regularization_method_}{A character string specifying the regularization method. Default is \code{"frobenius_norm"}.}
  \item{alpha_frob_grid}{A numeric vector specifying a grid of regularization parameter values for the Frobenius norm. Default is \code{c(0.01, 0.1)}.}
  \item{reduce_raw_dataset_size_}{A logical value indicating whether to reduce the size of the raw phenotype dataset for computational efficiency. Default is \code{TRUE}.}
  \item{nrow_lim_raw_dataset_}{An integer specifying the maximum number of rows to retain in the dataset after reduction. Default is \code{5e3}.}
  \item{parallelized_cholesky_}{A logical value indicating whether to perform Cholesky decomposition in parallel during whitening. Default is \code{TRUE}.}
  \item{k_folds_}{An integer specifying the number of folds to use in k-fold cross-validation. Default is \code{5}.}
}
\value{
  A list containing:
    \item{opt_results}{A data frame summarizing, for each prediction method, the optimal combination of whitening method and regularization parameter based on k-folds cv.}
    \item{opt_alpha_frob}{The regularization parameter associated to the statistical mode of the optimal combination of whitening method and regularization parameter based on k-folds cv.}
    \item{opt_whitening_method}{The whitening method associated to the statistical mode of the optimal combination of whitening method and regularization parameter based on k-folds cv..}
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage:

  # define inputs
  trait <- "Flowering_begin"
  whitening_method_grid <- c("ZCA-cor", "PCA-cor", "Cholesky")
  alpha_frob_grid <- c(0.01, 0.1)
  prediction_method <- c("rf", "svr", "gblup", "rkhs", "lasso")

  # optimize whitening, regularization, and prediction methods
  opt_results <- optimize_whitening_and_regularization(
    geno_df = refpop_geno_data_subset,
    raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
    trait_ = trait,
    whitening_method_grid = whitening_method_grid,
    alpha_frob_grid = alpha_frob_grid,
    prediction_method = prediction_method,
    k_folds_ = 5
  )

  # print the optimal results
  print(opt_results$opt_results)
  print(paste("optimal alpha_frob:", opt_results$opt_alpha_frob))
  print(paste("optimal whitening method:", opt_results$opt_whitening_method))
}
