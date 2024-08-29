\name{reduce_dataset_based_on_selected_whitening}
\alias{reduce_dataset_based_on_selected_whitening}
\title{Reduce dataset size based on selected whitening method}
\description{
  This function reduces the size of a dataset based on the specified whitening method. It selects a subset of the data based on predefined row limits associated with different whitening methods: ZCA-correlation, PCA-correlation, or Cholesky.
}
\usage{
  reduce_dataset_based_on_selected_whitening(
    whitening_method,
    raw_pheno_df,
    nrow_lim_raw_dataset_zca_cor,
    nrow_lim_raw_dataset_pca_cor,
    nrow_lim_raw_dataset_chol
  )
}
\arguments{
  \item{whitening_method}{A character string specifying the whitening method used to determine the dataset reduction approach. Possible values are \code{"ZCA-cor"}, \code{"PCA-cor"}, or \code{"Cholesky"}.}
  \item{raw_pheno_df}{A data frame containing the raw phenotypic data to be reduced.}
  \item{nrow_lim_raw_dataset_zca_cor}{An integer specifying the row limit for the dataset reduction when using the ZCA-correlation whitening method.}
  \item{nrow_lim_raw_dataset_pca_cor}{An integer specifying the row limit for the dataset reduction when using the PCA-correlation whitening method.}
  \item{nrow_lim_raw_dataset_chol}{An integer specifying the row limit for the dataset reduction when using the Cholesky whitening method.}
}
\value{
  A data frame containing the reduced dataset based on the selected whitening method. The number of rows in the reduced dataset is constrained by the row limit defined for the chosen whitening method.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  
  # load example data
  data("example_raw_pheno_df")
  
  # define parameters
  whitening_method <- "PCA-cor"
  nrow_lim_raw_dataset_zca_cor <- 1000
  nrow_lim_raw_dataset_pca_cor <- 500
  nrow_lim_raw_dataset_chol <- 750
  
  # reduce dataset based on the selected whitening method
  reduced_pheno_df <- reduce_dataset_based_on_selected_whitening(
    whitening_method = whitening_method,
    raw_pheno_df = example_raw_pheno_df,
    nrow_lim_raw_dataset_zca_cor = nrow_lim_raw_dataset_zca_cor,
    nrow_lim_raw_dataset_pca_cor = nrow_lim_raw_dataset_pca_cor,
    nrow_lim_raw_dataset_chol = nrow_lim_raw_dataset_chol
  )
  
  # print the dimensions of the reduced dataset
  print(dim(reduced_pheno_df))
}
