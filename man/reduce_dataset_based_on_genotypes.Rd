\name{reduce_dataset_based_on_genotypes}
\alias{reduce_dataset_based_on_genotypes}
\title{Reduce dataset based on genotype counts}
\description{
  This function reduces the size of a dataset by sampling rows within each genotype group. It ensures that each genotype is represented and that no sampling is applied if its count is under a limit, while also attempting to limit the overall dataset size based on a specified row limit.
}
\usage{
  reduce_dataset_based_on_genotypes(df_, nrow_approx_lim = 10e3, min_samples = 3)
}
\arguments{
  \item{df_}{A data frame containing the dataset to be reduced. The data frame must include a column named \code{Genotype} that identifies the genotype associated with each row.}
  \item{nrow_approx_lim}{An integer specifying the approximate maximum number of rows to retain in the dataset after reduction. Defaults to \code{10e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} is computed as \code{K = nrow(df_) / nrow_approx_lim}.}
  \item{min_samples}{An integer specifying the limit under which no sampling is applied for a genotype. Defaults to \code{30}.}
}
\value{
  A reduced data frame that retains a subset of rows for each genotype group, ensuring that the total number of rows is approximately limited by \code{nrow_approx_lim} and that each genotype is represented.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  # create a sample data frame
  df <- data.frame(
    Genotype = rep(c("G1", "G2", "G3"), each = 100),
    Trait1 = rnorm(300),
    Trait2 = rnorm(300)
  )

  # reduce the dataset size based on genotype counts
  df_reduced <- reduce_dataset_based_on_genotypes(df, nrow_approx_lim = 100, min_samples = 3)

  # print the reduced dataset
  print(df_reduced)
}
