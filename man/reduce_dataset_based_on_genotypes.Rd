\name{reduce_dataset_based_on_genotypes}
\alias{reduce_dataset_based_on_genotypes}
\title{Reduce dataset based on genotype counts}
\description{
  This function reduces the size of a dataset by sampling rows within each genotype group. It ensures that each genotype is represented by a minimum number of samples, while also attempting to limit the overall dataset size based on a specified row limit.
}
\usage{
  reduce_dataset_based_on_genotypes(df_, nrow_lim = 10e3, min_samples = 3)
}
\arguments{
  \item{df_}{A data frame containing the dataset to be reduced. The data frame must include a column named \code{Genotype} that identifies the genotype associated with each row.}
  \item{nrow_lim}{A numeric value specifying the maximum number of rows desired in the reduced dataset. Defaults to \code{10e3}.}
  \item{min_samples}{An integer specifying the minimum number of samples to retain for each genotype. Defaults to \code{3}.}
}
\value{
  A reduced data frame that retains a subset of rows for each genotype group, ensuring that the total number of rows is approximately limited by \code{nrow_lim} and that each genotype has at least \code{min_samples} rows.
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
  df_reduced <- reduce_dataset_based_on_genotypes(df, nrow_lim = 100, min_samples = 5)
  
  # print the reduced dataset
  print(df_reduced)
}
