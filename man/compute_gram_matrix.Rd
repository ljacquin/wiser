\name{compute_gram_matrix}
\alias{compute_gram_matrix}
\title{Compute the Gram matrix}
\description{
  This function computes the Gram matrix based on the given omic data. It supports different kernel types and ensures the resulting matrix is positive definite.
}
\usage{
  compute_gram_matrix(omic_df, kernel_type)
}
\arguments{
  \item{omic_df}{A data frame containing the omic data (e.g. genotypes as rows and SNPs, metabolites, or wavelength reflectances as columns), the row names of \code{omic_df} should be genotype identifiers.}

  \item{kernel_type}{A character string specifying the type of kernel to be used. Options are \code{"linear"} or \code{"identity"} (which is not recommended).}
}
\value{
  A matrix representing the Gram matrix. The matrix is positive definite, with genotype row names and column names assigned.
}
\details{
  The function scales the omic data and computes the Gram matrix based on the specified kernel type:
    \itemize{
      \item \code{"linear"}: Computes a linear kernel Gram matrix based on centered omic data.
      \item \code{"identity"}: Uses an identity matrix, which is not recommended due to its limited precision in capturing genotype relationships.
    }
  
  The function checks if the resulting matrix is positive definite and, if necessary, adjusts it to ensure positive definiteness using the \code{\link[Matrix]{nearPD}} function from the \code{Matrix} package.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage

  # load wiser library and example data subsets from the refpop dataset
  library(wiser)

  # load SNP marker data (for a subset of 80 genotypes from refpop data)
  data("refpop_geno_data_subset")
  head(refpop_geno_data_subset)[,1:10]
  
  # compute the Gram matrix with a linear kernel
  gram_matrix_linear <- compute_gram_matrix(
    omic_df = refpop_geno_data_subset,
    kernel_type = "linear"
  )
  
  # print the computed Gram matrix
  print(gram_matrix_linear)
  
}
