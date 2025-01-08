\name{compute_incidence_matrices_fixed_and_random_effects}
\alias{compute_incidence_matrices_fixed_and_random_effects}
\title{Compute incidence matrices for fixed and random effects}
\description{
  This function computes incidence matrices for both fixed and random effects based on the specified variables. It creates a matrix for each fixed effect, including an intercept column for the first fixed effect if applicable, and computes matrices for random effects without intercepts.
}
\usage{
  compute_incidence_matrices_fixed_and_random_effects(
    fixed_effect_vars,
    fixed_effect_vars_computed_as_factor,
    random_effects_vars,
    raw_pheno_df
  )
}
\arguments{
  \item{fixed_effect_vars}{A character vector of column names in \code{raw_pheno_df} that represent fixed effects variables.}
  \item{fixed_effect_vars_computed_as_factor}{A character vector of column names in \code{raw_pheno_df} that should be treated as factors when creating incidence matrices.}
  \item{random_effects_vars}{A character vector of column names in \code{raw_pheno_df} that represent random effects variables.}
  \item{raw_pheno_df}{A data frame containing the raw phenotypic data. The columns specified in \code{fixed_effect_vars} and \code{random_effects_vars} must be present in this data frame.}
}
\value{
  A list with two elements:
    \itemize{
      \item \code{x_mat}: A numeric matrix representing the incidence matrices for fixed effects, including an intercept column for the first fixed effect if specified.
      \item \code{z_mat}: A numeric matrix representing the incidence matrices for random effects.
    }
}
\details{
  The function constructs incidence matrices for fixed effects and random effects as follows:
    \itemize{
      \item For fixed effects variables, it includes an intercept column for the first fixed effect. If the first fixed effect is treated as a factor, it creates dummy variables for it. For other fixed effects, it either includes them directly or creates dummy variables if they are specified as factors.
      \item For random effects variables, it creates incidence matrices without an intercept column.
    }
  
  The resulting matrices are combined into single matrices for fixed effects (\code{x_mat}) and random effects (\code{z_mat}). Both matrices are converted to numeric format.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  library(wiser)
  
  # load required data
  data("refpop_raw_indiv_pheno_data_subset")
  
  # define fixed and random effects variables
  fixed_effect_vars <- c("Envir", "Country", "Year", "Row", "Position", "Management")
  fixed_effect_vars_computed_as_factor <- c("Envir", "Country", "Year", "Row", "Position", "Management")
  random_effects_vars <- c("Genotype")
  
  # compute incidence matrices
  incidence_matrices <- compute_incidence_matrices_fixed_and_random_effects(
    fixed_effect_vars = fixed_effect_vars,
    fixed_effect_vars_computed_as_factor = fixed_effect_vars_computed_as_factor,
    random_effects_vars = random_effects_vars,
    raw_pheno_df = refpop_raw_indiv_pheno_data_subset
  )
  
  # print the incidence matrices
  print(incidence_matrices$x_mat)
  print(incidence_matrices$z_mat)
}