\name{compute_fixed_effect_vars_declared_as_factors}
\alias{compute_fixed_effect_vars_declared_as_factors}
\title{Convert specified fixed effects variables to factors}
\description{
  This function converts specified fixed effects variables in a raw phenotypic data frame to factors. It handles variables declared as factors differently based on whether they are also declared to be computed as factors by environment. 
}
\usage{
  compute_fixed_effect_vars_declared_as_factors(
    raw_pheno_df,
    fixed_effect_vars,
    fixed_effect_vars_computed_as_factor,
    envir_var,
    fixed_effect_vars_computed_as_factor_by_envir
  )
}
\arguments{
  \item{raw_pheno_df}{A data frame containing the raw phenotypic data.}
  \item{fixed_effect_vars}{A character vector listing all fixed effects variables in the data frame.}
  \item{fixed_effect_vars_computed_as_factor}{A character vector specifying which fixed effects variables should be converted to factors.}
  \item{envir_var}{A character string specifying the environment variable to be used for creating environment-specific factors.}
  \item{fixed_effect_vars_computed_as_factor_by_envir}{A character vector listing the fixed effects variables that should be converted to factors by environment.}
}
\value{
  A data frame with the fixed effects variables converted to factors according to the specified criteria. The levels of factors that are not present in the data are also dropped.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  
  # load wiser library and example data subsets from the refpop dataset
  library(wiser)

  # load raw phenotypic measurements (for a subset of 80 genotypes from refpop data)
  data("refpop_raw_indiv_pheno_data_subset")
  head(refpop_raw_indiv_pheno_data_subset)
  
  # define parameters
  fixed_effect_vars <- c("Envir", "Country", "Year", "Row", "Position", "Management")
  fixed_effect_vars_computed_as_factor <- c("Envir", "Country", "Year", "Row", "Position", "Management")
  envir_var <- "Country"
  fixed_effect_vars_computed_as_factor_by_envir <- c("Row", "Position")
  
  # convert fixed effects variables to factors
  processed_data <- compute_fixed_effect_vars_declared_as_factors(
    raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
    fixed_effect_vars = fixed_effect_vars,
    fixed_effect_vars_computed_as_factor = fixed_effect_vars_computed_as_factor,
    envir_var = envir_var,
    fixed_effect_vars_computed_as_factor_by_envir = fixed_effect_vars_computed_as_factor_by_envir
  )
  
  # print the structure of the processed data
  str(processed_data)
}
