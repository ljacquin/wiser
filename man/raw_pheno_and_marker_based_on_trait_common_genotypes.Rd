\name{raw_pheno_and_marker_based_on_trait_common_genotypes}
\alias{raw_pheno_and_marker_based_on_trait_common_genotypes}
\title{Filter phenotypic and omic data based on common genotypes for a specific trait}
\description{
  This function filters the raw phenotypic data and the omic data based on common genotypes for a specified trait. It removes rows with missing values for the trait, retains only the relevant variables, and ensures that both datasets include only genotypes that are common between them.
}
\usage{
  raw_pheno_and_marker_based_on_trait_common_genotypes(
    raw_pheno_df,
    omic_df,
    trait_,
    fixed_effect_vars,
    random_effect_vars
  )
}
\arguments{
  \item{raw_pheno_df}{A data frame containing the raw phenotypic data. Each row should represent an observation with multiple columns, including the trait of interest.}
  \item{omic_df}{A data frame containing the omic data (e.g. genotypes as rows and SNPs, or metabolites, or wavelength reflectances/absorbances as columns).}
  \item{trait_}{A character string specifying the trait of interest to filter the phenotypic data.}
  \item{fixed_effect_vars}{A character vector listing the fixed effects variables to be retained in the filtered phenotypic data.}
  \item{random_effect_vars}{A character vector listing the random effects variables to be retained in the filtered phenotypic data.}
}
\value{
  A list with two elements:
    \item{raw_pheno_df}{A data frame containing the filtered phenotypic data, including only the rows with non-missing values for the specified trait and retaining only the relevant columns.}
  \item{omic_df}{A data frame containing the filtered omic data, including only the genotypes present in the filtered phenotypic data.}
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  
  # load example data
  data("example_raw_pheno_df")
  data("example_omic_df")
  
  # define parameters
  trait <- "Flowering_begin"
  fixed_effect_vars <- c("Envir", "Country", "Year", "Row", "Position", "Management")
  random_effect_vars <- "Genotype"
  
  # filter raw phenotypic and genomic data based on common genotypes
  filtered_data <- raw_pheno_and_marker_based_on_trait_common_genotypes(
    raw_pheno_df = example_raw_pheno_df,
    omic_df = example_omic_df,
    trait_ = trait,
    fixed_effect_vars = fixed_effect_vars,
    random_effect_vars = random_effect_vars
  )
  
  # print dimensions of the filtered datasets
  print(dim(filtered_data$raw_pheno_df))
  print(dim(filtered_data$omic_df))
}
