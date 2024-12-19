[<img src="img/wiser.png"/>]()

# WISER: WhItening and successive least Squares Estimation Refinement for phenotype estimation

##### Licence, status and metrics
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)]()
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub repo size](https://img.shields.io/github/repo-size/ljacquin/wiser)
![GitHub language count](https://img.shields.io/github/languages/count/ljacquin/wiser)
![GitHub top language](https://img.shields.io/github/languages/top/ljacquin/wiser)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/ljacquin/wiser)  
![GitHub all releases](https://img.shields.io/github/downloads/ljacquin/wiser/total)
![GitHub stars](https://img.shields.io/github/stars/ljacquin/wiser)  

##### Languages and technologies
[![R Badge](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://www.r-project.org/)

## Overview

The ```wiser``` package offers *user-friendly* advanced tools for more accurate phenotypic estimation, by leveraging genetic covariance information to correct for population structure in experimental designs. By employing a whitening transformation followed by successive ordinary least squares (OLS) estimation, ```wiser``` refines fixed effect estimates and eliminates biases, leading to improved phenotypic estimation. This approach is particularly beneficial in complex experimental designs where genetic and environmental factors are intricately linked. ```wiser``` includes methods for computing whitening matrices, fixed effects, residuals, and estimating phenotypes. Additionally, the package provides fast and stable variance component estimation using a parallelized approximate Bayesian computation (ABC) algorithm, specifically designed for large datasets associated with complex experimental designs.

For an experimental design, ```wiser``` solves the following model, in order to estimate a vector $v=(v_1,...,v_q)'$ of $q$ phenotypes which approximate the genetic values associated to $q$ genotypes: 
<br><br>

$$
\\
Y = WX\beta + Zv + \varepsilon = \tilde{X}\beta + \xi 
\\
$$

where :

* $Y \ (n \times 1)$ is the vector of individual phenotypic measurements, with values repeated for each genotype.

* $W \ (n \times n)$ is a whitening matrix built using omic data (e.g., SNP markers, metabolites, or wavelength reflectances), which corrects the fixed effects variables for population structure.

* $X \ (n \times l)$ is the incidence matrix linking fixed effects to individual phenotypic measurements.

* $\beta \ (l \times 1)$ is the vector of fixed effects.

* $v \ (q \times 1)$ corresponds to a vector of $q$ phenotypes approximating the genetic values associated with $q$ genotypes.

* $Z \ (n \times q)$ corresponds to the design matrix which links the q phenotypes to the individual phenotypic values in the experimental design.

* $\varepsilon \ (n \times 1)$ is the vector of residuals.

* $\tilde{X} = X W$ and $\xi = Zv + \varepsilon$.

<br>
In the ```wiser``` framework, the whitening matrix $W$ can be constructed using three different procedures: zero-phase component correlation analysis (```ZCA-cor```), principal component correlation analysis (```PCA-cor```) and cholesky (```Cholesky```). 
Let $\Sigma_u = \sigma^2_uZKZ'$ represent the genetic covariance matrix for the vector $U = Zu = (U_1,U_2, \ldots, U_n)$, which corresponds to the vector of individual genetic effects in the experimental design.
In  ```wiser```, the vector $u=(u_1,...,u_q)'$, representing the genetic values associated with $q$ genotypes, is assumed to follow the distribution $u \sim \mathcal{N}_q(0,\sigma^2_uK)$, where $K$ is the genetic covariance matrix between genotypes, estimated from omic data.<br>
Define $R_u = V^{-\frac{1}{2}} \Sigma_u V^{-\frac{1}{2}}$ as the correlation matrix associated with $\Sigma_u$, where $V=diag(\Sigma_u)$ is a diagonal matrix with the variances of $\Sigma_u$ on its diagonal. The matrices $\Sigma_u$ and $R_u$, which are symmetric and positive semi-definite, have the following spectral decompositions: $\Sigma_u = U \Lambda U'$ and $R_u = G \Theta G'$, where $U$ and $G$ are orthogonal eigenvectors matrices for $\Sigma_u$ and $R_u$, respectively, while $\Lambda$ and $\Theta$ are diagonal matrices of positive eigenvalues for $\Sigma_u$ and $R_u$, respectively. <br>
The inverse square root matrices of $\Sigma_u$ and $R_u$ are given by $\Sigma_u^{-\frac{1}{2}} = U \Lambda^{-\frac{1}{2}} U'$ and $R_u^{-\frac{1}{2}} = G \Theta^{-\frac{1}{2}} G'$. These satisfy $\Sigma_u^{-1}=\Sigma_u^{-\frac{1}{2}}\Sigma_u^{-\frac{1}{2}}$ and  $R_u^{-1} = R_u^{-\frac{1}{2}}R_u^{-\frac{1}{2}}$. 
According to Kessy $\textit{et al.}$ (2015), the whitening matrices associated to  ```ZCA-cor```, ```PCA-cor``` and ```Cholesky``` are given by :
<br><br>

$$
\begin{cases}
  W^{ZCA-cor} = R_u^{-\frac{1}{2}}V^{-\frac{1}{2}} = G \Theta^{-\frac{1}{2}} G'V^{-\frac{1}{2}} = G \big[ \Theta^{-\frac{1}{2}} G'V^{-\frac{1}{2}} \big]\\      
  W^{PCA-cor} = \Theta^{-\frac{1}{2}} G'V^{-\frac{1}{2}}\\
  W^{Cholesky} = L^{-1} \\
\end{cases}
$$

<br>
where $L$ is derived from the Cholesky decomposition of $\Sigma_u = LL'$. The ```PCA-cor``` whitening procedure can be seen as standardizing variables using $V^{-\frac{1}{2}}$, followed by a rotation using the transposed correlation eigenmatrix $G'$, and then scaling using the inverted correlation singular values matrix $\Theta^{-\frac{1}{2}}$. ```ZCA-cor``` whitening extends this by applying an additional rotation $G$ to revert to the original basis of the standardized variables. Each whitening method is optimal according to specific criteria. For instance, ```ZCA-cor``` is unique in ensuring that the whitened variables retain the maximum correlation with the original variables. Details of these criteria and the optimality of each method are discussed in Kessy $\textit{et al.}$ (2015).

The rationale for transforming $X$ into $\tilde{X}$ through whitening, to adjust fixed effect variables for population structure, is comprehensively addressed in Jacquin $\textit{et al.}$ (2025). In contrast, the chosen approach of successive ordinary least squares (OLS) in ```wiser``` avoids making assumptions about the properties of $\beta$ and $v$. Specifically, ```wiser``` does not assume that $v$ is a random vector drawn from a distribution with a specified covariance matrix. This approach prevents enforcing an unnecessary covariance structure during the estimation of $v$, which could be detrimental. <br>
For example, assuming $v \sim \mathcal{N}_q(0,\sigma^2_vI_q)$ is often unrealistic and would lead to using a decorrelated covariance structure in the best linear unbiased predictor (BLUP) of $v$, which can be highly undesirable. <br>
Crucially, the successive ordinary least squares (OLS) procedure implemented in ```wiser``` operates without such assumptions, ensuring that the estimation of $v$ remains independent of omic information or any imposed covariance structure. The only assumptions in the ```wiser``` framework are $u \sim \mathcal{N}_q(0,\sigma^2_u K)$ and $\varepsilon \sim \mathcal{N}_n(0,\sigma^2_\varepsilon I_n)$, which are necessary to estimate 
$\sigma^2_u$ for constructing the whitening matrix $W$. In this framework, the estimation of $\sigma^2_u$ and $\sigma^2_\varepsilon$ is performed using an ABC algorithm. 

In ```wiser```, two kernel functions are also provided to build $K$: ```linear``` and ```identity```. The ```identity``` kernel, which is **not estimated from omic data**, is generally discouraged due to its poor phenotypic predictive ability. The ```identity``` kernel is used solely to assess the impact of including or excluding the genetic covariance structure. The ```linear``` kernel is used by default in ```wiser```.

## Installation

You can install the latest version of the ```wiser``` package with :

```R
install.packages("devtools")
library(devtools)
install_github("ljacquin/wiser")
```

## Key Features

    ▸ Phenotype estimation: estimate phenotypic values that approximate genetic values, by applying whitening methods to correct for genetic covariance structure in experimental designs (i.e., population structure correction).
    ▸ Whitening methods: implement various whitening techniques, such as ZCA-cor, PCA-cor, and Cholesky, to effectively correct for population structure.
    ▸ Genetic covariance matrix regularization: ensure the stability of genetic covariance matrices, by applying regularization and shrinkage techniques which ensure positive definiteness.
    ▸ Optimal whitening and regularization: automatically determine the best whitening method and regularization parameter optimizing phenotypic predictive ability, through k-fold cross-validation using a subset of the data.
    ▸ Fast and stable variance component estimation: use a parallelized ABC algorithm to achieve fast and stable variance component estimation in large datasets associated to complex experimental designs.
  
## Main Functions

    ▸ estimate_wiser_phenotype: estimates phenotypic values approximating genetic values using whitening methods.
    ▸ optimize_whitening_and_regularization: finds the optimal combination of whitening method and regularization parameter through cross-validation for phenotype prediction.

## Example

### Estimating phenotypes with WISER

Here are straightforward examples demonstrating how to use the ```estimate_wiser_phenotype``` function for phenotype estimation:

```R
# ➡️ load wiser library, display package help ,and attach datasets associated to four species

# -- load library
library(wiser)

# -- show help for library
help(package = "wiser")

# -- attach datasets for random samples of 30 genotypes associated to apple, pine, maize
# and rice

# 📌 These datasets are small subsets derived from the original datasets used in Jacquin et al. (2025). 
# They are provided for illustrative purposes with the wiser package and are not intended to represent reference populations for genomic prediction or GWAS. These datasets include genomic data and raw individual phenotypic
# measurements for 30 randomly selected genotypes, associated with an experimental design for each of the 
# following species: apple, pine, maize and rice.

# apple
data("apple_raw_pheno_data")
data("apple_genomic_data")
# pine
data("pine_raw_pheno_data")
data("pine_genomic_data")
# maize
data("maize_raw_pheno_data")
data("maize_genomic_data")
# rice
data("rice_raw_pheno_data")
data("rice_genomic_data")

# ➡️ usage examples with apple data

# -- display apple data 
head(apple_raw_pheno_data)
head(apple_genomic_data)[, 1:10]

# -- for raw individual phenotypic data: create environment variable (combination 
# of country, year and management) 
apple_raw_pheno_data$Envir <- paste0(
  apple_raw_pheno_data$Country, "_",
  apple_raw_pheno_data$Year, "_",
  apple_raw_pheno_data$Management
)

# -- for genomic data: assign genotypes as rownames to the data (if necessary) to be
# compliant with wiser functions, and remove "Genotype" column (if present) 
rownames(apple_genomic_data) <- apple_genomic_data$Genotype
apple_genomic_data <- apple_genomic_data[, -match(
  "Genotype",
  colnames(apple_genomic_data)
)]

# -- define a trait for phenotypic estimation using wiser (note: some data should be
# available for the trait)
trait_ <- "Trunk_increment"

# -- apply wiser estimation function using default values for `whitening_method` (set to "ZCA-cor")
# and `alpha_` (set to 0.01)

# 📌 These default values typically yield satisfactory results for phenotypic predictive ability. Therefore,
# using `optimize_whitening_and_regularization()` may not always be necessary, especially for large datasets. 
# Nevertheless, these parameters should be optimized, when possible, for better results.

# 📌 Note that rows and positions are fitted as factors by environment to account for local spatial heterogeneity.

wiser_obj <- estimate_wiser_phenotype(
  omic_df = apple_genomic_data,
  raw_pheno_df = apple_raw_pheno_data,
  trait_ = trait_,
  fixed_effects_vars = c(
    "Envir", "Row", "Position"
  ),
  fixed_effects_vars_computed_as_factor = c(
    "Envir", "Row", "Position"
  ),
  envir_var = "Envir",
  fixed_effects_vars_computed_as_factor_by_envir = c("Row", "Position"),
  random_effects_vars = "Genotype"
)

# -- apply the WISER estimation function using the optimized whitening method
# and regularization parameter.

# 📌⚠️ highly recommended: increase memory size as specified below with options() before using optimize_whitening_and_regularization(). For optimal performance, it is strongly advised to use
# a high-performance computing cluster (HPC) when running this function.

run_example <- F
if (run_example){
  options(future.globals.maxSize = 16 * 1024^3)
  
  opt_white_reg_par <- optimize_whitening_and_regularization(
    omic_df = apple_genomic_data,
    raw_pheno_df = apple_raw_pheno_data,
    trait_ = trait_,
    whitening_method_grid = c("ZCA-cor", "Cholesky"),
    k_folds_ = 3,
    alpha_grid = c(0.01, 0.1)
  )
  print(opt_white_reg_par)
  opt_alpha_par_ <- as.numeric(opt_white_reg_par$opt_alpha_)
  opt_white_method_ <- as.character(opt_white_reg_par$opt_whitening_method)
  
  wiser_obj <- estimate_wiser_phenotype(
    omic_df = apple_genomic_data,
    raw_pheno_df = apple_raw_pheno_data,
    trait_ = trait_,
    fixed_effects_vars = c("Envir", "Row", "Position"),
    fixed_effects_vars_computed_as_factor = c("Envir", "Row", "Position"),
    envir_var = "Envir",
    fixed_effects_vars_computed_as_factor_by_envir = c("Row", "Position"),
    random_effects_vars = "Genotype",
    whitening_method = opt_white_method_,
    alpha_ = opt_alpha_par_
  )
}

# -- plot the density for the estimated phenotypes
dev.new()
plot(density(wiser_obj$wiser_phenotypes$v_hat), main = paste0(trait_, " v_hat"))

# get the estimated fixed effects from the whitening process and OLS, and variance components from ABC
# estimated fixed effects (from whitening and OLS)
print(wiser_obj$wiser_fixed_effect_estimates)

# estimated variance components (from ABC)
print(wiser_obj$wiser_abc_variance_component_estimates)

# print fundamental property of whitening, i.e. In = W*Σu*W'
id_mat <- wiser_obj$w_mat %*% wiser_obj$sig_mat_u %*% t(wiser_obj$w_mat)
print(id_mat[1:5,1:5])

# ➡️ usage examples with pine data

# -- display pine data
head(pine_raw_pheno_data)
head(pine_genomic_data)[, 1:10]

# -- for raw individual phenotypic data: create environment variable (combination
# of site, year and block)
pine_raw_pheno_data$Envir <- paste0(
  pine_raw_pheno_data$Site, "_",
  pine_raw_pheno_data$Year, "_",
  pine_raw_pheno_data$Block
)

# -- generate latitude and longitude variables per environment
pine_raw_pheno_data <- generate_latitude_longitude_variables_by_environment(
  pine_raw_pheno_data
)

# -- for genomic data: assign genotypes as rownames to the data (if necessary) to
# be compliant with wiser functions, and remove "Genotype" column (if present),
# here "V1" is "Genotype" column
rownames(pine_genomic_data) <- pine_genomic_data$V1
pine_genomic_data <- pine_genomic_data[, -match(
  "V1", colnames(pine_genomic_data)
)]

# -- define a trait for phenotypic estimation using wiser (note: some data should be
# available for the trait)
trait_ <- "H" # height

# -- get fixed effect vars where latitude and longitude are fitted as quantitative
# variables for each environment (i.e. combination of site, year and block). 

# 📌 Note that the latitude and longitude fixed-effect variables, which are quantitative and not
# considered as factors, are highly correlated to each environment for which they are fitted. 
# Therefore, the environment variable will not be fitted using wiser due to redundancy.

fixed_effect_vars_ <- grep("_latitude$|_longitude$", colnames(pine_raw_pheno_data),
  value = TRUE
)

# -- estimate wiser phenotype
start_time_ <- Sys.time()
wiser_obj <- estimate_wiser_phenotype(
  pine_genomic_data,
  pine_raw_pheno_data, 
  trait_,
  fixed_effects_vars = fixed_effect_vars_,
  fixed_effects_vars_computed_as_factor = NULL,
  envir_var = NULL,
  fixed_effects_vars_computed_as_factor_by_envir = NULL,
  random_effects_vars = "Genotype"
  )

# -- plot and print wiser objects

# plot the density for the estimated phenotypes
dev.new()
plot(density(wiser_obj$wiser_phenotypes$v_hat), main = paste0(trait_, " v_hat"))

# get the estimated fixed effects from the whitening process and OLS, and variance components from ABC
# estimated fixed effects (from whitening and OLS)
print(wiser_obj$wiser_fixed_effect_estimates)

# estimated variance components (from ABC)
print(wiser_obj$wiser_abc_variance_component_estimates)

# print fundamental property of whitening, i.e. In = W*Σu*W'
id_mat <- wiser_obj$w_mat %*% wiser_obj$sig_mat_u %*% t(wiser_obj$w_mat)
print(id_mat[1:5, 1:5])

# ➡️ usage examples with maize data


```

## Authors and References

* Author : Laval Jacquin
* Maintainer : Laval Jacquin jacquin.julien@gmail.com

## References
Kessy, A., Lewin, A., & Strimmer, K. (2018). Optimal whitening and decorrelation. The American Statistician, 72(4), 309-314.

