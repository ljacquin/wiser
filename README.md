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

In the ```wiser``` framework, the whitening matrix $W$ can be constructed using three different procedures: zero-phase component correlation analysis (```ZCA-cor```), principal component correlation analysis (```PCA-cor```) and cholesky (```Cholesky```). 
Let $\Sigma_u = \sigma^2_uZKZ'$ represent the genetic covariance matrix for the vector $U = Zu = (U_1,U_2, \ldots, U_n)$, which corresponds to the vector of individual genetic effects in the experimental design.
In  ```wiser```, the vector $u=(u_1,...,u_q)'$, representing the genetic values associated with $q$ genotypes, is assumed to follow the distribution $u \sim \mathcal{N}_q(0,\sigma^2_uK)$, where $K$ is the genetic covariance matrix between genotypes, estimated from omic data.
Define $R_u = V^{-\frac{1}{2}} \Sigma_u V^{-\frac{1}{2}}$ as the correlation matrix associated with $\Sigma_u$, where $V=diag(\Sigma_u)$ is a diagonal matrix with the variances of $\Sigma_u$ on its diagonal. The matrices $\Sigma_u$ and $R_u$, which are symmetric and positive semi-definite, have the following spectral decompositions: $\Sigma_u = U \Lambda U'$ and $R_u = G \Theta G'$, where $U$ and $G$ are orthogonal eigenvectors matrices for $\Sigma_u$ and $R_u$, respectively, while $\Lambda$ and $\Theta$ are diagonal matrices of positive eigenvalues for $\Sigma_u$ and $R_u$, respectively. 
The inverse square root matrices of $\Sigma_u$ and $R_u$ are given by $\Sigma_u^{-\frac{1}{2}} = U \Lambda^{-\frac{1}{2}} U'$ and $R_u^{-\frac{1}{2}} = G \Theta^{-\frac{1}{2}} G'$. These satisfy $\Sigma_u^{-1}=\Sigma_u^{-\frac{1}{2}}\Sigma_u^{-\frac{1}{2}}$ and  $R_u^{-1} = R_u^{-\frac{1}{2}}R_u^{-\frac{1}{2}}$. 
According to Kessy *et al.* (2015), the whitening matrices associated to  ```ZCA-cor```, ```PCA-cor``` and ```Cholesky``` are given by :
<br><br>

$$
\begin{cases}
  W^{ZCA-cor} = R_u^{-\frac{1}{2}}V^{-\frac{1}{2}} = G \Theta^{-\frac{1}{2}} G'V^{-\frac{1}{2}} = G \big[ \Theta^{-\frac{1}{2}} G'V^{-\frac{1}{2}} \big]\\      
  W^{PCA-cor} = \Theta^{-\frac{1}{2}} G'V^{-\frac{1}{2}}\\
  W^{Cholesky} = L^{-1} \\
\end{cases}
$$

where $L$ is derived from the Cholesky decomposition of $\Sigma_u = LL'$. The ```PCA-cor``` whitening procedure can be seen as standardizing variables using $V^{-\frac{1}{2}}$, followed by a rotation using the transposed correlation eigenmatrix $G'$, and then scaling using the inverted correlation singular values matrix $\Theta^{-\frac{1}{2}}$. ```ZCA-cor``` whitening extends this by applying an additional rotation $G$ to revert to the original basis of the standardized variables. Each whitening method is optimal according to specific criteria. For instance, ```ZCA-cor``` is unique in ensuring that the whitened variables retain the maximum correlation with the original variables. Details of these criteria and the optimality of each method are discussed in Kessy *et al.* (2015).

The rationale for transforming $X$ into $\tilde{X}$ through whitening, to adjust fixed effect variables for population structure, is comprehensively addressed in Jacquin *et al.* (year). In contrast, the chosen approach of successive ordinary least squares (OLS) in ```wiser``` avoids making assumptions about the properties of $\beta$ and $v$. Specifically, ```wiser``` does not assume that $v$ is a random vector drawn from a distribution with a specified covariance matrix. This approach prevents enforcing an unnecessary covariance structure during the estimation of $v$, which could be detrimental.
For example, assuming $v \sim \mathcal{N}_q(0,\sigma^2_vI_q)$ is often unrealistic and would lead to using a decorrelated covariance structure in the best linear unbiased predictor (BLUP) of $v$, which can be highly undesirable.
Crucially, the successive ordinary least squares (OLS) procedure implemented in ```wiser``` operates without such assumptions, ensuring that the estimation of $v$ remains independent of omic information or any imposed covariance structure. 
The only assumptions in the ```wiser``` framework are $u \sim \mathcal{N}_q(0,\sigma^2_u K)$ and $\varepsilon \sim \mathcal{N}_n(0,\sigma^2_ε I_n)$, which are necessary to estimate $\sigma^2_u$ for constructing the whitening matrix $W$. In this framework, the estimation of $\sigma^2_u$ and $\sigma^2_ε$ is performed using a parallelized ABC algorithm. 
In ```wiser```, two kernel functions are also provided to build $K$ : ```linear``` and ```identity```. The ```identity``` kernel, which is **not estimated from omic data**, is generally discouraged due to its poor performance associated with phenotypic predictive ability. The ```identity``` kernel is used solely to assess the impact of including or excluding the genetic covariance structure. The ```linear``` kernel is used by default in ```wiser```.

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

Here's a simple example illustrating the use of the ```estimate_wiser_phenotype``` function to estimate phenotypes :

```R
## ✅ load wiser library, display package help ,and attach data subsets from the refpop dataset

# load library 
library(wiser)

# show help for library 
help(package="wiser")
 
# load raw individual phenotypic measurements (for subset of 80 genotypes from refpop data)
data("refpop_raw_indiv_pheno_data_subset")
head(refpop_raw_indiv_pheno_data_subset)

# load SNP marker data (for subset of 80 genotypes from refpop data)
data("refpop_geno_data_subset")
head(refpop_geno_data_subset)[,1:10]

## ✅ estimate wiser phenotypes for a specified trait, i.e. Flowering_begin here

# 📌 define trait
trait_ <- "Flowering_begin"

# 📌 apply wiser estimation function using default values for `whitening_method` (set to "ZCA-cor") and `alpha_` (set to 0.01). These default values typically yield satisfactory results for phenotypic predictive ability. Therefore, using `optimize_whitening_and_regularization()` may not always be necessary, especially for large datasets. Nevertheless, these parameters should be optimized when low phenotype predictive ability is observed.
wiser_obj <- estimate_wiser_phenotype(
  omic_df = refpop_geno_data_subset,
  raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
  trait_ = trait_,
  fixed_effects_vars = c("Envir", "Country", "Year",
                         "Row", "Position", "Management"),
  fixed_effects_vars_computed_as_factor = c("Envir", "Country", "Year",
                                             "Row", "Position", "Management"),
  site_var = "Country",
  fixed_effects_vars_computed_as_factor_by_site = c("Row", "Position"),
  random_effects_vars = "Genotype",
  whitening_method = "ZCA-cor",
  alpha_ = 0.01
 )

# 📌 apply wiser estimation function using optimized whitening method and regularization parameter.
# ⚠️ highly recommended: increase memory size as follows before using optimize_whitening_and_regularization()
options(future.globals.maxSize = 16 * 1024^3)

opt_white_reg_par <- optimize_whitening_and_regularization(
  omic_df = refpop_geno_data_subset,
  raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
  trait_ = trait_,
  whitening_method_grid = c("ZCA-cor","Cholesky"),
  k_folds_ = 3,
  alpha_grid = c(0.01, 0.1)
)

wiser_obj <- estimate_wiser_phenotype(
  omic_df = refpop_geno_data_subset,
  raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
  trait_ = trait_,
  fixed_effects_vars = c("Envir", "Country", "Year",
                         "Row", "Position", "Management"),
  fixed_effects_vars_computed_as_factor = c("Envir", "Country", "Year",
                                             "Row", "Position", "Management"),
  site_var = "Country",
  fixed_effects_vars_computed_as_factor_by_site = c("Row", "Position"),
  random_effects_vars = "Genotype",
  whitening_method = as.character(opt_white_reg_par$opt_whitening_method),
  alpha_ = as.numeric(opt_white_reg_par$opt_alpha_)
 )

## ✅ plot and print wiser objects

# plot the density for the estimated phenotypes
dev.new()
plot(density(wiser_obj$wiser_phenotypes$v_hat), main = paste0(trait_, " v_hat"))

# get the estimated fixed effects from the whitening process and OLS, and variance components from ABC
# estimated fixed effects (from whitening and OLS)
print(wiser_obj$wiser_fixed_effect_estimates)

# estimated variance components (from ABC)
print(wiser_obj$wiser_abc_variance_component_estimates)

## ✅ print fundamental property of whitening, i.e. In = W*Σu*W'
id_mat <- wiser_obj$w_mat %*% wiser_obj$sig_mat_u %*% t(wiser_obj$w_mat)
print(id_mat[1:10,1:10])

```

## Authors and References

* Author : Laval Jacquin
* Maintainer : Laval Jacquin jacquin.julien@gmail.com

## References
Kessy, A., Lewin, A., & Strimmer, K. (2018). Optimal whitening and decorrelation. The American Statistician, 72(4), 309-314.

