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

The ```wiser``` package offers *user-friendly* advanced tools for more accurate phenotypic estimation, by leveraging genetic covariance information in the estimation of fixed effects. By employing a whitening transformation followed by successive ordinary least squares (OLS) estimation, ```wiser``` refines fixed effect estimates and eliminates biases, leading to improved phenotypic estimation. This approach is particularly beneficial in complex experimental designs where genetic and environmental factors are intricately linked. ```wiser``` includes methods for computing whitening matrices, fixed effects, residuals, and estimating phenotypes. Additionally, the package offers robust variance component estimation using approximate bayesian computation (ABC), ensuring stability and reliability of estimations. 

The phenotypic estimation performed by `wiser` for a vector $v = (v_1, v_2, \ldots, v_q)'$ of $q$ phenotypes, which approximates a vector $u = (u_1, u_2, \ldots, u_q)'$ of genetic values for $q$ genotypes, is carried out using ordinary least squares (OLS) as follows :

$$
\hat{v} = \underset{v \in \mathbb{R}^q}{argmin} || \hat{\xi} - Zv||^2_2 = (Z'Z)^{-1}Z'\hat{\xi} \ \ \ \
$$ 

where :

* $Z$ is the incidence matrix that links the $q$ phenotypes, which approximate the genetic values of $q$ genotypes, to the individual phenotypic measurements repeated for each genotype in the experimental design.

* $\overset{\wedge}{\xi}$ is the vector of estimated residuals obtained after removing the fixed effects, adjusted for the genetic covariance between individuals in the experimental design.

In this setting, the vector $\overset{\wedge}{\xi}$ is estimated as follows:

$$
\hat{\xi} = Y - \tilde{X}\hat{\beta} = Y - WX\hat{\beta} \ \ \ \
$$ 

where :

* $Y = (Y_1,Y_2, \dots, Y_n)'$ is the vector of $n$ individual phenotypic measurements in the experimental design, with repeated measurements for each genotype.

* $\tilde{X} (= WX)$ is the transformed (whitened) design matrix that links the estimated fixed effects, adjusted for genetic covariance in the experimental design, to the individual phenotypic measurements.

* $\hat{\beta} = (\hat{\beta}_1, \hat{\beta}_2, \ldots ,\hat{\beta}_l)'= (\tilde{X}'\tilde{X})^{-1}\tilde{X}'Y$ is the OLS vector of $l$ estimated fixed effects, adjusted for genetic covariance between individuals in the experimental design.
  
* $W$ is an $n$ x $n$ whitening matrix built using omic data (e.g. SNP markers, metabolites, or wavelength reflectances), which eliminates the genetic covariance structure between individuals during the fixed effects estimation process.

* $X$ is the original design matrix that links fixed effects to the individual phenotypic measurements in the experimental design.

In the ```wiser``` framework, the whitening matrix $W$ can be constructed using three different procedures: zero-phase component correlation analysis (```ZCA-cor```), principal component correlation analysis (```PCA-cor```) and cholesky (```Cholesky```). Let $\Sigma_u = \sigma^2_uZKZ'$ represent the genetic covariance matrix for the vector $U = Zu = (U_1,U_2, \ldots, U_n)$, which corresponds to the vector of individual genetic effects in the experimental design. It is assumed that $u \sim \mathcal{N}_q(0,\sigma^2_uK)$, where $K$ is the genetic covariance matrix between genotypes estimated from omic data. Define $R_u = V^{-\frac{1}{2}} \Sigma_u V^{-\frac{1}{2}}$ as the correlation matrix associated with $\Sigma_u$, where $V=diag(\Sigma_u)$ is a diagonal matrix with the variances of $\Sigma_u$ on its diagonal. Since $\Sigma_u$ and $R_u$ are symmetric and positive semi-definite, they have the following spectral decompositions: $\Sigma_u = U \Lambda U'$ and $R_u = G \Theta G'$, where $U$ and $G$ are orthogonal eigenvectors matrices for $\Sigma_u$ and $R_u$, respectively (i.e. $UU'=U'U=I_n$ and $GG'=G'G=I_n$), while $\Lambda$ and $\Theta$ are diagonal matrices of positive eigenvalues for $\Sigma_u$ and $R_u$, respectively. The inverse square root matrices of $\Sigma_u$ and $R_u$ are given by $\Sigma_u^{-\frac{1}{2}} = U \Lambda^{-\frac{1}{2}} U'$ and $R_u^{-\frac{1}{2}} = G \Theta^{-\frac{1}{2}} G'$. These satisfy $\Sigma_u^{-1}=\Sigma_u^{-\frac{1}{2}}\Sigma_u^{-\frac{1}{2}}$ and  $R_u^{-1} = R_u^{-\frac{1}{2}}R_u^{-\frac{1}{2}}$). According to Kessy $\textit{et al.}$ (2015), the whitening matrices associated to  ```ZCA-cor```, ```PCA-cor``` and ```Cholesky``` are given by :

$$
\begin{cases}
  W^{ZCA-cor} = R_u^{-\frac{1}{2}}V^{-\frac{1}{2}} = G \Theta^{-\frac{1}{2}} G'V^{-\frac{1}{2}} = G \big[ \Theta^{-\frac{1}{2}} G'V^{-\frac{1}{2}} \big]\\      
  W^{PCA-cor} = \Theta^{-\frac{1}{2}} G'V^{-\frac{1}{2}}\\
  W^{Cholesky} = L^{-1} \\
\end{cases}
$$

where $L$ is derived from the Cholesky decomposition of $\Sigma_u = LL'$. The ```PCA-cor``` whitening procedure can be seen as standardizing variables using $V^{-\frac{1}{2}}$, followed by a rotation using the transposed correlation eigenmatrix $G'$, and then scaling using the inverted correlation singular values matrix $\Theta^{-\frac{1}{2}}$. ```ZCA-cor``` whitening extends this by applying an additional rotation $G$ to revert to the original basis of the standardized variables. Each whitening method is optimal according to specific criteria. For instance, ```ZCA-cor``` is unique in ensuring that the whitened variables retain the maximum correlation with the original variables. Details of these criteria and the optimality of each method are discussed in Kessy $\textit{et al.}$ (2015).

In ```wiser```, two kernel functions are also provided to build $K$: ```linear``` and ```identity```. The ```identity``` kernel is generally discouraged due to its poor performance in genomic prediction. On the other hand, the ```linear``` kernel has demonstrated superior practical performance in this context. The ```identity``` kernel is mainly used for benchmarking, as detailed in Jacquin $\textit{et al.}$ (2024), to assess the impact of including versus excluding the covariance structure.

#### Remarks :

The rationale for transforming $X$ into $\tilde{X}$ through whitening, to adjust for genetic covariance among individuals when estimating fixed effects, is comprehensively addressed in Jacquin $\textit{et al.}$ (2024). In contrast, the chosen approach of successive ordinary least squares (OLS) in ```wiser``` avoids making assumptions about the properties of $\beta$ and $v$. Specifically, ```wiser``` does not assume that $v$ is a random vector drawn from a distribution with a specified covariance matrix. This approach prevents enforcing an unnecessary covariance structure during the estimation of $v$, which could be detrimental. For example, assuming $v \sim \mathcal{N}_q(0,\sigma^2_vI_q)$ is often unrealistic and would lead to using a decorrelated covariance structure in the best linear unbiased predictor (BLUP) of $v$. Furthermore, and crucially, the successive OLS procedure ensures that the estimated residuals or $v$ are not linked to the omic data as shown in Jacquin $\textit{et al.}$ (2024).

## Installation

You can install the latest version of the ```wiser``` package with :

```R
install.packages("devtools")
library(devtools)
install_github("ljacquin/wiser")
```

## Key Features

    ▸ Phenotype Estimation: Estimate phenotypic values that closely approximate genetic values by leveraging advanced whitening techniques. This approach enhances the modeling of both genetic and environmental effects, as well as their interactions, leading to more precise and reliable phenotype estimation for each genotype.
    ▸ Whitening Methods: Implement advanced techniques like ZCA-cor, PCA-cor, and Cholesky whitening to adjust for genetic covariance in fixed effects estimation. This process effectively accounts for genetic relationships between individuals in the experimental design.
    ▸ Fixed and Random Effects Modeling: Transform fixed effects variables using whitening techniques and obtain OLS fixed effect estimates, which fully accounting for genetic covariance in the experimental design.    
    ▸ Optimal Whitening and Regularization: Automatically select the best whitening method and regularization parameter using k-fold cross-validation for precise phenotype estimation.
    ▸ Genetic Covariance Matrix Regularization: Stabilize genetic covariance matrices through innovative regularization methods, ensuring they remain positive definite.
    ▸ Efficient Variance Component Estimation: Use Approximate Bayesian Computation (ABC) to accurately estimate variance components in complex experimental designs.
    ▸ Phenotype Simulation and Distance Metrics: Simulate phenotypic data based on estimated model parameters and compute distances like the squared L2 norm to compare observed and simulated outcomes.

## Main Functions

    ▸ estimate_wiser_phenotype: Estimates phenotypic values approximating genetic values using whitening methods.
    ▸ optimize_whitening_and_regularization: Finds the optimal combination of whitening method and regularization parameter through cross-validation for phenotype prediction.
    ▸ compute_transformed_vars_and_ols_estimates: Computes transformed (whitened) fixed effects variables and obtains OLS estimates for their associated parameters.
    ▸ regularize_covariance_frobenius_norm: Regularizes a covariance matrix based on Frobenius norm and an associated regularization parameter to ensure positive definiteness.
    ▸ regularize_covariance_mean_small_eigenvalues: Regularizes a covariance matrix by using the mean of the smallest eigenvalues, up to a specified percentage of the total number of eigenvalues, to ensure positive definiteness.
    ▸ regularize_covariance_mean_eigenvalues: Regularizes a covariance matrix using the mean of its eigenvalues to ensure positive definiteness.
    ▸ abc_variance_component_estimation: Implements ABC to estimate variance components, critical for genetic analysis.
    ▸ simulate_and_compute_squared_l2_norm: Simulates phenotypic data and evaluates the fit by computing the squared L2 norm distance between observed and simulated values.

## Example

### Estimating phenotypes with WISER

Here's a simple example illustrating the use of the ```estimate_wiser_phenotype``` function to estimate phenotypes :

```R
# ✅ load wiser library and data subsets from the refpop dataset

# load library 
library(wiser)

# load raw individual phenotypic measurements (for subset of 80 genotypes from refpop data)
data("refpop_raw_indiv_pheno_data_subset")
head(refpop_raw_indiv_pheno_data_subset)

# load SNP marker data (for subset of 80 genotypes from refpop data)
data("refpop_geno_data_subset")
head(refpop_geno_data_subset)[,1:10]

# ✅ estimate wiser phenotype for a specified trait,i.e. Flowering_begin here

# define trait
trait_ <- "Flowering_begin"

# increase memory size for future and get optimal whitening method and 
# regularization parameter using k-folds CV
options(future.globals.maxSize = 16 * 1024^3)

opt_white_reg_par <- optimize_whitening_and_regularization(
  omic_df = refpop_geno_data_subset,
  raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
  trait_ = trait_,
  whitening_method_grid = c("ZCA-cor", "PCA-cor", "Cholesky"),
  alpha_grid = c(0.01, 0.1)
)
print(opt_white_reg_par)
# 📌 whitening_method = "ZCA-cor" and alpha_ = 0.01 gives good results generally for 
# estimate_wiser_phenotype(), hence using optimize_whitening_and_regularization() is not 
# always necessary, specially for huge datasets

# apply wiser estimation function using optimized whitening method and regularization parameter
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
  kernel_type = "linear",
  whitening_method = as.character(opt_white_reg_par$opt_whitening_method),
  alpha_ = opt_white_reg_par$opt_alpha_
 )

# ✅ plot and print wiser objects

# plot the density for the estimated phenotypes
dev.new()
plot(density(wiser_obj$wiser_phenotypes$v_hat), main = paste0(trait_, " v_hat"))

# get the estimated fixed effects from the whitening process and OLS, and variance components from ABC
# estimated fixed effects (from whitening and OLS)
print(wiser_obj$wiser_fixed_effect_estimates)

# estimated variance components (from ABC)
print(wiser_obj$wiser_abc_variance_component_estimates)

# ✅ print fundamental property of whitening, i.e. In = W*Σu*W'
id_mat <- wiser_obj$w_mat %*% wiser_obj$sig_mat_u %*% t(wiser_obj$w_mat)
print(id_mat[1:10,1:10])

```

## Authors and References

* Author : Laval Jacquin
* Maintainer : Laval Jacquin jacquin.julien@gmail.com

## References
Kessy, A., Lewin, A., & Strimmer, K. (2018). Optimal whitening and decorrelation. The American Statistician, 72(4), 309-314.

