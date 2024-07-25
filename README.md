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

The ```wiser``` package offers *user-friendly* advanced tools for more accurate phenotypic estimation, by leveraging genetic covariance information in the estimation of fixed effects. By employing a whitening transformation followed by successive least squares estimation, ```wiser``` refines fixed effect estimates and eliminates biases, leading to improved phenotypic estimation. This approach is particularly beneficial in complex experimental designs where genetic and environmental factors are intricately linked. ```wiser``` includes methods for computing whitening matrices, fixed effects, residuals, and estimating phenotypes. Additionally, the package offers robust variance component estimation using approximate bayesian computation (ABC), ensuring stability and reliability of estimations. 

The phenotypic estimation performed by `wiser` for a vector $v = (v_1, v_2, \ldots, v_q)'$ of $q$ phenotypes, which approximates a vector $u = (u_1, u_2, \ldots, u_q)'$ of genetic values for $q$ genotypes, is carried out as follows:

$$
\hat{v} = \underset{v \in \mathbb{R}^q}{argmin} || \hat{\xi} - Zv||^2_2 = (Z'Z)^{-1}Z'\hat{\xi} \ \ \ \ (1)
$$ 

where :

* $Z$ is the incidence matrix that links the $q$ phenotypes, which approximate the genetic values of $q$ genotypes, to the individual phenotypic measurements repeated for each genotype in the experimental design.

* $\hat{\xi}$ is the vector of estimated residuals obtained after removing the fixed effects, adjusted for the genetic covariance between individuals in the experimental design.

In this setting, the vector $\hat{\xi}$ is estimated as follows :

$$
\hat{\xi} = Y - \tilde{X}\hat{\beta} = WX\hat{\beta} \ \ \ \ (2)
$$ 

where :

* $Y = (Y_1,Y_2, \dots, Y_n)'$ is the vector of $n$ individual phenotypic measurements in the experimental design, with repeated measurements for each genotype.

* $\tilde{X} (= WX)$ is the transformed (whitened) design matrix that links the estimated fixed effects, adjusted for genetic covariance, to the individual phenotypic measurements.

* $\hat{\beta} = (\hat{\beta}_1, \hat{\beta}_2, \ldots ,\hat{\beta}_l)' = (\tilde{X}'\tilde{X})^{-1}\tilde{X}'Y$ is the vector of $l$ estimated fixed effects, adjusted for genetic covariance between individuals in the experimental design.

* $W$ is an $n$ x $n$ whitening matrix. By default, $W = L^{-1}$ in ```wiser```, where $L$ is a full-rank lower triangular matrix obtained via the Cholesky decomposition of the genetic covariance matrix $\Sigma_u = \sigma^2_uZKZ'$ of the vector $U = Zu = (U_1,U_2, \ldots, U_n)$, which corresponds to the vector of individual genetic effects in the experimental design. The Cholesky decomposition is given by $\Sigma_u = LL'$. It is assumed that $u \sim \mathcal{N}_q(0,\sigma^2_uK)$, where $K$ is the genetic covariance matrix between the genotypes estimated using SNP marker data.

* $X$ is the original design matrix that links fixed effects to the individual phenotypic measurements in the experimental design.

In ```wiser```, the matrix $K$ is estimated using SNP markers and is also known as the genomic covariance matrix or the Gram matrix. Three kernel functions are provided to build $K$: ```linear```, ```gaussian```, and ```identity```. Among these, the `identity` kernel is not recommended due to its poor performance for genomic prediction. In contrast, the `linear` kernel has shown the best practical performance in this context. The `identity` kernel is primarily used for benchmarking purposes, as detailed in Jacquin $\textit{et al.}$ (2024). Specifically, the `identity` kernel facilitates the assessment of the impact of including versus excluding the covariance structure, as well as evaluating the effects of different methods for computing the covariance structure.

#### Remarks :

The rationale for transforming $X$ into $\tilde{X}$ through whitening, to correct for genetic covariance among individuals in the estimation of fixed effects, is comprehensively addressed in Jacquin $\textit{et al.}$ (2024). In contrast, the chosen approach of successive ordinary least squares (OLS) in ```wiser``` avoids making assumptions about the properties of $\beta$ and $v$.

Specifically, ```wiser``` does not assume that $v$ is a random vector drawn from a distribution with a specified covariance matrix. This approach prevents enforcing an unnecessary covariance structure during the estimation of $v$, which could be detrimental. For example, assuming $v \sim \mathcal{N}_q(0,\sigma^2_vI_q)$ is often unrealistic and would lead to using a decorrelated covariance structure in the best linear unbiased predictor (BLUP) of $v$. Furthermore, and crucially, the successive OLS procedure ensures that the estimated residuals or $v$ are not linked to the SNP marker data as shown in Jacquin $\textit{et al.}$ (2024).

## Installation

You can install the latest version of the ```wiser``` package with:

```R
install.packages("devtools")
library(devtools)
install_github("ljacquin/wiser")
```

## Key Features


## Main Functions


## Examples

### Estimating phenotypes with WISER

Here's a simple example illustrating the use of  

```R
```

## Authors and References

* Author: Laval Jacquin
* Maintainer: Laval Jacquin jacquin.julien@gmail.com

## References


