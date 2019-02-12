---
title: "The irlba package"
output:
html_document:
    theme: cerulean
---
<style>
body{ font-size: 15pt; }
blockquote{ font-size: 14pt; }
pre{ font-size: 14pt; }
h1{ font-size: 28pt; }
h2,h3,h4,h5,h6{ font-size: 24pt; }
.main-container{ max-width: 85%; }
</style>

The augmented implicitly restarted Lanczos bidiagonalization algorithm (IRLBA)
finds a few approximate largest singular values and corresponding singular
vectors of a sparse or dense matrix using a method of Baglama and Reichel.

> J. Baglama and L. Reichel, SIAM J. Sci. Comput. (2005). (http://www.math.uri.edu/~jbaglama/papers/paper14.pdf)

It is a fast and memory-efficient way to compute a partial SVD, principal
components, and some specialized partial eigenvalue decompositions.  I
introduced Baglama and Reichel's irlba algorithm to the R world in a talk at
the useR! conference in Dortmund way back in 2008.

The package provides the following functions (see help on each for details and
examples).

* `irlba()` partial SVD function
* `ssvd()` l1-penalized matrix decompoisition for sparse PCA (based on Shen and Huang's algorithm)--see https://bwlewis.github.io/irlba/ssvd.html  for more details
* `prcomp_irlba()`  principal components function similar to the `prcomp` function in stats package for computing the first few principal components of large matrices
* `svdr()` alternate partial SVD function based on randomized SVD
* `partial_eigen()` a very limited partial eigenvalue decomposition for symmetric matrices (see the [RSpectra](https://cran.r-project.org/package=RSpectra) package for more comprehensive truncated eigenvalue decomposition); see also https://bwlewis.github.io/irlba/comparison.html for more notes on RSpectra.

## Status

<a href="https://travis-ci.org/bwlewis/irlba">
<img src="https://travis-ci.org/bwlewis/irlba.svg?branch=master" alt="Travis CI status"/>
</a>
<a href="https://codecov.io/gh/bwlewis/irlba">
  <img src="https://codecov.io/gh/bwlewis/irlba/branch/master/graph/badge.svg" alt="Codecov" />
</a>
[![](https://cranlogs.r-pkg.org/badges/irlba)]("https://www.r-pkg.org/pkg/irlba")
<a href="https://cran.r-project.org/package=irlba">
  <img src="https://www.r-pkg.org/badges/version/irlba" />
</a>

## References

For your reading pleasure...

- Baglama, James, and Lothar Reichel. "Augmented implicitly restarted Lanczos bidiagonalization methods." SIAM Journal on Scientific Computing 27.1 (2005): 19-42.
- Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp. "Finding structure with randomness: Stochastic algorithms for constructing approximate matrix decompositions." (2009).
- Shen, Haipeng, and Jianhua Z. Huang. "Sparse principal component analysis via regularized low rank matrix approximation." Journal of multivariate analysis 99.6 (2008): 1015-1034.
- Witten, Daniela M., Robert Tibshirani, and Trevor Hastie. "A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis." Biostatistics 10.3 (2009): 515-534.

And, some interesting related recent work by Musco and colleagues:

- Jin, Chi, et al. "Robust shift-and-invert preconditioning: Faster and more sample efficient algorithms for eigenvector computation." arXiv preprint arXiv:1510.08896 (2015).
- Musco, Cameron, and Christopher Musco. "Randomized block krylov methods for stronger and faster approximate singular value decomposition." Advances in Neural Information Processing Systems. 2015.

## New in 2.3.3 (Februrary, 2019)

This is a bugfix release.

- Man contributions by <a href="https://github.com/LTLA">LTLA</a> fixing issues
in prcomp-style methods, among others.

## New in 2.3.2 (January, 2018)

- Fixed a regression in `prcomp_irlba()` discovered by Xiaojie Qiu, see https://github.com/bwlewis/irlba/issues/25, and other related problems reported in https://github.com/bwlewis/irlba/issues/32.
- Added rchk testing to pre-CRAN submission tests.
- Fixed a sign bug in `ssvd()` found by Alex Poliakov.


# What's new in Version 2.3.1 (October, 2017)?

- Fixed an `irlba()` bug associated with centering (PCA), see https://github.com/bwlewis/irlba/issues/21.
- Fixed `irlba()` scaling to conform to `scale`, see https://github.com/bwlewis/irlba/issues/22.
- Improved `prcomp_irlba()` from a suggestion by N. Benjamin Erichson, see https://github.com/bwlewis/irlba/issues/23.
- Significanty changed/improved `svdr()` convergence criterion.
- Added a version of Shen and Huang's Sparse PCA/SVD L1-penalized matrix decomposition (`ssvd()`).
- Fixed many valgrind errors inadvertently introduced in the C code a while ago.


# What's new in version 2.2.1 (May, 2017)

- IRLBA-based functions include stronger convergence criteria and a new argument called svtol for that. The new approach helps guarantee more accurate solutions for some difficult problems. The tradeoff is that the default behavior is a little slower than before because at least two Lanczos iterations are always run. The new convergence behavior can be disabled with svtol=Inf.
- The new package version includes a new function, svdr()--another state of the art truncated SVD method based on the randomized SVD algorithm of Gunnar Martinsson and others. Both irlba() and svdr() work well. Svdr uses a block method and may exhibit better convergence in problems where the largest singular values are clustered. See the documentation and examples in the package. (Block versions of irlba exists, but are not yet implemented by this R package--something coming in the future.)
- We re-introduced a solver for estimating the smallest singular values of a matrix and associated singular vector spaces. The solver is based on the oringial Harmonic Ritz vector augmentation method of Baglama and Reichel.  Beware that this method is somewhat experimental and may fail to converge, or may converge poorly, to estimated singular values for very ill-conditioned matrices. Along with block methods for irlba, this is an active area of work--feel free to contribute!
- The mult() argument is deprecated and will be removed in a future version. Instead, we now recommend simply defining a custom class with a custom multiplcation operator. The example below illustrates the old and new approaches.

```{r, eval=FALSE}
library(irlba)
set.seed(1)
A <- matrix(rnorm(100), 10)
# Define a custom matrix multiplication function that scales the columns of A
# to have unit norm (cf the scale option).

# ------------------ the new way ------------------------------------------
# Simply using an S4 class
setClass("scaled_matrix", contains="matrix", slots=c(scale="numeric"))
setMethod("%*%", signature(x="scaled_matrix", y="numeric"), function(x ,y) x@.Data %*% (y / x@scale))
setMethod("%*%", signature(x="numeric", y="scaled_matrix"), function(x ,y) (x %*% y@.Data) / y@scale)
a <- new("scaled_matrix", A, scale=col_scale)

irlba(a, 3)$d
## [1] 1.820227 1.622988 1.067185


# ------------------ the old way ------------------------------------------
col_scale <- sqrt(apply(A, 2, crossprod))
mult <- function(x, y)
        {
          # check if x is a  vector
          if (is.vector(x))
          {
            return((x %*% y) / col_scale)
          }
          # else x is the matrix
          x %*% (y / col_scale)
        }
irlba(A, 3, mult=mult)$d
## [1] 1.820227 1.622988 1.067185


# ------------------ compare ----------------------------------------------
irlba(A, 3, scale=col_scale)$d
## [1] 1.820227 1.622988 1.067185

# Or,
svd(sweep(A, 2, col_scale, FUN=`/`))$d[1:3]
## [1] 1.820227 1.622988 1.067185
```


# Vignettes

The package vignette (PDF): https://cran.r-project.org/web/packages/irlba/irlba.pdf.

This example uses a special one-sided basis option and a custom matrix product
to compute the top three principal components of the entire 1000 Genomes
Project variant data set (whole genome variants for 2,504 people):
http://bwlewis.github.io/1000_genomes_examples/PCA_whole_genome.html The
example optionally works in parallel and finishes in only 8 minutes on a 4 node
Linux cluster.

https://bwlewis.github.io/irlba/comparison.html compares the irlba package
with the RSpectra package, a high-quality eigenvalue solver for R.

## Other applications

The newe `tcor` algorithm (http://arxiv.org/abs/1512.07246) and R package
(https://github.com/bwlewis/tcor -- not yet on CRAN) use irlba to compute
thresholded correlation matrices very quickly. This example computes the most
highly correlated gene expressions from the Cancer Genome Atlas RNASeq gene
expression data for breast cancer:
https://github.com/bwlewis/tcor/blob/master/vignettes/brca.Rmd.

Many more applications will appear here soon, check back!




# Deprecated features

The `mult` argument is deprecated and will be removed in a future version. We
now recommend simply defining a custom class with a custom multiplcation
operator.  The example below illustrates the old and new approaches.

```{r}
library(irlba)
set.seed(1)
A <- matrix(rnorm(100), 10)

# ------------------ old way ----------------------------------------------
# A custom matrix multiplication function that scales the columns of A
# (cf the scale option). This function scales the columns of A to unit norm.
col_scale <- sqrt(apply(A, 2, crossprod))
mult <- function(x, y)
        {
          # check if x is a  vector
          if (is.vector(x))
          {
            return((x %*% y) / col_scale)
          }
          # else x is the matrix
          x %*% (y / col_scale)
        }
irlba(A, 3, mult=mult)$d

# Compare with:
irlba(A, 3, scale=col_scale)$d

# Compare with:
svd(sweep(A, 2, col_scale, FUN=`/`))$d[1:3]

# ------------------ new way ----------------------------------------------
setClass("scaled_matrix", contains="matrix", slots=c(scale="numeric"))
setMethod("%*%", signature(x="scaled_matrix", y="numeric"), function(x ,y) x@.Data %*% (y / x@scale))
setMethod("%*%", signature(x="numeric", y="scaled_matrix"), function(x ,y) (x %*% y@.Data) / y@scale)
a <- new("scaled_matrix", A, scale=col_scale)

irlba(a, 3)$d
```

We have learned that using R's existing S4 system is simpler, easier, and more
flexible than using custom arguments with idiosyncratic syntax and behavior.
We've even used the new approach to implement distributed parallel matrix
products for very large problems with amazingly little code.
