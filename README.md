# irlba

Implicitly-restarted Lanczos methods for fast truncated singular value decomposition
of sparse and dense matrices (also referred to as partial SVD).  IRLBA stands
for Augmented, <b>I</b>mplicitly <b>R</b>estarted <b>L</b>anczos
<b>B</b>idiagonalization <b>A</b>lgorithm. The package provides the following
functions (see help on each for details and examples).

* `irlba()` partial SVD function
* `svdr()` alternate partial SVD function based on randomized SVD
* `svds()` l1-penalized matrix decompoisition for sparse PCA (based on Shen and Huang's algorithm)
* `prcomp_irlba()`  principal components function similar to the `prcomp` function in stats package for computing the first few principal components of large matrices
* `partial_eigen()` a very limited partial eigenvalue decomposition for symmetric matrices (see the [RSpectra](https://cran.r-project.org/package=RSpectra) package for more comprehensive truncated eigenvalue decomposition)

Help documentation for each function includes extensive documentation and
examples. Also see the package vignette, `vignette("irlba", package="irlba")`.

An overview web page is here: https://bwlewis.github.io/irlba/.

## What's new in Version 2.3.1?

- Fixed an `irlba()` bug associated with centering (PCA), see https://github.com/bwlewis/irlba/issues/21.
- Fixed `irlba()` scaling to conform to `scale`, see https://github.com/bwlewis/irlba/issues/22.
- Improved `prcomp_irlba()` from a suggestion by N. Benjamin Erichson, see https://github.com/bwlewis/irlba/issues/23.
- Significanty changed/improved `svdr()` convergence criterion.
- Added a version of Shen and Huang's Sparse PCA/SVD L1-penalized matrix decomposition (`ssvd()`).
- Fixed valgrind errors.


## Deprecated features

I will remove `partial_eigen()` in a future version. As its documentation
states, users are better off using the RSpectra package for eigenvalue
computations (although not generally for singular value computations).

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
## [1] 1.820227 1.622988 1.067185

# Compare with:
irlba(A, 3, scale=col_scale)$d
## [1] 1.820227 1.622988 1.067185

# Compare with:
svd(sweep(A, 2, col_scale, FUN=`/`))$d[1:3]
## [1] 1.820227 1.622988 1.067185

# ------------------ new way ----------------------------------------------
setClass("scaled_matrix", contains="matrix", slots=c(scale="numeric"))
setMethod("%*%", signature(x="scaled_matrix", y="numeric"), function(x ,y) x@.Data %*% (y / x@scale))
setMethod("%*%", signature(x="numeric", y="scaled_matrix"), function(x ,y) (x %*% y@.Data) / y@scale)
a <- new("scaled_matrix", A, scale=col_scale)

irlba(a, 3)$d
## [1] 1.820227 1.622988 1.067185
```

We have learned that using R's existing S4 system is simpler, easier, and more
flexible than using custom arguments with idiosyncratic syntax and behavior.
We've even used the new approach to implement distributed parallel matrix
products for very large problems with amazingly little code.

## Wishlist / help wanted...

- Optional block implementation for some use cases
- More Matrix classes supported in the fast code path
- Help improving the solver for smallest singular values in tricky cases (basically, for ill-conditioned problems)

## References

* Baglama, James, and Lothar Reichel. "Augmented implicitly restarted Lanczos bidiagonalization methods." SIAM Journal on Scientific Computing 27.1 (2005): 19-42.
* Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp. "Finding structure with randomness: Stochastic algorithms for constructing approximate matrix decompositions." (2009).
* Shen, Haipeng, and Jianhua Z. Huang. "Sparse principal component analysis via regularized low rank matrix approximation." Journal of multivariate analysis 99.6 (2008): 1015-1034.
* Witten, Daniela M., Robert Tibshirani, and Trevor Hastie. "A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis." Biostatistics 10.3 (2009): 515-534.

## Status
<a href="https://travis-ci.org/bwlewis/irlba">
<img src="https://travis-ci.org/bwlewis/irlba.svg?branch=master" alt="Travis CI status"></img>
</a>
<a href="https://codecov.io/gh/bwlewis/irlba">
  <img src="https://codecov.io/gh/bwlewis/irlba/branch/master/graph/badge.svg" alt="Codecov" />
</a>
<a href="https://www.r-pkg.org/pkg/irlba">
  <img src="http://cranlogs.r-pkg.org/badges/irlba" />
</a>
<a href="https://cran.r-project.org/package=irlba">
  <img src="https://www.r-pkg.org/badges/version/irlba" />
</a>
