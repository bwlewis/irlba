# irlba

Implicitly-restarted Lanczos methods for fast truncated singular value decomposition
of sparse and dense matrices (also referred to as partial SVD).  IRLBA stands
for Augmented, <b>I</b>mplicitly <b>R</b>estarted <b>L</b>anczos
<b>B</b>idiagonalization <b>A</b>lgorithm. The package provides the following
functions (see help on each for details and examples).

* `irlba` main partial SVD function
* `prcomp_irlba`  principal components function similar to the `prcomp` function in stats package for computing the first few principal components of large matrices
* `partial_eigen` a very limited partial eigenvalue decomposition for symmetric matrices (see the [RSpectra](https://cran.r-project.org/package=RSpectra) package for more comprehensive truncated eigenvalue decomposition)

Help documentation for each function includes examples. Also see the package
vignette, `vignette("irlba", package="irlba")`, and demo,
`demo("custom_matrix_multiply", package="irlba")`.

## What's new?

Version 2.2.0 includes stronger convergence criteria and a new argument `svtol`
associated with that. The new approach helps guarantee more accurate solutions
for some difficult problems. The tradeoff is that the default behavior is a
little slower than before because at least two Lanczos iterations are always
run. The new convergence behavior can be disabled with `svtol=Inf`.

## Deprecated features

The `mult` argument is deprecated and will be removed in a future version. We
now recommend simply defining a custom class with a custom multiplcation
operator.  The example below illustrates the old and new approaches, as well as
a few other equivalent approaches.

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

## Wishlist

- Optional block implementation for some use cases
- Re-introduce a solver for smallest singular values

## References

* Augmented Implicitly Restarted Lanczos Bidiagonalization Methods, J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005. (http://www.math.uri.edu/~jbaglama/papers/paper14.pdf)


## Status
<a href="https://travis-ci.org/bwlewis/irlba">
<img src="https://travis-ci.org/bwlewis/irlba.svg?branch=master" alt="Travis CI status"></img>
</a>
[![codecov.io](https://codecov.io/github/bwlewis/irlba/coverage.svg?branch=master)](https://codecov.io/github/bwlewis/irlba?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/irlba)](https://cran.r-project.org/package=irlba)
![](http://cranlogs.r-pkg.org/badges/irlba)
