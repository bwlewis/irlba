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

Version 2.1.2 of the package includes a convenience `prcomp`-like function for
computing principal components and numerous bug fixes and numerical stability
improvements in edge cases.


## References

* Augmented Implicitly Restarted Lanczos Bidiagonalization Methods, J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005. (http://www.math.uri.edu/~jbaglama/papers/paper14.pdf)


## Status
<a href="https://travis-ci.org/bwlewis/irlba">
<img src="https://travis-ci.org/bwlewis/irlba.svg?branch=master" alt="Travis CI status"></img>
</a>
[![codecov.io](https://codecov.io/github/bwlewis/irlba/coverage.svg?branch=master)](https://codecov.io/github/bwlewis/irlba?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/irlba)](https://cran.r-project.org/package=irlba)
![](http://cranlogs.r-pkg.org/badges/irlba)
