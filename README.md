# irlba

Implicitly-restarted Lanczos methods for fast truncated singular value and
symmetric eigenvalue decompositions of sparse and dense matrices.  IRLBA stands
for Augmented, <b>I</b>mplicitly <b>R</b>estarted <b>L</b>anczos
<b>B</b>idiagonalization <b>A</b>lgorithm.

Version 2.1.0 of the package includes a convenience `prcomp`-like function for
computing principal components and a fast C-language implementation for
improved computational speed. The original R algorithm implementation is
maintained for some special cases and for refernce.

## TODO

Implement Leja-point polynomial acceleration for even more efficient memory use
for extremely large problems.

## References

* Augmented Implicitly Restarted Lanczos Bidiagonalization Methods, J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005. (http://www.math.uri.edu/~jbaglama/papers/paper14.pdf)


## Status
<a href="https://travis-ci.org/bwlewis/irlba">
<img src="https://travis-ci.org/bwlewis/irlba.svg?branch=master" alt="Travis CI status"></img>
</a>
[![codecov.io](https://codecov.io/github/bwlewis/irlba/coverage.svg?branch=master)](https://codecov.io/github/bwlewis/irlba?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/irlba)](https://cran.r-project.org/package=irlba)
![](http://cranlogs.r-pkg.org/badges/irlba)
