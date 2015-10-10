# irlba

Implicitly-restarted Lanczos methods for fast truncated singular value and
symmetric eigenvalue decompositions of sparse and dense matrices. IRLBA stands
for <b>I</b>mplicitly <b>R</b>estarted <b>L</b>anczos <b>B</b>idiagonalization
<b>A</b>lgorithms.

Version 2.0.0 of the package is a major update that includes many changes, new
features, and removal of some old features that did not work well. In
particular, use of harmonic Ritz vector augmentation and the ability to
estimate the smallest singular values of a matrix WAS REMOVED. That method
didn't work very well, and sometimes not at all. It suffered from poor
performance and rarely converged to a solution for ill-conditioned matrices.
Replacements are under consideration but did not make it in to the update.


### Package features that are mostly the same
- Fast truncated singular value decomposition

### Package features that were removed
- Support for estimating smallest singular values
- Harmonic Ritz vector augmentation

### New features
- Support for fast symmetric partial eigenvalue decompositions of real-valued matrices
- Efficient principal components decomposition
- Restarting (picking up where you left off to add more vectors/values)
- Nice syntax for centering and scaling
- Efficient generic subspace deflation (used by PCA)
- Fixed support for complex-valued SVD of dense matrices

## What's still missing

The following new feature didn't make it in but will appear in version 2.1.0:
an implementation of the IRBLB method that trickles out singular triplets in
bounded memory using fast Leja-based accelerating polynomials.

## References

* Augmented Implicitly Restarted Lanczos Bidiagonalization Methods, J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005. (http://www.math.uri.edu/~jbaglama/papers/paper14.pdf)


## Status
<a href="https://travis-ci.org/bwlewis/irlba">
<img src="https://travis-ci.org/bwlewis/irlba.svg?branch=master" alt="Travis CI status"></img>
</a>
[![codecov.io](https://codecov.io/github/bwlewis/irlba/coverage.svg?branch=master)](https://codecov.io/github/bwlewis/irlba?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/irlba)](http://cran.rstudio.com/web/packages/irlba/index.html)
![](http://cranlogs-dev.r-pkg.org/badges/irlba)
