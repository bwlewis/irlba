# irlba

Implicitly-restarted Lanczos methods for fast truncated singular value
decomposition of sparse and dense matrices (also referred to as partial SVD).
IRLBA stands for Augmented, <b>I</b>mplicitly <b>R</b>estarted <b>L</b>anczos
<b>B</b>idiagonalization <b>A</b>lgorithm. The package provides the following
functions (see help on each for details and examples).

* `irlba()` partial SVD function
* `ssvd()` l1-penalized matrix decomposition for sparse PCA (based on Shen and Huang's algorithm)
* `prcomp_irlba()`  principal components function similar to the `prcomp` function in stats package for computing the first few principal components of large matrices
* `svdr()` alternate partial SVD function based on randomized SVD (see also the [rsvd](https://cran.r-project.org/package=rsvd) package by N. Benjamin Erichson for an alternative implementation)
* `partial_eigen()` a very limited partial eigenvalue decomposition for symmetric matrices (see the [RSpectra](https://cran.r-project.org/package=RSpectra) package for more comprehensive truncated eigenvalue decomposition)

Help documentation for each function includes extensive documentation and
examples. Also see the package vignette, `vignette("irlba", package="irlba")`.

An overview web page is here: https://bwlewis.github.io/irlba/.

## New in 2.4.0

- Re-factored and minimized C code, limiting its use to fast dense matrix multiplication. This change simplifies the package, largely preserves performance for the dense case, and expands Matrix/sparse support to more Matrix classes. It also eliminates direct use of internal hard-to-support SuiteSparse methods. The `fastpath` option is deprecated.
- Finally added support for efficient model deflation (see examples). That means that if your partial SVD subspace isn't big enough, you can efficiently carry on the algorithm from where you left off. This works for smallest and largest portions of the subspace. User-facing behavior is unchanged, but runs more efficiently now.
- Test coverage is expanded.
- Default tolerance slightly reduced and work dimension slightly increased for better default accuracy with minor performance cost.
- Despite the significant internal changes, everything should just work.

## New in 2.3.3

- Several important reported problems with sparse matrices other than class
  dgCMatrix and prcomp_irlba, and many bug fixes contributed by Aaron Lun, see for
  instance: https://github.com/bwlewis/irlba/issues/47.

## New in 2.3.2

- Fixed a regression in `prcomp_irlba()` discovered by Xiaojie Qiu, see https://github.com/bwlewis/irlba/issues/25, and other related problems reported in https://github.com/bwlewis/irlba/issues/32.
- Added rchk testing to pre-CRAN submission tests.
- Fixed a sign bug in `ssvd()` found by Alex Poliakov.

## New in Version 2.3.1

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

## References

* Baglama, James, and Lothar Reichel. "Augmented implicitly restarted Lanczos bidiagonalization methods." SIAM Journal on Scientific Computing 27.1 (2005): 19-42.
* Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp. "Finding structure with randomness: Stochastic algorithms for constructing approximate matrix decompositions." (2009).
* Shen, Haipeng, and Jianhua Z. Huang. "Sparse principal component analysis via regularized low rank matrix approximation." Journal of multivariate analysis 99.6 (2008): 1015-1034.
* Witten, Daniela M., Robert Tibshirani, and Trevor Hastie. "A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis." Biostatistics 10.3 (2009): 515-534.
