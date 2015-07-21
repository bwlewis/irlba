# IRL

Augmented, implicitly restarted Lanczos methods for fast truncated singular
value and symmetric eigenvalue decompositions of sparse and dense real-valued
matrices.

Version 2.0.0 of the package is a major update that includes many changes, new
features, and removal of some old features that did not work well. In
particular, use of harmonic Ritz vector augmentation and the ability to
estimate the smallest singular values of a matrix was removed. That method
doesn't work very well, and sometimes not at all. It suffered from poor
performance and rarely converged to a solution for ill-conditioned matrices.
Replacements are under consideration but did not make it in to the update.

### Mostly the same
- Fast truncated singular value decomposition

### Removed
- Support for estimating smallest singular values
- Harmonic Ritz vector augmentation

### New
- Support for fast symmetric eigenvalue decompositions
- Efficient principal components decomposition
- Restarting
- Efficient subspace deflation (used by PCA)
- Limited support for solution of some large-scale network problems
- Fast estimates of large sparse thresholded correlation matrices

## What's still missing

The following new feature didn't make it in but will appear in version 2.1.0:
an implementation of the IRBLB method that trickles out singular triplets in
bounded memory using fast Leja-based accelerating polynomials.

## References

This work is almost totally based on work of Baglama and Reichel, and also Fenu
and Rodriguez.  Here are links to important papers. Read them!

* Augmented Implicitly Restarted Lanczos Bidiagonalization Methods http://www.math.uri.edu/~jbaglama/papers/paper14.pdf
* Network analysis via partial spectral factorization and Gauss quadrature http://www.math.kent.edu/~reichel/publications/netwrk.pdf
* Quadrature Rule-Based Bounds for Functions of Adjacency Matrices http://www.mathcs.emory.edu/~benzi/Web_papers/adjacency_paper.pdf
