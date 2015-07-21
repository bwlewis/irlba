# IRL

Augmented, implicitly restarted Lanczos methods for fast truncated singular
value and symmetric eigenvalue decompositions of sparse and dense real-valued
matrices.

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
- Support for fast symmetric partial eigenvalue decompositions
- Efficient principal components decomposition
- Restarting
- Efficient subspace deflation (used by PCA)
- Nice syntax for centering and scaling

## What's still missing

The following new feature didn't make it in but will appear in version 2.1.0:
an implementation of the IRBLB method that trickles out singular triplets in
bounded memory using fast Leja-based accelerating polynomials.

## References

* Jim Baglama and Lothar Reichel, Augmented Implicitly Restarted Lanczos Bidiagonalization Methods, http://www.math.uri.edu/~jbaglama/papers/paper14.pdf
