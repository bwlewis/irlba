IRL
===

Implicitly restarted Lanczos methods for R

## NOTE

The package name is "irlba" for backward compatibility with the old irlba
package on CRAN.


## WARNING

This is an "alpha" version of the new irlba package, which remains under
substantial development. Not all planned features are implemented or even
working. This package will replace the "irlba" package on CRAN sometime
over the autumn of 2014.

Check back here periodically for updates.

## Plan for the package

The package contains an updated irlba algorithm for truncated SVD that supports
restarting and other new features.

It will soon contain:

* "Top M" style network analysis methods for centrality measures of
directed and undirected networks.
* Golub-Meurant Gauss/Gauss-Radau quadrature bounds for ordering the
"top M" results, based on the work of Benzi, Boito and Estrada.
* An implementation of the IRBLB method that trickles out singular
triplets in bounded memory using fast Leja-based accelerating polynomials.
* Some kind of symmetric implementation, either IRBLeigs or adaptations
of IRBLB/IRLBA to symmetric problems.

## References

This work is almost totally based on work of Baglama, Reichel, Fenu, Rodriguez,
and Benzi, Boito and Estrada. Here are links to important papers. Read them!

* Network analysis via partial spectral factorization and Gauss quadrature http://www.math.kent.edu/~reichel/publications/netwrk.pdf
* Quadrature Rule-Based Bounds for Functions of Adjacency Matrices http://www.mathcs.emory.edu/~benzi/Web_papers/adjacency_paper.pdf
* Augmented Implicitly Restarted Lanczos Bidiagonalization Methods http://www.math.uri.edu/~jbaglama/papers/paper14.pdf
