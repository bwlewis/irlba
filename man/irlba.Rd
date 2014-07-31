\name{irlba}
\alias{irlba}
\title{Find a few approximate singular values and singular vectors of a matrix. }
\description{
The augmented implicitly restarted Lanczos bi-diagonalization (IRLBA) algorithm
finds a few approximate singular values and corresponding singular vectors of a
matrix using a method of Baglama and Reichel.  It is a fast and
memory-efficient way to compute a partial SVD.
}
\usage{
irlba(A, nv = 5, nu = nv, aug = c("ritz","harm"),
      sigma = c("ls","ss"), maxit = 1000, m_b = nu+5, reorth = TRUE,
      tol = 1e-3, V = NULL, d, U, right_only=FALSE, verbose=FALSE,
      dU, ds, dV)
}
\arguments{
  \item{A}{A double-precision real or complex matrix or real sparse matrix}
  \item{nv}{Number of desired right singular vectors}
  \item{nu}{Number of desired left singular vectors}
  \item{adjust}{Number of extra approximate singular values to compute to enhance convergence}
  \item{aug}{"ritz" for Ritz "harm" for harmonic Ritz vector augmentation}
  \item{sigma}{"ls" for largest few singular values, "ss" for smallest}
  \item{maxit}{Maximum number of iterations}
  \item{m_b}{Dimension of the working subspace}
  \item{reorth}{Logical value full (TRUE) or one-sided (FALSE) reorthogonalization}
  \item{tol}{Convergence is determined when ||A*V - U*S|| <= tol*||A||, where ||A|| is approximated by the largest singular value of all projection matrices. }
  \item{V}{Optional starting vector or right restart subspace}
  \item{d}{Optional restart singular values}
  \item{U}{Optional left restart subspace}
  \item{right_only}{TRUE means only return the right singular vectors (V)}
  \item{verbose}{When TRUE, emit progress messages to stdout}
  \item{dU}{Optional deflation vector (see details)}
  \item{ds}{Optional deflation scalar (see details)}
  \item{dV}{Optional deflation vector (see details)}
}

\value{
Returns a list with entries:
  \item{d}{max (nu, nv) approximate singular values}
  \item{u}{nu approximate left singular vectors (only when right_only=FALSE)}
  \item{v}{nv approximate right singular vectors}
  \item{iter}{The number of Lanczos iterations carried out}
  \item{mprod}{The total number of matrix vector products carried out}
}
\details{
The syntax of \code{irlba} partially follows \code{svd}, with an important
exception. The usual R \code{svd} function always returns a complete set of
singular values, even if the number of singular vectors \code{nu} or \code{nv}
is set less than the maximum. The \code{irlba} function returns a number of
estimated singular values equal to the maximum of the number of specified
singular vectors \code{nu} and \code{nv}.

Use the optional deflation parameters to compute the truncated SVD of
the deflated matrix \eqn{A - dU ds\cdot dV^t}, for example, to efficiently
compute principal components without explicitly forming the centered
matrix. This is important to preserve sparsity in the matrix.

Document restarting XXX


}
\references{"Augmented Implicitly Restarted Lanczos Bidiagonalization Methods",
     J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005.  }
\author{ Adapted for R by B. W. Lewis <blewis@illposed.net>}
\examples{
set.seed(1)

A <- matrix(runif(200),nrow=20)
S <- irlba(A)
S$d

# Compare with svd
svd(A)$d[1:5]

# Principal components
P <- irlba(A, nv=1, tol=1e-5, dU=rep(1,nrow(A)), ds=1, dV=colMeans(A))

# Compare with singular vectors from
Q <- svd(scale(A, center=TRUE, scale=FALSE), nv=1, nu=1)

# The vectors are only guaranteed to be identical up to sign:
cbind(P$v, Q$v)

}
\keyword{matrix}
\keyword{svd}
