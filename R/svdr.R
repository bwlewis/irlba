#' Find a few approximate largest singular values and corresponding
#' singular vectors of a matrix.
#'
#' The randomized method for truncated SVD by P. G. Martinsson and colleagues
#' finds a few approximate largest singular values and corresponding
#' singular vectors of a sparse or dense matrix. It is a fast and
#' memory-efficient way to compute a partial SVD, similar in performance
#' for many problems to \code{\link{irlba}}. The \code{svdr} method
#' is a block method and may produce more accurate estimations with
#' less work for problems with clustered large singular values (see
#' the examples). In other problems, \code{irlba} may exhibit faster
#' convergence.
#'
#' Also see an alternate implementation (\code{rsvd}) of this method by N. Benjamin Erichson
#' in the https://cran.r-project.org/package=rsvd package.
#'
#' @param x numeric real- or complex-valued matrix or real-valued sparse matrix.
#' @param k dimension of subspace to estimate (number of approximate singular values to compute).
#' @param tol stop iteration when the largest absolute relative change in estimated singular
#'   values from one iteration to the next falls below this value.
#' @param it maximum number of algorithm iterations.
#' @param extra number of extra vectors of dimension \code{ncol(x)}, larger values generally improve accuracy (with increased
#' computational cost).
#' @param center optional column centering vector whose values are implicitly subtracted from each
#'   column of \code{A} without explicitly forming the centered matrix (preserving sparsity).
#'   Optionally specify \code{center=TRUE} as shorthand for \code{center=colMeans(x)}.
#'   Use for efficient principal components computation.
#' @param Q optional initial random matrix, defaults to a matrix of size \code{ncol(x)} by \code{k + extra} with
#' entries sampled from a normal random distribution.
#' @param return.Q if \code{TRUE} return the \code{Q} matrix for restarting (see examples).
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{d:}{ k approximate singular values}
#'   \item{u:}{ k approximate left singular vectors}
#'   \item{v:}{ k approximate right singular vectors}
#'   \item{mprod:}{ total number of matrix products carried out}
#'   \item{Q:}{ optional subspace matrix (when \code{return.Q=TRUE})}
#' }
#' @seealso \code{\link{irlba}}, \code{\link{svd}}, \code{rsvd} in the rsvd package
#' @references
#' Finding structure with randomness: Stochastic algorithms for constructing
#' approximate matrix decompositions N. Halko, P. G. Martinsson, J. Tropp. Sep. 2009.
#' @examples
#' set.seed(1)
#'
#' A <- matrix(runif(400), nrow=20)
#' svdr(A, 3)$d
#'
#' # Compare with svd
#' svd(A)$d[1:3]
#'
#' # Compare with irlba
#' irlba(A, 3)$d
#'
#' \dontrun{
#' # A problem with clustered large singular values where svdr out-performs irlba.
#' tprolate <- function(n, w=0.25)
#' {
#'   a <- rep(0, n)
#'   a[1] <- 2 * w
#'   a[2:n] <- sin( 2 * pi * w * (1:(n-1)) ) / ( pi * (1:(n-1)) )
#'   toeplitz(a)
#' }
#'
#' x <- tprolate(512)
#' set.seed(1)
#' tL <- system.time(L <- irlba(x, 20))
#' tR <- system.time(R <- svdr(x, 20))
#' S <- svd(x)
#' plot(S$d)
#' data.frame(time=c(tL[3], tR[3]),
#'            error=sqrt(c(crossprod(L$d - S$d[1:20]), crossprod(R$d - S$d[1:20]))),
#'            row.names=c("IRLBA", "Randomized SVD"))
#'
#' # But, here is a similar problem with clustered singular values where svdr
#' # doesn't out-perform irlba as easily...clusters of singular values are,
#' # in general, very hard to deal with!
#' # (This example based on https://github.com/bwlewis/irlba/issues/16.)
#' set.seed(1)
#' s <- svd(matrix(rnorm(200 * 200), 200))
#' x <- s$u %*% (c(exp(-(1:100)^0.3) * 1e-12 + 1, rep(0.5, 100)) * t(s$v))
#' tL <- system.time(L <- irlba(x, 5))
#' tR <- system.time(R <- svdr(x, 5))
#' S <- svd(x)
#' plot(S$d)
#' data.frame(time=c(tL[3], tR[3]),
#'            error=sqrt(c(crossprod(L$d - S$d[1:5]), crossprod(R$d - S$d[1:5]))),
#'            row.names=c("IRLBA", "Randomized SVD"))
#' }
#' @export
svdr <- function(x, k, tol=1e-5, it=100L, extra=min(10L, dim(x) - k), center=NULL, Q=NULL, return.Q=FALSE)
{
  eps2 <- .Machine$double.eps ^ (4 / 5)
  n <- min(ncol(x), k + extra)
  if (isTRUE(center)) center <- colMeans(x)
  if (is.null(Q)) Q <-  matrix(rnorm(ncol(x) * n), ncol(x))
  d <- rep(0, k)
  for (j in 1:it)
  {
    if (is.null(center))
    {
      Q <- qr.Q(qr(x %*% Q))
      B <- crossprod(x, Q)
      Q <- qr.Q(qr(B))
    } else
    {
      Q <- qr.Q(qr(x %*% Q - rep(1, nrow(x)) %*% crossprod(center, Q)))
      B <- crossprod(Q, x) - tcrossprod(crossprod(Q, rep(1, nrow(x))), center)
      Q <- qr.Q(qr(t(B)))
    }
    d1 <- svd(B, nu=0, nv=0)$d[1:k]
    idx <- d1 > eps2
    if (all(! idx)) break
    if (max(abs((d1[idx] - d[idx]) / d[idx])) < tol) break
    d <- d1
  }
  if (return.Q) Q1 <- Q
  if (is.null(center))
  {
    Q <- qr.Q(qr(x %*% Q))
    B <- crossprod(Q, x)
  } else
  {
    Q <- qr.Q(qr(x %*% Q - rep(1, nrow(x)) %*% crossprod(center, Q)))
    B <- crossprod(Q, x) - tcrossprod(crossprod(Q, rep(1, nrow(x))), center)
  }
  s <- svd(B)
  s$u <- Q %*% s$u
  s$u <- s$u[, 1:k]
  s$d <- s$d[1:k]
  s$v <- s$v[, 1:k]
  s$mprod <- 2 * j + 1
  if (return.Q) s$Q <- Q1
  s
}
