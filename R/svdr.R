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
#' @param x numeric real- or complex-valued matrix or real-valued sparse matrix.
#' @param k dimension of subspace to estimate (number of approximate singular values to compute).
#' @param it number of algorithm iterations
#' @param extra number of extra vectors of dimension \code{ncol(x)}, large values generally improve accuracy and performance.
#' @param center optional column centering vector whose values are subtracted from each
#'   column of \code{A} or, optionally, use \code{center=TRUE} as shorthand for \code{center=colMeans(x)}.
#'   Used for efficient principal components computation.
#' @param Q optiona initial random matrix, defaults to a matrix of size \code{ncol(x)} by \code{k + extra} with
#' entries sampled from a normal random distribution.
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{d:}{ k approximate singular values}
#'   \item{u:}{ k approximate left singular vectors}
#'   \item{v:}{ k approximate right singular vectors}
#'   \item{mprod:}{ The total number of matrix vector products carried out}
#' }
#' @seealso \code{\link{irlba}}, \code{\link{svd}}
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
svdr <- function(x, k, it=3, extra=10, center=NULL, Q=NULL)
{
  n <- min(ncol(x), k + extra)
  if (isTRUE(center)) center <- colMeans(x)
  if (is.null(Q)) Q <-  matrix(rnorm(ncol(x) * n), ncol(x))
  for (j in 1:it)
  {
    if (is.null(center))
    {
      Q <- qr.Q(qr(x %*% Q))
      Q <- qr.Q(qr(t(crossprod(Q, x))))   # a.k.a. Q <- qr.Q(qr(t(x) %*% Q)), but avoiding t(x)
    } else
    {
      Q <- qr.Q(qr(x %*% Q - rep(1, nrow(x)) %*% crossprod(center, Q)))
      Q <- qr.Q(qr(t(crossprod(Q, x) - tcrossprod(crossprod(Q, rep(1, nrow(x))), center))))
    }
  }
  if (is.null(center))
  {
    Q <- qr.Q(qr(x %*% Q))
    B <- t(Q) %*% x
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
  s$mprod <- 2 * it + 1
  s
}
