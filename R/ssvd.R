#' Sparse regularized low-rank matrix approximation.
#' 
#' Estimate an \eqn{{\ell}1}{l1}-penalized
#' singular value or principal components decomposition (SVD or PCA) that introduces sparsity in the
#' regularized right singular vectors (the PCA loading vectors), based on the fast and memory-efficient
#' sPCA-rSVD algorithm of Haipeng Shen and Jianhua Huang.
#' @param x A numeric real- or complex-valued matrix or real-valued sparse matrix.
#' @param k Matrix rank of the computed decomposition (see the Details section below).
#' @param n Number of nonzero components in the right singular vectors. If \code{k > 1},
#'        then a single value of \code{n} specifies the number of nonzero components
#'        in each regularized right singular vector. Or, specify a vector of length
#'        \code{k} indicating the number of desired nonzero components in each
#'        returned vector. See the examples.
#' @param maxit Maximum number of soft-thresholding iterations.
#' @param tol Convergence is determined when \eqn{\|U_j - U_{j-1}\|_F < tol}{||U_j - U_{j-1}||_F < tol}, where \eqn{U_j} is the matrix of estimated left regularized singular vectors at iteration \eqn{j}.
#' @param center a logical value indicating whether the variables should be
#'        shifted to be zero centered. Alternately, a centering vector of length
#'        equal the number of columns of \code{x} can be supplied. Use \code{center=TRUE}
#'        to perform a regularized sparse PCA.
#' @param scale. a logical value indicating whether the variables should be
#'        scaled to have unit variance before the analysis takes place.
#'        Alternatively, a vector of length equal the number of columns of \code{x} can be supplied.
#'
#'        The value of \code{scale} determines how column scaling is performed
#'        (after centering).  If \code{scale} is a numeric vector with length
#'        equal to the number of columns of \code{x}, then each column of \code{x} is
#'        divided by the corresponding value from \code{scale}.  If \code{scale} is
#'        \code{TRUE} then scaling is done by dividing the (centered) columns of
#'        \code{x} by their standard deviations if \code{center=TRUE}, and the
#'        root mean square otherwise.  If \code{scale} is \code{FALSE}, no scaling is done.
#'        See \code{\link{scale}} for more details.
#' @param alpha Optional  scalar regularization parameter between zero and one (see Details below).
#' @param  tsvd Optional initial rank-k truncated SVD or PCA (skips computation if supplied).
#' @param ... Additional arguments passed to \code{\link{irlba}}.
#' @details
#' The \code{ssvd} function implements a version of an algorithm by
#' Shen and Huang that computes a penalized SVD or PCA that introduces
#' sparsity in the loading vectors by solving a penalized least squares problem.
#' The algorithm in the rank 1 case finds vectors \eqn{u, w}{u, w} that minimize
#' \deqn{\|x - u w^T\|_F^2 + \lambda \|w\|_1}{||x - u w^T||_F^2 + lambda||w||_1}
#' such that \eqn{\|u\| = 1}{||u|| = 1},
#' and then sets \eqn{v = w / \|w\|}{v = w / ||w||} and
#' \eqn{d = u^T x v}{d = u^T x v};
#' see the referenced paper for details. The penalty \eqn{\lambda}{lambda} is
#' implicitly determined from the specified desired number of nonzero values \code{n}.
#' Higher rank output is determined similarly
#' but using a sequence of \eqn{\lambda}{lambda} values determined to maintain the desired number
#' of nonzero elements in each column of \code{v} specified by \code{n}.
#' Unlike standard SVD or PCA, the columns of the returned \code{v} when \code{k > 1} may not be orthogonal.
#'
#' @return
#' A list containing the following components:
#' \itemize{
#'    \item{u} {regularized left singular vectors with orthonormal columns}
#'    \item{d} {regularized upper-triangluar projection matrix so that \code{x \%*\% v == u \%*\% d}} 
#'    \item{v} {regularized, sparse right singular vectors with columns of unit norm}
#'    \item{center, scale} {the centering and scaling used, if any}
#'    \item{lambda} {the per-column regularization parameter found to obtain the desired sparsity}
#'    \item{iter} {number of soft thresholding iterations}
#'    \item{n} {value of input parameter \code{n}}
#'    \item{alpha} {value of input parameter \code{alpha}}
#' }
#' @note
#' Our \code{ssvd} implementation of the Shen-Huang method makes the following choices:
#' \enumerate{
#' \item{The l1 penalty is the only available penalty function. Other penalties may appear in the future.}
#' \item{Given a desired number of nonzero elements in \code{v}, value(s) for the \eqn{\lambda}{lambda}
#'       penalty are determined to achieve the sparsity goal subject to the parameter \code{alpha}.}
#' \item{An experimental block implementation is used for results with rank greater than 1 (when \code{k > 1})
#'       instead of the deflation method described in the reference.}
#' \item{The choice of a penalty lambda associated with a given number of desired nonzero
#'       components is not unique. The \code{alpha} parameter, a scalar between zero and one,
#'       selects any possible value of lambda that produces the desired number of
#'       nonzero entries.}
#' \item{Our method returns an upper-triangular matrix \code{d} when \code{k > 1} so
#'       that \code{x \%*\% v == u \%*\% d}. Non-zero
#'       elements above the diagonal result from non-orthogonality of the \code{v} matrix,
#'       providing a simple interpretation of cumulative information, or explained variance
#'       in the PCA case, via the singular value decomposition of \code{d}.}
#' }
#'
#' What if you have no idea for values of the argument \code{n} (the desired sparsity)?
#' The reference describes a cross-validation and an ad-hoc approach; neither of which are
#' in the package yet. Both are prohibitively computationally expensive for matrices with a huge
#' number of columns. A future version of this package will include a revised approach to
#' automatically selecting a reasonable sparsity constraint.
#'
#' Compare with the similar but more general functions \code{SPC} and \code{PMD} in the \code{PMA} package
#' by Daniela M. Witten, Robert Tibshirani, Sam Gross, and Balasubramanian Narasimhan.
#' The \code{PMD} function can compute low-rank regularized matrix decompositions with sparsity penalties
#' on both the \code{u} and \code{v} vectors. The \code{ssvd} function is
#' similar to the PMD(*, L1) method invocation of \code{PMD} or alternatively the \code{SPC} function.
#' Although less general than \code{PMD},
#' the \code{ssvd} function can be faster and more memory efficient for the problems
#' that it can solve. See the examples below for more information.
#'
#' @references
#' \itemize{
#'   \item{Shen, Haipeng, and Jianhua Z. Huang. "Sparse principal component analysis via regularized low rank matrix approximation." Journal of multivariate analysis 99.6 (2008): 1015-1034.}
#'   \item{Witten, Tibshirani and Hastie (2009) A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis. _Biostatistics_ 10(3): 515-534.}
#' }
#' @examples
#'
#' set.seed(1)
#' u <- matrix(rnorm(200), ncol=1)
#' v <- matrix(c(runif(50, min=0.1), rep(0,250)), ncol=1)
#' u <- u / drop(sqrt(crossprod(u)))
#' v <- v / drop(sqrt(crossprod(v)))
#' x <- u %*% t(v) + 0.001 * matrix(rnorm(200*300), ncol=300)
#' s <- ssvd(x, n=50)
#' table(actual=v[, 1] != 0, estimated=s$v[, 1] != 0)
#' oldpar <- par(mfrow=c(2, 1))
#' plot(u, cex=2, main="u (black circles), Estimated u (blue discs)")
#' points(s$u, pch=19, col=4)
#' plot(v, cex=2, main="v (black circles), Estimated v (blue discs)")
#' points(s$v, pch=19, col=4)
#'
#' # Compare with SPC from the PMA package, regularizing only the v vector and choosing
#' # the regularization constraint `sum(abs(s$v))` computed above by ssvd
#' # (they find the about same solution in this "sparse SVD" case):
#' if (requireNamespace("PMA", quietly = TRUE)) {
#'   p <- PMA::SPC(x, sumabsv=sum(abs(s$v)), center=FALSE)
#'   table(actual=v[, 1] != 0, estimated=p$v[, 1] != 0)
#'   # compare optimized values
#'   print(c(ssvd=s$d, SPC=p$d))
#'
#'   # Same example, but computing a "sparse PCA", again about the same results:
#'   sp <- ssvd(x, n=50, center=TRUE)
#'   pp <- PMA::SPC(x, sumabsv=sum(abs(sp$v)), center=TRUE)
#'   print(c(ssvd=sp$d, SPC=pp$d))
#' }
#'
#'
#' # Let's consider a trivial rank-2 example (k=2) with noise. Like the
#' # last example, we know the exact number of nonzero elements in each
#' # solution vector of the noise-free matrix.
#' set.seed(1)
#' u <- qr.Q(qr(matrix(rnorm(400), ncol=2)))
#' v <- matrix(0, ncol=2, nrow=300)
#' v[sample(300, 15), 1] <- runif(15, min=0.1)
#' v[sample(300, 50), 2] <- runif(50, min=0.1)
#' v <- qr.Q(qr(v))
#' x <- u %*% (c(2, 1) * t(v)) + .01 * matrix(rnorm(200 * 300), 200)
#' s <- ssvd(x, k=2, n=c(15, sum(v[, 2] != 0)))
#'
#' # Compare actual and estimated vectors:
#' table(actual=v[, 1] != 0, estimated=s$v[, 1] != 0)
#' table(actual=v[, 2] != 0, estimated=s$v[, 2] != 0)
#' plot(v[, 1], cex=2, main="True v1 (black circles), Estimated v1 (blue discs)")
#' points(s$v[, 1], pch=19, col=4)
#' plot(v[, 2], cex=2, main="True v2 (black circles), Estimated v2 (blue discs)")
#' points(s$v[, 2], pch=19, col=4)
#' par(oldpar)
#'
#' @export
ssvd <- function(x, k=1, n=2, maxit=500, tol=1e-3, center=FALSE, scale.=FALSE, alpha=0, tsvd=NULL, ...)
{
  if (alpha < 0  || alpha >= 1) stop("0 <= alpha < 1")
  if (is.logical(center) && center) center <- colMeans(x)
  if (is.logical(scale.))
  {
    if (scale.)
    {
      if (is.numeric(args$center))
      { 
        f <- function(i) sqrt(sum((x[, i] - center[i]) ^ 2) / (nrow(x) - 1L))
        scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE)
      } else scale. <- apply(x, 2L, function(v) sqrt(sum(v ^ 2) / max(1, length(v) - 1L)))
    }
  }
  if (all(n > ncol(x) - 1))
  {
    warning("no sparsity constraints specified")
    return(irlba(x, k, ...))
  }
  n <- ncol(x) - n
  if (length(n) != k) n <- rep(n, length.out=k) # warn?
  s <- tsvd
  if (is.null(tsvd)) s <- irlba(x, k, scale=scale., center=center, ...)
  lambda <- c()
  soft <- function(x, u, p)
  {
    y <- crossprod(x, u)
    if (is.numeric(center)) y <- y - sum(u) * center
    if (is.numeric(scale.)) y <- y / scale.
    # apply a column-wise penalty
    a <- abs(y)
    z <- apply(a, 2, sort)
    lambda <<- vapply(seq(length(p)), function(j) (1 - alpha) * z[p[j], j] + alpha * z[p[j] + 1, j], pi, USE.NAMES=FALSE)
    sign(y) * pmax(sweep(a, 2, lambda, `-`), 0)
  }
  s$v <- s$d * s$v
  iter <- 1
  delta_u <- Inf
  while(delta_u > tol && iter < maxit)
  {
    u <- s$u
    s$v <- soft(x, s$u, n)
    if (is.numeric(scale.)) s$v <- s$v / scale.
    if (is.numeric(center)) s$u <- qr.Q(qr(x %*% s$v - drop(crossprod(center, s$v))))
    else s$u <- qr.Q(qr(x %*% s$v))
    delta_u <- sqrt(sum(apply(u - s$u, 2, crossprod)))
    iter <- iter + 1
  }
  if (iter >= maxit) warning("Maximum number of iterations reached before convergence: solution may not be optimal. Consider increasing 'maxit'.")
  s$v <- s$v %*% diag(1 / sqrt(apply(s$v, 2, crossprod)), ncol(s$v), ncol(s$v))
  d <- s$v
  if (is.numeric(scale.)) d <- d / scale.
  d1 <- x %*% d
  if (is.numeric(center)) d1 <- d1 - drop(crossprod(center, d))
  d <- crossprod(s$u, d1)
  list(u = s$u, v = s$v, d = d, iter = iter, lambda = lambda, center=center, scale=scale., n=n, alpha=alpha)
}


adhoc <- function(f, a, b, tol)
{
  c <- (b - a) / 2
  ans <- f(c)
  while(abs(ans) > tol)
  {
  }
}
