#' Find a few approximate largest eigenvalues and corresponding eigenvectors of a symmetric matrix.
#'
#' Use \code{partial_eigen} to estimate a subset of the largest (most positive)
#' eigenvalues and corresponding eigenvectors of a symmetric dense or sparse
#' real-valued matrix.
#'
#' @param x numeric real-valued dense or sparse matrix.
#' @param n number of largest eigenvalues and corresponding eigenvectors to compute.
#' @param symmetric \code{TRUE} indicates \code{x} is a symmetric matrix (the default); specify \code{symmetric=FALSE} to compute the largest eigenvalues and corresponding eigenvectors of \code{t(x) \%*\% x} instead.
#' @param ... optional additional parameters passed to the \code{irlba} function.
#'
#' @return
#' Returns a list with entries:
#' \itemize{
#'   \item{values}{ n approximate largest eigenvalues}
#'   \item{vectors}{ n approximate corresponding eigenvectors}
#' }
#'
#' @note
#' Specify \code{symmetric=FALSE} to compute the largest \code{n} eigenvalues
#' and corresponding eigenvectors of the symmetric matrix cross-product
#' \code{t(x) \%*\% x}.
#'
#' This function uses the \code{irlba} function under the hood. See \code{?irlba}
#' for description of additional options, especially the \code{tol} parameter.
#'
#' @references
#' Augmented Implicitly Restarted Lanczos Bidiagonalization Methods, J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005.
#'
#' @examples
#' set.seed(1)
#' # Construct a symmetric matrix with some positive and negative eigenvalues:
#' V <- qr.Q(qr(matrix(runif(100),nrow=10)))
#' x <- V %*% diag(c(10, -9, 8, -7, 6, -5, 4, -3, 2, -1)) %*% t(V)
#' partial_eigen(x, 3)$values
#'
#' # Compare with eigen
#' eigen(x)$values[1:3]
#'
#' @seealso \code{\link{eigen}}, \code{\link{irlba}}
#' @export
partial_eigen <- function(x, n=5, symmetric=TRUE, ...)
{
  if (n > 0.5 * min(nrow(x),ncol(x)))
  {
    warning("You're computing a large percentage of total eigenvalues, the standard eigen function will likely work better!")
  }
  if (!symmetric)
  {
    L <- irlba(x, n, ...)
    return(list(vectors=L$v, values=L$d ^ 2))
  }
  L <- irlba(x, n, ...)
  s <- sign(L$u[1,] * L$v[1,])
  if (all(s > 0))
  {
    return(list(vectors=L$u, values=L$d))
  }
  i <- min(which(s < 0))
  shift <- L$d[i]
  L <- irlba(x, n, shift=shift, ...)
  return(list(vectors=L$u, values=L$d - shift))
}
