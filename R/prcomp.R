#' Principal Components Analysis
#'
#' Efficient computation of a truncated principal components analysis of a given data matrix
#' using an implicitly restarted Lanczos method from the \code{\link{irlba}} package.
#'
#' @param x a numeric or complex matrix (or data frame) which provides
#'          the data for the principal components analysis.
#' @param retx a logical value indicating whether the rotated variables should be returned.
#' @param center a logical value indicating whether the variables should be
#'          shifted to be zero centered. Alternately, a centering vector of length
#'          equal the number of columns of \code{x} can be supplied.
#' @param scale. a logical value indicating whether the variables should be
#'          scaled to have unit variance before the analysis takes place.
#'          The default is \code{FALSE} for consistency with S, but scaling is often advisable.
#'          Alternatively, a vector of length equal the number of columns of \code{x} can be supplied.
#'
#'          The value of \code{scale} determines how column scaling is performed
#'          (after centering).  If \code{scale} is a numeric vector with length
#'          equal to the number of columns of \code{x}, then each column of \code{x} is
#'          divided by the corresponding value from \code{scale}.  If \code{scale} is
#'          \code{TRUE} then scaling is done by dividing the (centered) columns of
#'          \code{x} by their standard deviations if \code{center=TRUE}, and the
#'          root mean square otherwise.  If \code{scale} is \code{FALSE}, no scaling is done.
#'          See \code{\link{scale}} for more details.
#' @param n integer number of principal component vectors to return, must be less than
#' \code{min(dim(x))}.
#' @param ... additional arguments passed to \code{\link{irlba}}.
#'
#' @return
#' A list with class "prcomp" containing the following components:
#' \itemize{
#'    \item{sdev} {the standard deviations of the principal components (i.e.,
#'          the square roots of the eigenvalues of the
#'          covariance/correlation matrix, though the calculation is
#'          actually done with the singular values of the data matrix).}
#'   \item{rotation} {the matrix of variable loadings (i.e., a matrix whose columns
#'          contain the eigenvectors).}
#'   \item {x} {if \code{retx} is \code{TRUE} the value of the rotated data (the centred
#'          (and scaled if requested) data multiplied by the \code{rotation}
#'         matrix) is returned.  Hence, \code{cov(x)} is the diagonal matrix
#'          \code{diag(sdev^2)}.}
#'   \item{center, scale} {the centering and scaling used, or \code{FALSE}.}
#' }
#'
#' @note
#' The signs of the columns of the rotation matrix are arbitrary, and
#' so may differ between different programs for PCA, and even between
#' different builds of R.
#'
#' NOTE DIFFERENCES WITH THE DEFAULT \code{\link{prcomp}} FUNCTION!
#' The \code{tol} truncation argument found in \code{prcomp} is not supported.
#' In place of the truncation tolerance in the original function, the
#' \code{prcomp_irlba}  function has the argument \code{n} explicitly giving the
#' number of principal components to return. A warning is generated if the
#' argument \code{tol} is used, which is interpreted differently between
#' the two functions.
#'
#' @examples
#' set.seed(1)
#' x  <- matrix(rnorm(200), nrow=20)
#' p1 <- prcomp_irlba(x, n=3)
#' summary(p1)
#'
#' # Compare with
#' p2 <- prcomp(x, tol=0.7)
#' summary(p2)
#'
#'
#' @seealso \code{\link{prcomp}}
#' @import Matrix
#' @importFrom stats rnorm prcomp sd var
#' @importFrom methods slotNames slot
#' @export
prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE, ...)
{
  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  if ("tol" %in% a)
    warning("The `tol` truncation argument from `prcomp` is not supported by
`prcomp_irlba`. If specified, `tol` is passed to the `irlba` function to
control that algorithm's convergence tolerance. See `?prcomp_irlba` for help.")
# Try to convert data frame to matrix...
  if (is.data.frame(x)) x <- as.matrix(x)
  args <- list(A=x, nv=n)
  if (is.logical(center))
  {
    if (center) args$center <- colMeans(x)
  } else args$center <- center
  if (is.logical(scale.))
  {
      if (is.numeric(args$center))
      {
        f <- function(i) sqrt(sum((x[, i] - args$center[i]) ^ 2) / (nrow(x) - 1L))
        scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE)
        if (ans$scale) ans$totalvar <- ncol(x)
        else ans$totalvar <- sum(scale. ^ 2)
      } else
      {
        if (ans$scale)
        {
          scale. <- apply(x, 2L, function(v) sqrt(sum(v ^ 2) / max(1, length(v) - 1L)))
          f <- function(i) sqrt(sum((x[, i] / scale.[i]) ^ 2) / (nrow(x) - 1L))
          ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE) ^ 2)
        } else
        {
          f <- function(i) sum(x[, i] ^ 2) / (nrow(x) - 1L)
          ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE))
        }
      }
      if (ans$scale) args$scale <- scale.
  } else
  {
    args$scale <- scale.
    f <- function(i) sqrt(sum((x[, i] / scale.[i]) ^ 2) / (nrow(x) - 1L))
    ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE))
  }
  if (!missing(...)) args <- c(args, list(...))

  s <- do.call(irlba, args=args)
  ans$sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  ans$center <- args$center
  if (retx)
  {
    ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN=`*`)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}

#' Summary method for truncated pca objects computed by \code{prcomp_irlba}.
#' @param object An object returned by \code{prcomp_irlba}.
#' @param ... Optional arguments passed to \code{summary}.
#' @method summary irlba_prcomp
#' @export
summary.irlba_prcomp <- function(object, ...)
{
  chkDots(...)
  vars <- object$sdev ^ 2
  vars <- vars / object$totalvar
  importance <- rbind("Standard deviation" = object$sdev,
                      "Proportion of Variance" = round(vars, 5),
                      "Cumulative Proportion" = round(cumsum(vars), 5))
  k <- ncol(object$rotation)
  colnames(importance) <- c(colnames(object$rotation), rep("", length(vars) - k))
  object$importance <- importance
  class(object) <- "summary.prcomp"
  object
}
