require("irlba")

VERB = FALSE
TOL = 1e-5
set.seed(1)
N = 20

norm2 = irlba:::norm2
cross = irlba:::cross

check <- function(x, s, N, ht=head) {
  i <- ht(seq(1, ncol(s$u)), N)
  stopifnot(norm2(x$d - ht(s$d, N)) < TOL)
  stopifnot(norm2(abs(diag(cross(x$u, s$u[,i, drop=FALSE]))) - 1) < TOL)
  stopifnot(norm2(abs(diag(cross(x$v, s$v[,i, drop=FALSE]))) - 1) < TOL)
  cat("OK\n", file=stderr())
}

x = matrix(rnorm(N*(N-1)), ncol=N)
s = svd(x)
cat("Testing: dense real largest          \t", file=stderr())
l = irlba:::irlba(x, nv=3, verbose=VERB)
check(l, s, 3, head)

# deflate
cat("Testing: dense real largest deflate \t", file=stderr())
l2 = irlba:::irlba(x, nv=3, v=l, verbose=VERB)
check(l2, s, 6, head)

# smallest
cat("Testing: dense real smallest        \t", file=stderr())
l = irlba:::irlba(x, nv=3, smallest=TRUE, verbose=VERB)
check(l, s, 3, tail)

# deflate
cat("Testing: dense real smallest deflate\t", file=stderr())
l2 = irlba:::irlba(x, nv=3, v=l, smallest=TRUE, verbose=VERB)
check(l2, s, 6, tail)

cat("Testing: dense complex largest      \t", file=stderr())
x = x + matrix(rnorm(N*(N-1)) + 1i, ncol=N)
s = svd(x)
l = irlba:::irlba(x, nv=3, verbose=VERB)
check(l, s, 3, head)

cat("Testing: dense complex largest defl \t", file=stderr())
l2 = irlba:::irlba(x, nv=3, v=l, verbose=VERB)
check(l2, s, 6, head)


# Matrix-classed
require(Matrix)
N = 20
cat("Testing: Matrix real largest        \t", file=stderr())
x = Matrix(matrix(rnorm(N*(N-1)), ncol=N))
s = svd(x)
l = irlba:::irlba(x, nv=3, verbose=VERB)
check(l, s, 3, head)

cat("Testing: Matrix real largest defl   \t", file=stderr())
l2 = irlba:::irlba(x, nv=3, v=l, verbose=VERB)
check(l2, s, 6, head)

# smallest
cat("Testing: Matrix real smallest       \t", file=stderr())
l = irlba:::irlba(x, nv=3, smallest=TRUE, verbose=VERB)
check(l, s, 3, tail)

cat("Testing: sparse real smallest defl  \t", file=stderr())
l2 = irlba:::irlba(x, nv=3, v=l, smallest=TRUE, verbose=VERB)
check(l2, s, 6, tail)

# Matrix does not (yet) support complex values :(

# test for https://github.com/bwlewis/irlba/issues/47
cat("GitHub issue #47                    \t", file=stderr())
set.seed(2345)
a <- spMatrix(50, 40, x=runif(200), i=sample(50, 200, replace=TRUE), j=sample(40, 200, replace=TRUE))
center <- runif(ncol(a))
scale <- runif(ncol(a))
L <- irlba(a, 5, scale=scale, center=center)
S <- svd(scale(a, center=center, scale=scale))
check(L, S, 5, head)

cat("Testing center=TRUE                \t", file=stderr())
x = matrix(rnorm(N*N), ncol=N)
s = svd(sweep(x, 2, colMeans(x), FUN=`-`))
l = irlba:::irlba(x, nv=3, nu=NULL, verbose=VERB, center=TRUE)
check(l, s, 3, head)

