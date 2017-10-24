require("irlba")

# prcomp convenience function
x  <- matrix(rnorm(200), nrow=20)
p1 <- prcomp_irlba(x, n=3)
p2 <- prcomp(x, tol=0.7)
if (!isTRUE(all.equal(p1$sdev[1:2], p2$sdev[1:2])))
{
  stop("Failed basic prcomp test")
}

s <- summary(p1)

# scaling bug identified in issue #21
normalize_signs <- function(X, Y) {
  for (i in 1:ncol(X)) {
    if (sign(X[1, i]) != sign(Y[1, i])) {
      Y[, i] <- -Y[, i]
    }
  }
  return(Y)
}

all.equal_pca <- function(X, Y) {
  Y <- normalize_signs(X, Y)
  return(all.equal(X, Y, check.attributes=F, tolerance=1e-4))
}

set.seed(1)
X <- matrix(rnorm(2000), ncol=40)
M <- 5 # number of PCA components
centers <- colMeans(X)
sds <- apply(X, 2, sd)
rms <- apply(X, 2, function(x) sqrt(sum(x^2) / (length(x) - 1)))
Xc <- sweep(X, 2, centers, `-`)
Xs <- sweep(X, 2, sds, `/`)
Xcs <- sweep(Xc, 2, sds, `/`)
Xrms <- sweep(X, 2, rms, `/`)

# unscaled
scaled <- FALSE
centered <- FALSE
pca <- prcomp(X, center=centered, scale.=scaled)
sv <- svd(X)
svir <- irlba(X, nv=M, nu=M)
pcair <- prcomp_irlba(X, n=M, center=centered, scale.=scaled)
Xpca <- predict(pca)[, 1:M]
Xsvl <- sv$u[, 1:M] %*% diag(sv$d[1:M])
Xsvr <- X %*% sv$v[, 1:M]
Xsvirl <- svir$u %*% diag(svir$d)
Xsvirr <- X %*% svir$v
Xpcair <- predict(pcair)
Xpcair2 <- X %*% pcair$rotation

if (! isTRUE(all.equal_pca(Xsvl, Xsvr)) &&
     isTRUE(all.equal_pca(Xpca, Xsvl)) &&
     isTRUE(all.equal_pca(Xsvirl, Xsvirr)) &&
     isTRUE(all.equal_pca(Xpca, Xsvirl)) &&
     isTRUE(all.equal_pca(Xpcair, Xpcair2)) &&
     isTRUE(all.equal_pca(Xpca, Xpcair)) &&
     isTRUE(all.equal_pca(Xpcair, Xsvirl)))
{
  stop("failed unscaled, uncentered prcomp")
}

# scaled, uncentered
scaled <- TRUE
centered <- FALSE
pca <- prcomp(X, center=centered, scale.=scaled)
sv <- svd(Xrms)
svir <- irlba(X, nv=M, nu=M, scale=rms)
pcair <- prcomp_irlba(X, n=M, center=centered, scale.=scaled)

Xpca <- predict(pca)[, 1:M]
Xsvl <- sv$u[, 1:M] %*% diag(sv$d[1:M])
Xsvr <- Xrms %*% sv$v[, 1:M]
Xsvirl <- svir$u %*% diag(svir$d)
Xsvirr <- Xrms %*% svir$v
Xpcair <- predict(pcair)
Xpcair2 <- Xrms %*% pcair$rotation

if (! isTRUE(all.equal_pca(Xsvl, Xsvr)) &&
      isTRUE(all.equal_pca(Xpca, Xsvl)) &&
      isTRUE(all.equal_pca(Xsvirl, Xsvirr)) &&
      isTRUE(all.equal_pca(Xpca, Xsvirl)) &&
      isTRUE(all.equal_pca(Xpcair, Xpcair2)) &&
      isTRUE(all.equal_pca(Xpca, Xpcair)) &&
      isTRUE(all.equal_pca(Xpcair, Xsvirl)))
{
  stop("failed scaled, uncentered prcomp")
}


# issue #25 prcomp_irlba regression (error in scale. handling)
set.seed(1)
x <- matrix(rnorm(100), 10)
p <- prcomp_irlba(x, 3, scale.=TRUE, fastpath=FALSE)
p <- prcomp_irlba(x, 3, scale.=TRUE, fastpath=TRUE)
