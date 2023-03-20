library(MASS)
n = 100
Xx = mvrnorm(n,rep(0,100),diag(100))
hist(Xx[,1])
cor(Xx[,2],Xx[,6])
cor(Xx[2,],Xx[3,])
min(cor(Xx)-diag(100))
max(cor(t(Xx))-diag(100))

Xy = matrix(rnorm(n*100),n,100)
hist(Xy[,1])
hist(Xy[1,])
cor(Xy[,2],Xy[,3])
cor(Xy[1,],Xy[3,])
min(cor(Xy)-diag(100))
max(cor(t(Xy))-diag(100))

sourcecode_mvrnorm = function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE, 
                               EISPACK = FALSE) 
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  if (EISPACK) 
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  print(eS$vectors)
  print(diag(sqrt(pmax(ev, 0)), p))
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  print(X)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1) 
    drop(X)
  else t(X)
}
Xx = sourcecode_mvrnorm(3,rep(0,4),diag(4))
print(Xx)
