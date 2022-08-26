library(gsvb)

# ----------------------------------------
# Example 1
# ----------------------------------------
n <- 100
p <- 1000
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rnorm(gsize, 4, 2), rep(0, gsize * 4),
       rnorm(gsize, -4, 2), rep(0, p - 7 * gsize))
y <- X %*% b + rnorm(n, 0, 1)


f <- gsvb.fit(y, X, groups, lambda=0.5, tol=1e-3)
gsvb.elbo(f, y, X, groups)
f <- gsvb.fit(y, X, groups, diag_covariance=FALSE, lambda=0.5, tol=1e-3, track_elbo=TRUE, track_elbo_every=1)
gsvb.elbo(f, y, X, groups)

plot(f$beta_hat, col=4, ylab=expression(hat(beta)))
points(b, pch=20)
plot(f$elbo, type="l", lwd=2)


# ----------------------------------------
# Example 2
# ----------------------------------------
b <- 15

x <- runif(100, -1, 1)
y <- cos(x)  - 2 * sin(2*x) + rnorm(100, 0, 0.5)

plot(x, y)
curve(cos(x)  -2 * sin(2*x), -1, 1, add=T)

X <- t(sapply(x, function(x) c(cos(1:b * x), sin(1:b * x))))
groups <- rep(1:30, each=1)
f <- gsvb.fit(y, X, groups, intercept=TRUE, diag_covariance=FALSE, niter=150, tau_a0=1, tau_b0=1)
f <- gsvb.fit(y, X, groups, intercept=FALSE, diag_covariance=TRUE, niter=150, tau_a0=1, tau_b0=1)
f <- gsvb.fit(y, X, groups, intercept=FALSE, diag_covariance=FALSE, niter=150, tau_a0=1, tau_b0=1)

fx <- function(x) sum(f$b[1:b] * cos(1:b * x) + f$b[(b+1):(2*b)] * sin(1:b * x))

xs <- seq(-1, 1, by=0.01)
plot(xs, sapply(xs, fx), type="l", ylim=c(-2, 3), lwd=3, lty=2)
curve(cos(x) - 2 * sin(2*x), -1, 1, add=T, lwd=2, col=2)
points(x, y)


# ----------------------------------------
# Example 3
# ----------------------------------------
# no longer group sparse
n <- 100
p <- 1000
gsize <- 1
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rnorm(gsize, 4, 2), rep(0, gsize * 4),
       rnorm(gsize, -4, 2), rep(0, p - 7 * gsize))
y <- X %*% b + rnorm(n, 0, 1)


f <- gsvb.fit(y, X, groups, FALSE, lambda=0.5, tol=1e-3)
f <- gsvb.fit(y, X, groups, TRUE, diag_covariance=TRUE, lambda=0.5, tol=1e-3, 
	      track_elbo=TRUE, track_elbo_every=1)
plot(f$beta_hat, col=4, ylab=expression(hat(beta)))
points(b, pch=20)
