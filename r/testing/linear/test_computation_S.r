# --- Test the computation of S --------
library(Rcpp)
Rcpp::sourceCpp("../src/fit.cpp")


# ----------------------------------------
# Test data
# ----------------------------------------
set.seed(1)
n <- 100
p <- 1000
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-4, gsize), rep(8, gsize), rep(0, p - 3 * gsize))
y <- X %*% b + rnorm(n, 0, 1)


# ----------------------------------------
# R implementation
# ----------------------------------------
r_compute_S <- function(yty, yx, xtx, groups, mu, s, g, p)
{
    # S := E[ || y - Xb ||^2 ]
    #    = <y, y> - 2 <y, X E[b] > + E[<Xb, Xb>]
    res <- 0
    cij <- 0
    cijg <- 0
    cijng <- 0
    for (i in 1:p) {
	for (j in 1:p) {
	    if (i == j) {
		res <- res + xtx[i, i] * g[i] * (s[i]*s[i] + mu[i]*mu[i])
	    }

	    if ((i != j) && (groups[i] == groups[j])) {
		res <- res + xtx[i, j] * g[i] * mu[i] * mu[j]
	    }

	    if ((i != j) && (groups[i] != groups[j])) {
		res <- res + xtx[i, j] * g[i] * g[j] * mu[i] * mu[j]
	    }
	}
    }
    res <- yty - 2 * sum(yx * mu * g) + res
    return(res)
}


# ----------------------------------------
# Compare against C++ imp
# ----------------------------------------
yty <- sum(y * y)
yx <- t(y) %*% X
xtx <- t(X) %*% X

# generate random mu
set.seed(1)
mu <- rnorm(p)
s <- rgamma(p, 1, 1)
g <- rep(runif(p/5), each=5)


working <- approxEq(
    r_compute_S(yty, yx, xtx, groups, mu, s, g, p),
    compute_S(yty, yx, xtx, groups, mu, s, g, p)
)


if (working) message("TEST PASSED") else message("TEST FAILED")


# ----------------------------------------
# Simulation with MCI
# ----------------------------------------
reps <- replicate(1e5, {
    # sample from b ~ SpSL
    b <- (runif(p/5) < g[unique(groups) * 5])[groups] * rnorm(p, mu, s)

    l2(y - X %*% b)^2
})

approxEq(
    compute_S(yty, yx, xtx, groups, mu, s, g, p),
    mean(reps),
    e = 1e3
)
