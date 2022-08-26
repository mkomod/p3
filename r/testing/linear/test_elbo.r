library(Rcpp)

Rcpp::sourceCpp("../src/fit.cpp")

# ----------------------------------------
# Test data
# ----------------------------------------
n <- 100
p <- 1000
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-4, gsize), rep(8, gsize), rep(0, p - 3 * gsize))
y <- X %*% b + rnorm(n, 0, 1)

yty <- sum(y * y)
yx <- t(X) %*% y
xtx <- t(X) %*% X

# ----------------------------------------
# Fit model
# ----------------------------------------
f <- fit(y, X, groups, 1, 1, p/gsize, 1, 1, runif(p, -0.2, 0.2), rep(0.2, p),
	 rep(0.5, p), T, 1, 500, 50, 1e-3, T)

elbo(yty, yx, xtx, groups, n, p, f$m, f$s, f$g, f$tau_a, f$tau_b, 
     1, 1, p/gsize, 1, 1, 500, F, 1e-3)

f <- gsvb.fit(y, X, groups, intercept=TRUE, track_elbo_every=1, track_elbo=F)

gsvb.elbo(f, y, X, groups)


lmvde <- function(x, lambda, m) 
{
    # log multivariate double exp
    -m*log(2) - 0.5*(m-1)*log(pi) - lgamma(0.5*(m+1)) + m*log(lambda) -
	lambda * sqrt(sum(x^2))
}

# ----------------------------------------
# ELBO test func
# ----------------------------------------
elbor <- function(f, y, X, groups, mcn=1e3) 
{
    n <- nrow(X)
    p <- ncol(X)

    xtx <- t(X) %*% X

    attach(f); attach(parameters)
    sigma_nsq <- sigma^-2

    res <- - 0.5 * n * log(2 * pi * sigma^2) -
	0.5 * sigma_nsq * sum(y^2) +
	sigma_nsq * t(y) %*% (X %*% (g * mu)) +
	0.5 * sum(g * log(2 * pi * s^2))

    a <- 0
    for (i in 1:p) {
	for (j in 1:p) {
	    if (i == j) {
		a <- a + xtx[j, i] * g[j] * (s[j]^2 + mu[j]^2) 
	    }
	    if (i != j && groups[i] == groups[j]) {
		a <- a + xtx[j, i] * g[j] * mu[j] * mu[i]
	    }
	    if (i != j && groups[i] != groups[j]) {
		a <- a + xtx[j, i] * g[j] * g[i] * mu[j] * mu[i]
	    }
	}
    }
    
    r <- 0
    for (gi in unique(groups)) {
	G <- which(groups == gi)
	mk <- length(G)
	gk <- g[G[1]]

	mci <- mean(replicate(mcn, lambda * sqrt(sum(rnorm(mk, mu[G], s[G])^2))))

	Ck = -mk*log(2) - 0.5*(mk-1)*log(pi) - lgamma(0.5*(mk+1))
	r <- r + gk * (
	    Ck + 
	    0.5*mk + 
	    mk*log(lambda) - 
	    mci - 
	    log((gk + 1e-8)/(w +1e-8))
	)
	r <- r - (1-gk) * log((1-gk+1e-8)/(1-w+1e-8))
    }

    detach(f); detach(parameters)
    res - 0.5 * sigma_nsq * a  + r
}

# ----------------------------------------
# ELBO test func 2
# ----------------------------------------
elbo.mci <- function(f, y, X, groups, mcn=1e3)
{
    n <- nrow(X)
    p <- ncol(X)
    xtx <- crossprod(X, X)

    attach(f); attach(parameters)
    w <- a0 / (a0 + b0)

    # non stochasitc part of the ELBO
    res <- - 0.5*n*log(2*pi*sigma^2) - 0.5 * sigma^-2 * sum(y^2)
    
    gs <- sapply(unique(groups), function(i) which(groups == i)[1])
    ugroups <- unique(groups)

    sigma_nsq <- sigma^-2
    
    mci <- replicate(mcn, {
	gg <- (runif(length(gs)) < g[gs])[groups]
	b <- rnorm(p, mu, s)
	bg <- b * gg

	r <- - 0.5 * sigma_nsq * sum(xtx * outer(bg, bg)) +
	    sigma_nsq * t(y) %*% (X %*% bg)
	
	a <- 0.0
	for (gi in gs) {
	    if (runif(1) > g[gi]) {
		a <- a + log((1-g[gi]) / (1-w))
	    } else {
		G <- which(groups[gi] == groups)
		b <- rnorm(length(G), mu[G], s[G])
		a <- a + log(w) + lmvde(b, lambda, length(G))
		a <- a - log(g[gi]) - sum(dnorm(b, mu[G], s[G], log=T))
	    }
	}
	as.vector(r) + a
    })

    detach(f); detach(parameters) 
    
    res + mean(mci)
}


# ----------------------------------------
# Testing the different ELBOs
# ----------------------------------------
elbo(y, X, groups, f$m, f$s, f$g, 1, 1, 200, 1, 10000)
elbor(f, y, X, groups, 1000)
elbo.mci(f, y, X, groups, 1e4)

