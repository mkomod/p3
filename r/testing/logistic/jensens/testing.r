library(Rcpp)

Rcpp::sourceCpp("./jensens.cpp")

approxEq <- function(a, b, e=1e-3) all(abs(a - b) < e)

# ----------------------------------------
# Test data
# ----------------------------------------
set.seed(1)
n <- 350
p <- 500
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-1, gsize), rep(1, gsize), rep(0, p - 3 * gsize))
xb <- X %*% b
prob <- 1 / (1 + exp(-xb))
y <- rbinom(n, 1, prob)


# ----------------------------------------
# Test alg
# ----------------------------------------
gsvb.logistic <- function(y, X, groups, m, s, g, niter=500, fit=NULL)
{
    n <- nrow(X) 
    p <- ncol(X)
    w <- 1 / p * length(unique(groups))
    
    P <- compute_P(X, m, s, g, groups)

    for (iter in 1:niter)
    {
	m.old <- m; s.old <- s; g.old <- g

	for (group in unique(groups))
	{
	    G <- which(groups == group) - 1

	    P <- P / compute_P_G(X, m, s, g, G)

	    m[G+1] <- jen_update_mu(y, X, m, s, 1, G, P)
	    s[G+1] <- jen_update_s( y, X, m, s, 1, G, P)
	    g[G+1] <- jen_update_g( y, X, m, s, 1, w, G, P)

	    P <- P * compute_P_G(X, m, s, g, G)
	}
	cat(iter)

	if (all(c(abs(m.old - m) < 1e-3, 
		  abs(s.old - s) < 1e-3, 
		  abs(g.old - g) < 1e-3))) 
	    break
    }
    
    return(list(m=m, s=s, g=g))
}

set.seed(1)
mu <- rnorm(p); s <- rep(0.1, p); g <- runif(p/gsize, 0.1, 0.2)[groups]
fit <- gsvb.logistic(y, X, groups-1, mu, s, g, niter=300)

plot(fit$m)
plot(fit$s)
plot(fit$g)


# ----------------------------------------
# Test update mu
# ----------------------------------------
mu <- rnorm(p); s <- rep(0.1, p); g <- runif(p/gsize)[groups]
G <- which(groups == 1) - 1

P <- compute_P(X, mu, s, g, groups) 
P <- P / compute_P_G(X, mu, s, g, G)

jen_update_mu(y, X, mu, s, 1, G, P)

optim(mu[G+1], 
    fn=function(mG) opt_mu(mG, y, X, mu, s, g, G+1, 1, P), 
    control=list(maxit=200),
    method="BFGS")$par


# ----------------------------------------
# Test update s
# ----------------------------------------
mu <- rnorm(p); s <- rep(0.1, p); g <- runif(p/gsize, 0.1, 0.2)[groups]
G <- which(groups == 1) - 1

P <- compute_P(X, mu, s, g, groups) 
P <- P / compute_P_G(X, mu, s, g, G)

jen_update_s(y, X, mu, s, 1, G, P)

optim(s[G+1], 
    fn=function(sG) opt_s(sG, y, X, mu, s, g, G+1, 1, P),
    control=list(maxit=200), method="BFGS")


# ----------------------------------------
# Test update g
# ----------------------------------------
mu <- rnorm(p); s <- rep(0.1, p); g <- runif(p/gsize, 0.1, 0.2)[groups]
G <- which(groups == 1) - 1

P <- compute_P(X, mu, s, g, groups) 
P <- P / compute_P_G(X, mu, s, g, G)

jen_update_g(y, X, mu, s, 1, 1/200, G, P)
opt_g(y, X, mu, s, g, G+1, 1, 1/200, P)


# ----------------------------------------
# Test compute P
# ----------------------------------------
mu <- rnorm(p); s <- rep(0.1, p); g <- runif(p/gsize)[groups]
G <- which(groups == 1)

compute_P_G(X, mu, s, g, G-1)[1]
compute_S_G(X, mu, s, g, G)[1]

approxEq(compute_P(X, b, rep(.1, p), rep(0.5, p), groups),
	 compute_S(X, b, rep(.1, p), rep(0.5, p), groups))

microbenchmark::microbenchmark(
    compute_P(X, b, rep(.1, p), rep(0.5, p), groups)[1],
    compute_S(X, b, rep(.1, p), rep(0.5, p), groups)[1]
)


# ----------------------------------------
# MVN MGF
# ----------------------------------------
n_mgf <- function(X, m, s)
{
    apply(X, 1, function(x) {
	exp(sum(x * m + 0.5 * x^2 * s^2))
    })
}

approxEq(mvnMGF(X, b, rep(.1, p)), n_mgf(X, b, rep(.1, p)))


test_minus(matrix(1:10, ncol=2), 1:5)


