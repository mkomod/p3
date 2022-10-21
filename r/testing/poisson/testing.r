# ----------------------------------------
# Test data
# ----------------------------------------
set.seed(1)
n <- 250
p <- 500
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-0.5, gsize), rep(0.2, gsize), rep(0, p - 3 * gsize))
xb <- X %*% b
y <- rpois(n, exp(xb))

# glm(y ~ X[ , !!b], family=poisson)

f <- gsvb.pois(y, X, groups, niter=20)
plot(f$m * f$g)
plot(f$g)
plot(f$s)


# ----------------------------------------
# algorithm
# ----------------------------------------
gsvb.pois <- function(y, X, groups, niter=500, fit=NULL, lambda=1)
{
    n <- nrow(X) 
    p <- ncol(X)
    ugroups <- unique(groups)

    # initialize
    m <- rnorm(p, 0, 0.5)
    s <- runif(p, 0.1, 0.2)
    # g <- +!!b
    g <- rep(0.1, p)
    w <- 1/(1+length(ugroups))

    P <- compute_P(X, m, s, g, ugroups)

    # main loop
    for (iter in 1:niter)
    {
	for (gi in seq_along(ugroups))
	{
	    G <- which(groups == ugroups[gi])
	    mk <- length(G)

	    P <- P / compute_P_G(X, m, s, g, G)

	    m[G] <- optim(m[G], 
		fn=function(m_G) opt_mu(m_G, y, X, s, G, lambda, P),
		control=list(maxit=20),
		method="BFGS")$par
	    
	    s[G] <- optim(s[G],
		fn=function(s_G) opt_s(s_G, X, m, G, lambda, P),
		control=list(maxit=20),
		method="L-BFGS-B", lower=s[G], upper=s[G]+0.2)$par

	    g[G] <- opt_g(y, X, m, s, g, G, lambda, w, P)
	    
	    P <- P * compute_P_G(X, m, s, g, G)
	}
	cat("\niter: ", iter)
	plot(g)
    }

    return(list(m=m, s=s, g=g))
}


# ----------------------------------------
# Update mu
# ----------------------------------------
opt_mu <- function(m_G, y, X, s, G, lambda, P) 
{
    PP <- sum(P * mvmgf(X[ , G], m_G, s[G]))

    res <- - sum(y * X[ , G] %*% m_G) +
    PP +
    lambda * (sum(diag(s[G]) + m_G^2))^(1/2)

    return(res)
}


# ----------------------------------------
# ypdate S
# ----------------------------------------
opt_s <- function(s_G, X, m, G, lambda, P)
{
    PP <- sum(P * mvmgf(X[ , G], m[G], s_G))

    res <- PP - 
    sum(log(s_G)) +
    lambda * (sum(s_G^2 + m[G]^2))^(1/2)

    return(res)
}


# ----------------------------------------
# update g
# ----------------------------------------
opt_g <- function(y, X, m, s, g, G, lambda, w, P)
{
    mk <- length(G)
    Ck <- mk * log(2) + (mk -1)/2 * log(pi) + lgamma( (mk + 1) / 2)
    P1 <- mvmgf(X[ , G], m[G], s[G])

    res <- 
	log(w / (1- w)) + 
	0.5 * mk - 
	Ck +
	mk * log(lambda) +
	0.5 * sum(log(2 * pi * s[G]^2)) -
	lambda * sqrt(sum(s[G]^2) + sum(m[G]^2)) +
	sum(y * X[ , G] %*% m[G]) -
	sum(P*(P1 - 1))

    sigmoid(res)
}


# ---------------------------------------
# misc funcs
# ---------------------------------------
compute_P <- function(X, m, s, g, ugroups) 
{
    P <- rep(1, nrow(X))

    for (group in ugroups) {
	G <- which(groups == group)
	P <- P * compute_P_G(X, m, s, g, G)
    }
    return(P)
}


compute_P_G <- function(X, m, s, g, G)
{
    (1-g[G[1]]) + g[G[1]] * mvmgf(X[, G], m[G], s[G])
}


mvmgf <- function(X, mu, s) 
{
    exp(X %*% mu + 0.5 * diag(X %*% diag(s^2) %*% t(X)))
}

sigmoid <- function(x) 1/(1+exp(-x))
