# ----------------------------------------
# Test data
# ----------------------------------------
set.seed(1)
n <- 150
p <- 500
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-0.5, gsize), rep(0.2, gsize), rep(0, p - 3 * gsize))
xb <- X %*% b
y <- rpois(n, exp(xb))

# glm(y ~ X[ , !!b], family=poisson)

gsvb.pois(y, X, groups)


# ----------------------------------------
# algorithm
# ----------------------------------------
gsvb.pois <- function(y, X, groups, niter=500, fit=NULL, lambda=1)
{
    n <- nrow(X) 
    p <- ncol(X)
    ugroups <- unique(groups)

    # initialize
    m <- rnorm(p)
    theta <- list()
    S <- list()
    g <- +!!b

    for (i in seq_along(ugroups)) 
    {
	G <- which(groups == ugroups[i])
	mk <- length(G)
	pars <- runif(n + 1, 0.1, 0.2)

	theta[[i]] <- pars
	S[[i]] <- solve(t(X[ , G]) %*% diag(pars[1:n]) %*% X[ , G] + 
			    pars[n + 1] * diag(rep(1, mk)))
    }

    P <- compute_P(X, m, S, g, ugroups)

    # main loop
    for (iter in 1:niter)
    {
	for (gi in seq_along(ugroups))
	{
	    G <- which(groups == ugroups[gi])
	    mk <- length(G)

	    P <- P / compute_P_G(X, m, S, g, G, gi)

	    m[G] <- optim(m[G], 
		fn=function(m_G) opt_mu(m_G, y, X, S[[gi]], G, lambda, P),
		control=list(maxit=20),
		method="BFGS")$par
	    
	    theta[[gi]] <- optim(theta[[gi]], 
		fn=function(theta) opt_S(theta, X, m, G, lambda, P),
		control=list(maxit=20),
		method="BFGS")$par
	    
	    pars <- theta[[gi]]
	    S[[gi]] <- solve(t(X[ , G]) %*% diag(as.numeric(P) * exp(pars[1:n])) %*% X[ , G] + 
				pars[n + 1] * diag(rep(1, mk)))

	    print(log(det(S[[gi]])))

	    P <- P * compute_P_G(X, m, S, g, G, gi)
	}
	cat("\niter: ", iter)
    }

    return(list(m=m, S=S, g=g))
}


# ----------------------------------------
# Update mu
# ----------------------------------------
opt_mu <- function(m_G, y, X, S_G, G, lambda, P) 
{
    PP <- sum(P * mvmgf(X[ , G], m_G, S_G))

    res <- - sum(y * X[ , G] %*% m_G) +
    PP +
    lambda * (sum(diag(S_G) + m_G^2))^(1/2)

    return(res)
}


# ----------------------------------------
# ypdate S
# ----------------------------------------
opt_S <- function(theta, X, m, G, lambda, P)
{
    mk <- length(G)
    A <- diag(as.numeric(P) * exp(theta[1:n]))
    gam <- theta[n+1]

    S <- solve(t(X[ , G]) %*% A %*% X[ , G] + gam * diag(rep(1, mk)))

    PP <- sum(P * mvmgf(X[ , G], m[G], S))

    res <- PP - 
    0.5 * log(det(S)) +
    lambda * (sum(diag(S)+ m[G]^2))^(1/2)

    return(res)
}


# ----------------------------------------
# update g
# ----------------------------------------
opt_g <- function(y, X, XAX, m, s, g, G, Gc, lambda, w) 
{
    mk <- length(G)
    Ck <- mk * log(2) + (mk -1)/2 * log(pi) + lgamma( (mk + 1) / 2)

    res <- 
	log(w / (1 - w)) + 
	0.5 * mk -
	Ck +
	mk * log(lambda) +
	0.5 * sum(log(2 * pi * s[G]^2)) -
	lambda * sqrt(sum(s[G]^2) + sum(m[G]^2)) +
	sum((y - 0.5) * X[ , G] %*% m[G]) -
	0.5 * t(m[G]) %*% XAX[G, G] %*% m[G] -
	0.5 * sum(diag(XAX[G, G]) * s[G]^2) -
	t(m[G]) %*% XAX[G, Gc] %*% (g[Gc] * m[Gc])

    return(sigmoid(res))
}


# ---------------------------------------
# misc funcs
# ---------------------------------------
compute_P <- function(X, m, S, g, ugroups) 
{
    P <- rep(1, nrow(X))

    for (i in seq_along(ugroups)) {
	group <- ugroups[i]
	G <- which(groups == group)
	P <- P * compute_P_G(X, m, S, g, G, i)
    }
    return(P)
}


compute_P_G <- function(X, m, S, g, G, gi)
{
    (1-g[G[1]]) + g[G[1]] * mvmgf(X[, G], m[G], S[[gi]])
}


mvmgf <- function(X, mu, S) 
{
    exp(X %*% mu + 0.5 * diag(X %*% S %*% t(X)))
}

