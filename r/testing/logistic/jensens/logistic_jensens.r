# -------------------------------------------------------------------------------
# Logistic Regression with a direct application of Jensens
# -------------------------------------------------------------------------------
set.seed(1)
n <- 350
f <- 500
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-1, gsize), rep(1, gsize), rep(0, p - 3 * gsize))
xb <- X %*% b
prob <- 1 / (1 + exp(-xb))
y <- rbinom(n, 1, prob)

lambda <- 1
a0 <- 1
b0 <- 200
w <- a0 / (a0 + b0)

f <- gsvb.logistic(y, X, groups, niter=50)
f <- gsvb.logistic(y, X, groups, niter=50, fit=f)

plot(f$m * f$g)
plot(f$g)

ff <- sparseGAM::SSGL(y, X, X, groups, family="binomial", lambda0=2, lambda1=1)
plot(ff$beta)

gsvb.logistic <- function(y, X, groups, niter=500, fit=NULL)
{
    n <- nrow(X) 
    p <- ncol(X)
    
    # initialize
    if (is.null(fit)) {
	m <- rnorm(p)
	s <- runif(p, 0.1, 0.2)
	g <- rep(0.1, p)
    } else {
	m <- f$m
	s <- f$s
	g <- f$g
    }
    
    S <- compute_S(X, m, s, g, groups)

    M <- matrix(ncol=0, nrow=p)
    Sig <- matrix(ncol=0, nrow=p)
    Gpr <- matrix(ncol=0, nrow=p)

    # main loop
    for (iter in 1:niter)
    {
	for (group in unique(groups))
	{
	    G <- which(groups == group)
	    Gc <- which(groups != group)

	    # rm G from S
	    S <- S / compute_S_G(X, m, s, g, G)
	    
	    m[G] <- optim(m[G], 
		fn=function(mG) opt_mu(mG, y, X, m, s, g, G, lambda, S), 
		control=list(maxit=20),
		method="BFGS")$par

	    s[G] <- optim(s[G], 
		fn=function(sG) opt_s(sG, y, X, m, s, g, G, lambda, S), 
		control=list(maxit=20),
		method="L-BFGS-B", lower=1e-3, upper=s[G][1] + 0.2)$par

	    g[G] <- opt_g(y, X, m, s, g, G, lambda, S)

	    # add G to S
	    S <- S * compute_S_G(X, m, s, g, G)

	    # cat("\n")
	    # print(m[G])
	    # print(s[G])
	    # print(g[G])
	    # cat(group, "\n")
	}

	M <- cbind(M, m)
	Sig <- cbind(Sig, s)
	Gpr <- cbind(Gpr, g)

	matplot(t(M), type="l")

	cat("\niter: ", iter)
    }
    
    return(list(m=m, s=s, g=g, M=M, Sig=Sig, Gpr=Gpr))

}


compute_S <- function(X, m, s, g, groups) 
{
    S <- rep(1, nrow(X))

    for (group in unique(groups)) {
	G <- which(groups == group)
	S <- S * compute_S_G(X, m, s, g, G)
    }
    return(S)
}


compute_S_G <- function(X, m, s, g, G)
{
    apply(X[ , G], 1, function(x) {
	(1 - g[G][1]) + g[G][1] * exp(sum(x * m[G] + 0.5 * x^2 * s[G]^2))
    })
}

n_mgf <- function(X, m, s)
{
    apply(X, 1, function(x) {
	exp(sum(x * m + 0.5 * x^2 * s^2))
    })
}


opt_mu <- function(m_G, y, X, m, s, g, G, lambda, S) 
{
    # maybe combine in a Monte Carlo step rather than use
    # Jensen's for this part?
    S <- S * n_mgf(X[ , G], m_G, s[G])

    sum(log1p(S) - y * (X[ , G] %*% m_G)) +
    lambda * sqrt(sum(s[G]^2) + sum(m_G^2))
}



opt_s <- function(s_G, y, X, m, s, g, G, lambda, S) 
{
    S <- S * n_mgf(X[ , G], m[g], s_G)

    sum(log1p(S)) -
    sum(log(s_G)) +
    lambda * sqrt(sum(s_G^2) + sum(m[G]^2))
}


opt_g <- function(y, X, m, s, g, G, lambda, S) 
{
    mk <- length(G)
    Ck <- mk * log(2) + (mk -1)/2 * log(pi) + lgamma( (mk + 1) / 2)
    S1 <- S * n_mgf(X[ , G], m[G], s[G])

    res <- 
	log(w / (1- w)) + 
	0.5 * mk - 
	Ck +
	mk * log(lambda) +
	0.5 * sum(log(2 * pi * s[G]^2)) -
	lambda * sqrt(sum(s[G]^2) + sum(m[G]^2)) +
	sum(y * X[ , G] %*% m[G]) -
	sum(log1p(S1)) +
	sum(log1p(S))

    sigmoid(res)
}


sigmoid <- function(x) 1.0 / (1.0 + exp(-x))

