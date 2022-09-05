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

lambda <- 1
tau = 0.5
a0 <- 1
b0 <- 200
w <- a0 / (a0 + b0)

# gsvb.logistic <- function(y, X, groups, niter=500, fit=NULL)
# {
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
    
    M <- matrix(ncol=0, nrow=p)
    Sig <- matrix(ncol=0, nrow=p)
    # Gpr <- matrix(ncol=0, nrow=p)

    # main loop
    for (iter in 1:niter)
    {
	for (group in unique(groups))
	{
	    G <- which(groups == group)
	    Gc <- which(groups != group)

	    m[G] <- optim(m[G], 
		fn=function(mG) opt_mu(mG, y, X, m, s, g, G, lambda, tau), 
		control=list(maxit=20),
		method="BFGS")$par

	    s[G] <- optim(s[G], 
		fn=function(sG) opt_s(sG, y, X, m, s, g, G, lambda, tau), 
		control=list(maxit=20),
		method="L-BFGS-B", lower=1e-3, upper=s[G][1] + 0.2)$par

	    # g[G] <- opt_g(y, X, m, s, g, G, lambda, S)
	}

	M <- cbind(M, m)
	Sig <- cbind(Sig, s)
	# Gpr <- cbind(Gpr, g)

	matplot(t(M), type="l")

	cat("\niter: ", iter)
    }
    
    # return(list(m=m, s=s, g=g, M=M, Sig=Sig, Gpr=Gpr))

# }


nb3 <- function(mu, sig) {
    sig / sqrt(2 * pi) * exp(-mu^2 / (2 * sig^2)) + mu * (1 - pnorm(-mu/sig)) +
    exp(mu + 0.5*sig^2) * pnorm(-mu/sig - sig) +
    exp(-mu + 0.5*sig^2) * pnorm(mu/sig - sig)
}


opt_mu <- function(m_G, y, X, m, s, g, G, lambda, tau) 
{
    m[G] <- m_G 
    J <- union(G, which(g > tau))

    sum(
	nb3(X[ , J] %*% m[J], X[ , J]^2 %*% s[J]^2) -
	y * (X[ , G] %*% m[G])
    ) +
    lambda * sqrt(sum(s[G]^2) + sum(m_G^2))
}



opt_s <- function(s_G, y, X, m, s, g, G, lambda, tau) 
{
    s[G] <- s_G 
    J <- union(G, which(g > tau))

    sum(
	nb3(X[ , J] %*% m[J], X[ , J]^2 %*% s[J]^2)
    ) -
    sum(log(s_G)) +
    lambda * sqrt(sum(s_G^2) + sum(m[G]^2))
}


# opt_g <- function(y, X, m, s, g, G, lambda, S) 
# {
#     mk <- length(G)
#     Ck <- mk * log(2) + (mk -1)/2 * log(pi) + lgamma( (mk + 1) / 2)
#     S1 <- S * n_mgf(X[ , G], m[G], s[G])

#     res <- 
# 	log(w / (1- w)) + 
# 	0.5 * mk + 
# 	Ck +
# 	mk * log(lambda) +
# 	0.5 * sum(log(2 * pi * s[G]^2)) -
# 	lambda * sqrt(sum(s[G]^2) + sum(m[G]^2)) +
# 	sum(y * X[ , G] %*% m[G]) -
# 	sum(log1p(S1)) +
# 	sum(log1p(S))

#     sigmoid(res)
# }
