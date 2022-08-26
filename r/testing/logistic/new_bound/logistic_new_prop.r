# Impelementation GSVB for logistic regression

# ----------------------------------------
# Generate test data
# ----------------------------------------
set.seed(1)
n <- 300
p <- 150
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


f <- gsvb.logistic(y, X, groups, niter=20)

# ----------------------------------------
# Logistic GSVB
# ----------------------------------------
gsvb.logistic <- function(y, X, groups, niter=500)
{
    n <- nrow(X) 
    p <- ncol(X)
    
    # initialize
    m <- rnorm(p)
    s <- runif(p, 0.1, 0.2)
    g <- rep(0.5, p)
    
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
	    S <- S - compute_S_G(X, m, s, g, G, Gc)
	    
	    # opt mu
	    m[G] <- optim(m[G], 
		fn=function(mG) opt_mu(mG, y, X, m, s, g, G, Gc, lambda, S), 
		control=list(maxit=20),
		method="BFGS")$par

	    # opt s
	    s[G] <- optim(s[G], 
		fn=function(sG) opt_s(sG, y, X, m, s, g, G, Gc, lambda, S), 
		control=list(maxit=20),
		method="L-BFGS-B", lower=1e-3, upper=5)$par

	    g[G] <- opt_g(y, X, m, s, g, G, Gc, lambda, S)

	    cat("\n")

	    print(m[G])
	    print(s[G])
	    print(g[G])

	    # add G to S
	    S <- S + compute_S_G(X, m, s, g, G, Gc)
	    cat(group, "\n")
	}

	M <- cbind(M, m)
	Sig <- cbind(Sig, s)
	Gpr <- cbind(Gpr, g)

	matplot(t(M), type="l")

	cat("\n\niter: ", iter, "\n\n\n")
    }
    
    return(list(m=m, s=s, g=g, M=M, Sig=Sig, Gpr=Gpr))

}


compute_S <- function(X, m, s, g, groups) 
{
    gg <- g[!duplicated(groups)]
    tmp <- matrix(1, nrow=length(gg), ncol=length(gg))
    diag(tmp) <- ifelse(gg, gg, 1)
    G <- outer(gg, gg) / tmp

    E_bi_bj <- 
	(outer(m , m) + diag(s^2)) * G[groups, groups]

    apply(X, 1, function(x) {
	sum(outer(x, x) * E_bi_bj)
    })
}

compute_S_G <- function(X, m, s, g, G, Gc)
{
    A <- outer(m[G], m[Gc]) * outer(g[G], g[Gc])
    B <- (outer(m[G], m[G]) + diag(s[G]^2)) * g[G][1] 
    apply(X, 1, function(x) {
	2 * sum(outer(x[G], x[Gc]) * A) + sum(outer(x[G], x[G]) * B)
    })
}


compute_S_G_K1 <- function(X, m, s, g, G, Gc)
{
    A <- outer(m[G], m[Gc]) * outer(rep(1, length(G)), g[Gc])
    B <- (outer(m[G], m[G]) + diag(s[G]^2))
    apply(X, 1, function(x) {
	2 * sum(outer(x[G], x[Gc]) * A) + sum(outer(x[G], x[G]) * B)
    })
}



opt_mu <- function(m_G, y, X, m, s, g, G, Gc, lambda, S) 
{
    m[G] <- m_G
    S <- S + compute_S_G_K1(X, m, s, g, G, Gc)

    sum((0.5 - y) * (X[ , G] %*% m_G) + 0.5 * (2 + S)^0.5) +
    lambda * sqrt(sum(s[G]^2) + sum(m_G^2))
}


opt_s <- function(s_G, y, X, m, s, g, G, Gc, lambda, S) 
{
    s[G] <- s_G
    S <- S + compute_S_G_K1(X, m, s, g, G, Gc)

    sum(0.5 * sqrt(2 + S)) -
    sum(log(s[G])) +
    lambda * sqrt(sum(s[G]^2) + sum(m[G]^2))
}

opt_g <- function(y, X, m, s, g, G, Gc, lambda, S) 
{
    mk <- length(G)
    Ck <- mk * log(2) + (mk -1) * log(pi) + lgamma( (mk + 1) / 2)
    S1 <- S + compute_S_G_K1(X, m, s, g, G, Gc)

    res <- log(w / (1- w)) + 0.5 * mk + 
	0.5 * sum(log(2 * pi * s[G]^2)) -
	Ck -
	lambda * sqrt(sum(s[G]^2) + sum(m[G]^2)) +
	mk * log(lambda) -
	sum(
	    (0.5 - y) * X[ , G] %*% m[G] +
	    0.5 * sqrt(2 + S1) -
	    0.5 * sqrt(2 + S)
	)
    sigmoid(res)
}

sigmoid <- function(x) 1.0 / (1.0 + exp(-x))

