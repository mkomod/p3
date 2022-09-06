# ------------------------------------------------------------------------------
# Jaakkola & Jordan Bound
# ------------------------------------------------------------------------------
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
a0 <- 1
b0 <- 200
w <- a0 / (a0 + b0)

f <- gsvb.logistic(y, X, groups, niter=50)
f <- gsvb.logistic(y, X, groups, niter=50, fit=f)

plot(f$m * f$g)
plot(f$g)

gsvb.logistic <- function(y, X, groups, niter=500, fit=NULL)
{
    # X <- cbind(1, X)
    # groups <- c(0, groups)

    n <- nrow(X) 
    p <- ncol(X)

    # initialize
    if (is.null(fit)) {
	m <- rnorm(p)
	s <- runif(p, 0.1, 0.2)
	g <- rep(0.05, p)
	l <- rgamma(n, 1, 1)
    } else {
	m <- f$m
	s <- f$s
	g <- f$g
	l <- f$l
    }


    M <- matrix(ncol=0, nrow=p)
    Sig <- matrix(ncol=0, nrow=p)
    Gpr <- matrix(ncol=0, nrow=p)
    L <- matrix(ncol=0, nrow=n)

    # main loop
    for (iter in 1:niter)
    {
	l <- opt_l(X, m, s, g)
	A <- diag(as.numeric(a(l)))
	XAX <- t(X) %*% A %*% X

	for (group in unique(groups))
	{
	    G <- which(groups == group)
	    Gc <- which(groups != group)

	    m[G] <- optim(m[G], 
		fn=function(mG) opt_mu(mG, y, X, XAX, m, s, g, G, Gc, lambda),
		control=list(maxit=20),
		method="BFGS")$par

	    s[G] <- optim(s[G], 
		fn=function(sG) opt_s(sG, XAX, m, G, lambda),
		control=list(maxit=20),
		method="L-BFGS-B", lower=1e-3, upper=s[G][1] + 0.2)$par

	    # g[G] <- opt_g(y, X, XAX, m, s, g, G, Gc, lambda)

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

    return(list(m=m, s=s, g=g, l=l, M=M, Sig=Sig, Gpr=Gpr))
} 


opt_mu <- function(m_G, y, X, XAX, m, s, g, G, Gc, lambda)
{
    res <- 0.5 * t(m_G) %*% XAX[G, G] %*% m_G +
    t(m_G) %*% XAX[G, Gc] %*% (g[Gc] * m[Gc]) +
    sum((0.5 - y) * X[ , G] %*% m_G) +
    lambda * (sum(s[G]^2 + m_G^2))^(1/2)

    print(res)
    res
}


opt_s <- function(s_G, XAX, m, G, lambda)
{
    res <- 0.5 * sum(diag(XAX[G, G]) * s_G^2) -
    sum(log(s_G)) +
    lambda * (sum(s_G^2 + m[G]^2))^(1/2)

    print(res)
    res
}

opt_g <- function(y, X, XAX, m, s, g, G, Gc, lambda) 
{
    mk <- length(G)
    Ck <- mk * log(2) + (mk -1)/2 * log(pi) + lgamma( (mk + 1) / 2)

    res <- 
	log(w / (1- w)) + 
	0.5 * mk + 
	0.5 * sum(log(2 * pi * s[G]^2)) +
	Ck +
	mk * log(lambda) -
	lambda * sqrt(sum(s[G]^2) + sum(m[G]^2)) +
	sum((y - 0.5) * X[ , G] %*% m[G]) -
	0.5 * t(m[G]) %*% XAX[G, G] %*% m[G] -
	0.5 * sum(diag(XAX[G, G]) * s[G]^2) -
	t(m[G]) %*% XAX[G, Gc] %*% (g[Gc] * m[Gc])

    return(sigmoid(res))
}

opt_l <- function(X, m, s, g) 
{
    sqrt((X %*% (g * m))^2 + apply(X, 1, function(x) x^2 %*% (g * s^2)))
}

a <- function(x) (sigmoid(x) - 0.5) / x

sigmoid <- function(x) 1/(1 + exp(-x))

