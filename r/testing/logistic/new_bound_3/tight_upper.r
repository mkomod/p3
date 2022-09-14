set.seed(1)
n <- 350
p <- 500
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-1, gsize), rep(1, gsize), 
       rep(0, p - 3 * gsize))
xb <- X %*% b
prob <- 1 / (1 + exp(-xb))
y <- rbinom(n, 1, prob)

lambda <- 1
tau = 1 - 1e-4
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
	g <- rep(1, p)
	g <- 0 + !!b
    } else {
	m <- f$m
	s <- f$s
	g <- f$g
    }

    ug <- g[!duplicated(groups)]

    M <- matrix(ncol=0, nrow=p)
    Sig <- matrix(ncol=0, nrow=p)
    Gpr <- matrix(ncol=0, nrow=p)

    X.m <- sapply(unique(groups), function(group) {
	G <- which(groups == group)
	X[ , G] %*% m[G]
    })

    X.s <- sapply(unique(groups), function(group) {
	G <- which(groups == group)
	X[ , G]^2 %*% s[G]^2
    })

    # main loop
    for (iter in 1:niter)
    {
	for (group in unique(groups))
	{
	    cat(group)

	    G <- which(groups == group)
	    Gc <- which(groups != group)

	    # m[G] <- optim(m[G], 
		# fn=function(mG) opt_mu(mG, y, X, m, s, ug, lambda, group,
			# G, X.m, X.s),
		# control=list(maxit=10),
		# method="BFGS")$par
	    m[G] <- opt_m_cpp(y, X, m, s, ug, lambda, group-1, G-1, X.m, X.s, 0.02, 20)
	    X.m[ , group] <- X[ , G] %*% m[G]
	    
	    plot(m)
	    # s[G] <- optim(s[G], 
		# fn=function(sG) opt_s(sG, X, m, s, ug, lambda, group,
			# G, X.m, X.s),
		# control=list(maxit=10),
		# method="L-BFGS-B", lower=1e-3, upper=s[G][1] + 0.2)$par

	    s[G] <- opt_s_cpp(y, X, m, s, ug, lambda, group-1, G-1, X.m, X.s, 0.02, 20)
	    X.s[ , group] <- X[ , G]^2 %*% s[G]^2

	    # g[G] <- ug[group] <- opt_g(y, X, m, s, ug, lambda, group, G,
		# X.m, X.s)

	    g[G] <- ug[group] <- opt_g_cpp(y, X, m, s, ug, lambda, group-1,
		G-1, X.m, X.s, 0.02, 20, w)
	}

	M <- cbind(M, m)
	Sig <- cbind(Sig, s)
	Gpr <- cbind(Gpr, g)

	# matplot(t(M), type="l")
	cat("\niter: ", iter)
    }
    
    # return(list(m=m, s=s, g=g, M=M, Sig=Sig, Gpr=Gpr))

# }

opt_mu <- function(m_G, y, X, m, s, ug, lambda, group, G,
    X.m, X.s, thresh=0.02, l=50) 
{
    ug[group] <- 1
    m[G] <- m_G

    xm <- X[ , G] %*% m[G]
    X.m[ , group] <- xm

    # print(e_ll(X.m, X.s, ug, thresh, l))
    # print(sum(y * xm))
    # print(lambda * sqrt(sum(s[G]^2) + sum(m_G^2)))

    ell(X.m, X.s, ug, thresh, l) -
    sum(y * xm) +
    lambda * sqrt(sum(s[G]^2) + sum(m_G^2))
}


opt_s <- function(s_G, X, m, s, ug, lambda, group, G,
    X.m, X.s, thresh=0.02, l=50) 
{
    ug[group] <- 1
    s[G] <- s_G 

    xs <- X[ , G]^2 %*% s[G]^2
    X.s[ , group] <- xs

    ell(X.m, X.s, ug, thresh, l) -
    sum(log(s_G)) +
    lambda * sqrt(sum(s_G^2) + sum(m[G]^2))
}


opt_g <- function(y, X, m, s, ug, lambda, group, G, X.m, X.s,
	tau=0.02, l=50)
{
    mk <- length(G)
    Ck <- mk * log(2) + (mk -1)/2 * log(pi) + lgamma( (mk + 1) / 2)

    ug[group] <- 1
    S1 <- e_ll(X.m, X.s, ug, tau, l)

    ug[group] <- 0
    S0 <- e_ll(X.m, X.s, ug, tau, l)

    res <- 
	log(w / (1- w)) + 
	0.5 * mk - 
	Ck +
	mk * log(lambda) +
	0.5 * sum(log(2 * pi * s[G]^2)) -
	lambda * sqrt(sum(s[G]^2) + sum(m[G]^2)) -
	S1 + S0 + sum(y * (X[ , G] %*% m[G]))

    sigmoid(res)
}


Rcpp::sourceCpp("../new_bound_3/bound.cpp")

