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

	    m[G] <- optim(m[G], 
		fn=function(mG) opt_mu(mG, y, X, m, s, ug, lambda, group,
			G, X.m, X.s),
		control=list(maxit=10),
		method="BFGS")$par

	    X.m[ , group] <- X[ , G] %*% m[G]

	    s[G] <- optim(s[G], 
		fn=function(sG) opt_s(sG, X, m, s, ug, lambda, group,
			G, X.m, X.s),
		control=list(maxit=10),
		method="L-BFGS-B", lower=1e-3, upper=s[G][1] + 0.2)$par

	    X.s[ , group] <- X[ , G]^2 %*% s[G]^2

	    g[G] <- ug[group] <- opt_g(y, X, m, s, ug, lambda, group, G,
		X.m, X.s)
	}

	M <- cbind(M, m)
	Sig <- cbind(Sig, s)
	Gpr <- cbind(Gpr, g)

	matplot(t(M), type="l")
	cat("\niter: ", iter)
    }
    
    # return(list(m=m, s=s, g=g, M=M, Sig=Sig, Gpr=Gpr))

# }


opt_mu <- function(m_G, y, X, m, s, ug, lambda, group, G,
    X.m, X.s, thresh=0.02, l=20) 
{
    ug[group] <- 1
    m[G] <- m_G

    xm <- X[ , G] %*% m[G]
    X.m[ , group] <- xm

    e_ll(X.m, X.s, ug, thresh, l) -
    sum(y * xm) +
    lambda * sqrt(sum(s[G]^2) + sum(m_G^2))
}


opt_s <- function(s_G, X, m, s, ug, lambda, group, G,
    X.m, X.s, thresh=0.02, l=20) 
{
    ug[group] <- 1
    s[G] <- s_G 

    xs <- sqrt(X[ , G]^2 %*% s[G]^2)
    X.s[ , group] <- xs

    e_ll(X.m, X.s, ug, thresh, l) -
    sum(log(s_G)) +
    lambda * sqrt(sum(s_G^2) + sum(m[G]^2))
}


opt_g <- function(y, X, m, s, ug, lambda, group, G, X.m, X.s,
	tau=0.02, l=20)
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

e_ll <- function(X.m, X.s, ug, tau, l=20)
{
    mid <- which(ug >= tau & ug <= (1-tau))
    big <- which(ug > (1-tau))
    gk <- ug[mid]
    tot <- 0

    if (length(mid) >= 10) print("shittles")
    
    if (length(mid) == 0) {
	if (length(big) == 1) {
	    mu <- X.m[ , big]
	    sig <- sqrt(X.s[ , big])
	} else {
	    mu <- apply(X.m[ , big], 1, sum)
	    sig <- sqrt(apply(X.s[ , big], 1, sum))
	}
	return(tll(mu, sig, l))
    }

    for (i in 0:(2^length(mid)-1)) 
    {
	sk <- as.integer(intToBits(i)[1:length(mid)])
	J <- c(mid[!!sk], big)
	
	if (length(J) == 1) {
	    mu <- X.m[ , J]
	    sig <- sqrt(X.s[ , J])
	} else {
	    mu <- apply(X.m[ , J], 1, sum)
	    sig <- sqrt(apply(X.s[ , J], 1, sum))
	}

	res <- tll(mu, sig, l)
	tot <- tot + prod(gk^sk * (1 - gk)^(1 - sk)) *  sum(res)
    }

    tot
}

Rcpp::sourceCpp("../new_bound_3/bound.cpp")

