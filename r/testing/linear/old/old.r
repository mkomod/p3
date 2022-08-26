# ----------------------------------------
# Component-wise updates for mu
# ----------------------------------------

f1 <- function(m, i,  mu, s, g, xtx, yx, lambda, sigma, G, Gc) {
    mu[i] <- m
    0.5 * sigma^-2 * (sum(xtx[G, i] * mu[G] * m)) +
	sigma^-2 * sum(xtx[Gc, i] * g[Gc] * mu[Gc] * m) -
	sigma^-2 * m * yx[i] + lambda * sqrt(sum(s[G]^2 + mu[G]^2))
}


f2 <- function(i, mu, s, g, xtx, yx, lambda, sigma, G, Gc) {
    gni <- setdiff(G, i)
    if (length(gni) == 0) {
	res <- - sum(xtx[Gc, i]*mu[Gc]*g[Gc]) + yx[i]
    } else {
	res <- - 0.5 * sum(xtx[gni, i] * mu[gni]) - sum(xtx[Gc, i]*mu[Gc]*g[Gc]) + yx[i]
    }
    res / (xtx[i, i] + 2 * lambda * sigma^2)
}


# --- Implementation of Group sparse variational Bayes ---
source("funcs.r")


# ----------------------------------------
# Main function 
# ----------------------------------------
gsvb.fit <- function(y, X, groups, mu, s, g, a0, b0, lambda,
		     niter=1e4, tol=1e-4, verbose=T, sigma=1, level=2) 
{
    n <- nrow(X)
    p <- ncol(X)

    xtx <- t(X) %*% X
    yx <- t(t(y) %*% X)

    # initialize parameters
    s_old  <- s  <- sapply(1:p, function(i) g1(i, xtx, lambda, sigma))
    w <- a0 /(a0 + b0)

    current_level <- 1
    
    for (. in 1:niter) 
    {
	mu_old <- mu; s_old <- s; g_old <- g;

	for (gi in unique(groups))
	{
	    G <- which(groups == gi)
	    Gc <- which(groups != gi)
	   
	    if (current_level == 1)
		mu[G] <- f1(mu, s, g, xtx, yx, lambda, sigma, G, Gc)

	    if (current_level == 2)
		mu[G] <- optim(mu[G], fn = function(m) 
		    f2(m, mu, s, g, xtx, yx, lambda, sigma, G, Gc),
		    method="CG")$par

	    if (current_level == 3)
		mu[G] <- optim(mu[G], fn = function(m) 
		    f3(m, mu, s, g, xtx, yx, lambda, sigma, G, Gc, 5e2),
		    method="CG", control=list(abstol=1e-3))$par

	    for (i in G) 
		s[i]  <- optim(s[i], fn = function(si) 
		    g2(si, i, mu, s, g, xtx, yx, lambda, sigma, G, Gc),
		    method="Brent", lower=1e-4, upper=s[i] * 4)$par
	    
	    # g[G] <- h1(mu, s, g, xtx, yx, lambda, sigma, G, Gc, w)
	    g[G] <- update_g(G-1, Gc-1, xtx, yx, mu, s, g, sigma, lambda, w)
	    
	}
	
	if (verbose)
	    cat(.)

	if ((sum(abs(s - s_old)) < tol) &&
	    (sum(abs(mu - mu_old)) < tol) &&
	    sum(abs(g - g_old)) < tol)
	    if (level == current_level) {
		break;
	    } else {
		current_level <- current_level + 1
	    }
    }
    return(list(m=mu, s=s, g=g))
}


# ----------------------------------------
# Optimizers for mu / sigma / gamma
# ----------------------------------------

# ---------------- mu ----------------
f1 <- function(mu, s, g, xtx, yx, lambda, sigma, G, Gc) {
    psi <- solve(xtx[G, G] + 2 * lambda * sigma^2 * diag(length(G)))
    res <- psi %*% (yx[G] - xtx[G, Gc] %*% (g[Gc] *  mu[Gc]))
    return(res)
}


f2 <- function(m, mu, s, g, xtx, yx, lambda, sigma, G, Gc) {
    xtxm <- crossprod(xtx[G, G], m)
    res <- 0.5 * sigma^-2 * crossprod(xtxm, m) +
	sigma^-2 * crossprod(crossprod(xtx[G, Gc], m), g[Gc] * mu[Gc]) -
	sigma^-2 * crossprod(yx[G], m) +
	lambda * sum(s[G]^2 + m^2)^0.5
    return(res)
}


f3 <- function(m, mu, s, g, xtx, yx, lambda, sigma, G, Gc, mcn=5e2) {
    xtxm <- crossprod(xtx[G, G], m)
    # mci <- norm(matrix(rnorm(mcn * length(G), mean=m, sd=s[G]), nrow=length(m)),
	# type="F") / sqrt(mcn)
    mci <- mean(replicate(mcn, l2(rnorm(length(G), mean=m, sd=s[G]))))
    
    res <- 0.5 * sigma^-2 * crossprod(xtxm, m) +
	sigma^-2 * crossprod(crossprod(xtx[G, Gc], m), g[Gc] * mu[Gc]) -
	sigma^-2 * crossprod(yx[G], m) +
	lambda * mci
    return(res)
}

# --------------- sigma ---------------
g1 <- function(i, xtx, lambda, sigma) {
    1/sqrt(xtx[i, i] / sigma^2 + 2 * lambda)
}


g2 <- function(si, i, mu, s, g, xtx, yx, lambda, sigma, G, Gc) {
    s[i] <- si
    0.5 * sigma^-2 * xtx[i, i] * si^2 - log(si) + lambda * sqrt(sum(s[G]^2 + mu[G]^2))
}



# --------------- gamma ---------------
h1 <- function(mu, s, g, xtx, yx, lambda, sigma, G, Gc, w=0.5) 
{
    mk <- length(G)
    con <- 0.5 * sigma^-2
    res <- log(w/(1-w)) + mk/2 + sigma^-2 * t(yx[G]) %*% mu[G] +
	mk * log(lambda) - 
    (
	2 * con * sum(t(xtx[Gc, G] * g[Gc] * mu[Gc]) * mu[G]) +		# [t]
	con * sum(diag(xtx)[G] * (s[G]^2)) +				# [t]
	con * sum(xtx[G, G] * mu[G] %*% t(mu[G])) -			# [t]
	0.5 * sum(log(2 * pi * s[G]^2)) +
	mk * log(2) + 
	((mk - 1)/2) * log(pi) +
     	lgamma((mk+1)/2) +
	lambda * sqrt(sum(s[G]^2) + sum(mu[G]^2))		# [t]
    )
    return(sigmoid(res))
}


#
#
#
n <- 100
p <- 1000
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-4, gsize), rep(8, gsize), rep(0, p - 3 * gsize))
y <- X %*% b + rnorm(n, 0, 1)


# ----------------------------------------
# Comparison of C++ version to R version
# ----------------------------------------
# init params
mu <- runif(p, -0.02, 0.2)
g <- rep(0.5, p)

m1 <-      fit(y, X, groups, mu, rep(0, p), g, 1, 200, .5, 200, 5e-2, T, 1)
m2 <- gsvb.fit(y, X, groups, mu, rep(0, p), g, 1, 200, .5, 200, 5e-2, T, 1, 2)


plot(m1$m)
plot(m1$m * m1$g)
plot(m2$m)
plot(m2$m * m2$g)



a <- matrix(1:6, nrow=2)
aa <- c(0.3, 0.9, 0.7)
a
aa * a






# -------------------------------------------------- 
# Old optimizers
# -------------------------------------------------- 
# ---------------- mu ----------------
f1 <- function(mu, s, g, xtx, yx, lambda, sigma, G, Gc) {
    psi <- solve(xtx[G, G] + 2 * lambda * sigma^2 * diag(length(G)))
    res <- psi %*% yx[G] - psi %*% xtx[G, Gc] %*% (g[Gc] * mu[Gc])
    return(res)
}


f2 <- function(m, mu, s, g, xtx, yx, lambda, sigma, G, Gc) {
    xtxm <- crossprod(xtx[G, G], m)
    res <- 0.5 * sigma^-2 * crossprod(xtxm, m) +
	sigma^-2 * crossprod(crossprod(xtx[G, Gc], m), g[Gc] * mu[Gc]) -
	sigma^-2 * crossprod(yx[G], m) +
	lambda * sum(s[G]^2 + m^2)^0.5
    return(res)
}


f3 <- function(m, mu, s, g, xtx, yx, lambda, sigma, G, Gc, mcn=5e2) {
    xtxm <- crossprod(xtx[G, G], m)
    # mci <- norm(matrix(rnorm(mcn * length(G), mean=m, sd=s[G]), nrow=length(m)),
	# type="F") / sqrt(mcn)
    mci <- mean(replicate(mcn, l2(rnorm(length(G), mean=m, sd=s[G]))))
    
    res <- 0.5 * sigma^-2 * crossprod(xtxm, m) +
	sigma^-2 * crossprod(crossprod(xtx[G, Gc], m), g[Gc] * mu[Gc]) -
	sigma^-2 * crossprod(yx[G], m) +
	lambda * mci
    return(res)
}

# --------------- sigma ---------------
g1 <- function(i, xtx, lambda, sigma) {
    1/sqrt(xtx[i, i] / sigma^2 + 2 * lambda)
}


g2 <- function(si, i, mu, s, g, xtx, yx, lambda, sigma, G, Gc) {
    s[i] <- si
    0.5 * sigma^-2 * xtx[i, i] * si^2 - log(si) + lambda * sqrt(sum(s[G]^2 + mu[G]^2))
}



# --------------- gamma ---------------
h1 <- function(mu, s, g, xtx, yx, lambda, sigma, G, Gc, w=0.5) 
{
    mk <- length(G)
    con <- 0.5 * sigma^-2
    res <- log(w/(1-w)) + mk/2 + sigma^-2 * t(yx[G]) %*% mu[G] +
	mk * log(lambda) - 
    (
	2 * con * sum(t(xtx[Gc, G] * g[Gc] * mu[Gc]) * mu[G]) +		# [t]
	con * sum(diag(xtx)[G] * (s[G]^2)) +				# [t]
	con * sum(xtx[G, G] * mu[G] %*% t(mu[G])) -			# [t]
	0.5 * sum(log(2 * pi * s[G]^2)) +
	mk * log(2) + 
	((mk - 1)/2) * log(pi) +
     	lgamma((mk+1)/2) +
	lambda * sqrt(sum(s[G]^2) + sum(mu[G]^2))		# [t]
    )
    return(sigmoid(res))
}


# ---------------- elbo ----------------
elbo.mci <- function(f, y, X, groups, mcn=1e3)
{
    n <- nrow(X)
    p <- ncol(X)
    xtx <- crossprod(X, X)

    attach(f); attach(parameters)
    w <- a0 / (a0 + b0)

    # non stochasitc part of the ELBO
    res <- - 0.5 * n * log(2 * pi * sigma^2) -
	0.5 * sigma^-2 * crossprod(y, y)
    
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
	    if (runif(1) < g[gi]) {
		a <- a + log(1-g[gi]+1e-8) - log(1-w +1e-8)
	    } else {
		G <- which(groups[gi] == groups)
		b <- rnorm(length(G), mu[G], s[G])
		a <- a + log(w) + lmvde(b, lambda, length(G)) - log(g[gi]) - sum(dnorm(b, mu[G], s[G], log=T))
	    }
	}
	as.vector(r) + a
    })

    detach(f); detach(parameters) 

    res + mean(mci, na.rm=T)
}

