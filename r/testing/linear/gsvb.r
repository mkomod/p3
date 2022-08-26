# --- Implementation of Group sparse variational Bayes ---
source("funcs.r")
Rcpp::sourceCpp("../src/fit.cpp")


# ----------------------------------------
# Main function 
# ----------------------------------------
# Wrapper code that uses the C++ functions directly
#
gsvb.fit <- function(y, X, groups, lambda, a0=1, b0=length(groups), 
	sigma=1, mu=runif(ncol(X), -0.2, 0.2), g=rep(0.5, ncol(X)), 
	niter=1e4, tol=1e-4, verbose=T)
{
    n <- nrow(X)
    p <- ncol(X)

    xtx <- t(X) %*% X
    yx <- t(t(y) %*% X)
    w <- a0 /(a0 + b0)

    # initialize s
    s <- sapply(1:p, function(i) 1/sqrt(xtx[i, i] / sigma^2 + 2 * lambda))

    for (. in 1:niter) 
    {
	mu_old <- mu; s_old <- s; g_old <- g;

	for (gi in unique(groups))
	{
	    G <- which(groups == gi)
	    Gc <- which(groups != gi)
	   
	    mu[G] <- update_mu(G-1, Gc-1, xtx, yx, mu, s, g, sigma, lambda)
	    s[G]  <- update_s(G-1, xtx, mu, s, sigma, lambda)
	    g[G]  <- update_g(G-1, Gc-1, xtx, yx, mu, s, g, sigma, lambda, w)
	}
	
	if (verbose)
	    cat(.)

	if ((sum(abs(s - s_old)) < tol) &&
	    (sum(abs(mu - mu_old)) < tol) &&
	    sum(abs(g - g_old)) < tol)
		break;
    }
    return(list(m=mu, s=s, g=g))
}



