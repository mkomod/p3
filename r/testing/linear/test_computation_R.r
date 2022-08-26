library(Rcpp)
Rcpp::sourceCpp("../src/linear_u.cpp")

# ----------------------------------------
# R implementation
# ----------------------------------------
r_compute_R <- function(yty, yx, xtx, groups, mu, S, g, p)
{
    res <- 0
    cij <- 0
    cijn <- 0

    for (i in 1:p) {
	mindx <- min(which(groups == groups[i])) - 1
	for (j in 1:p) {
	    if (groups[i] == groups[j]) {
		grp <- groups[[i]]
		cij <- cij + xtx[i, j] * g[i] * (S[[grp]][i - mindx, j - mindx] + mu[i]*mu[j])
	    } else {
		cijn <- cijn + xtx[i, j] * g[i] * g[j] * mu[i] * mu[j]
	    }
	}
    }
    print(cij)
    print(cijn)
    res <- cij + cijn
    print(res)
    res <- yty - 2 * sum(yx * mu * g) + res
    return(res)
}


# ----------------------------------------
# Compare against C++ imp
# ----------------------------------------
yty <- sum(y * y)
yx <- t(t(y) %*% X)
xtx <- t(X) %*% X

# generate random mu
set.seed(1)
mu <- rnorm(p)
s <- rgamma(p, 1, 1)
g <- rep(runif(p/5), each=5)


working <- approxEq(
    r_compute_S(yty, yx, xtx, groups, mu, s, g, p),
    compute_S(yty, yx, xtx, groups, mu, s, g, p)
)


Rcpp::sourceCpp("../src/linear_u.cpp")
mu <- rnorm(p, 0, 0.2)
f <- fit_linear_u(y, X, groups-1, 1, 1, 200, 1e-3, 1e-3, mu,
	rep(0.2, p), rep(0.5, p), T, 1, 500, 60, 1e-3, T)

S <- lapply(f$s, function(s) matrix(s, nrow=5))
r_compute_R(yty, as.vector(yx), xtx, groups, f$m, S, f$g, 1000)
solve(0.02 * xtx[1:5, 1:5] + diag(1:5))


ss <- expm::sqrtm(S[[25]])
xmvt <- replicate(10000, (ss %*% rnorm(5) + 1:5)[1])
xmv <- mvtnorm::rmvnorm(10000, 1:5, S[[25]])
hist(xmvt)
hist(xmv[ , 1])




