library(Rcpp)

source("./gsvb.r")
Rcpp::sourceCpp("../src/fit.cpp")

# ----------------------------------------
# Test data
# ----------------------------------------
n <- 100
p <- 1000
gsize <- 5
groups <- c(rep(1:(p/gsize), each=gsize))

X <- matrix(rnorm(n * p), nrow=n, ncol=p)
b <- c(rep(0, gsize), rep(-4, gsize), rep(8, gsize), rep(0, p - 3 * gsize))
y <- X %*% b + rnorm(n, 0, 1)

lambda <- 1
sigma <- 1
s  <- sapply(1:p, function(i) g1(i, xtx, lambda, sigma))

# ----------------------------------------
# Test func
# ----------------------------------------
xtx <- t(X) %*% X
yx <- t(t(y) %*% X)

gi <- 1
G <- which(groups == gi) 
Gc <- which(groups != gi)

h1(mu, s, g, xtx, yx, lambda, sigma, G, Gc, 0.5)
update_g(G - 1, Gc - 1, xtx, yx, mu, s, g, sigma, lambda, 0.5)


# ----------------------------------------
# Timing
# ----------------------------------------
microbenchmark::microbenchmark(
    h1(mu, s, g, xtx, yx, lambda, sigma, G, Gc, 0.5),
    update_g(G - 1, Gc - 1, xtx, yx, mu, s, g, sigma, lambda, 0.5),
    setup = {
	X <- matrix(rnorm(n * p), nrow=n, ncol=p)
	y <- X %*% b + rnorm(n, 0, 1)
	s  <- sapply(1:p, function(i) g1(i, xtx, lambda, sigma))
	xtx <- t(X) %*% X
	yx <- t(t(y) %*% X)
    }
)

