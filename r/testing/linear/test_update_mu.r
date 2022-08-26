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

xtx <- crossprod(X, X)
yx <- t(X) %*% y

lambda <- 1
sigma <- 1
s  <- sapply(1:p, function(i) g1(i, xtx, lambda, sigma))

# ----------------------------------------
# Test func
# ----------------------------------------
gi <- 1
G <- which(groups == gi) 
Gc <- which(groups != gi)

g <- rep(1e-8, p)

f1(mu, s, g, xtx, yx, lambda, sigma, G, Gc)
update_mu(G - 1, Gc - 1, xtx, yx, mu, s, g, sigma, lambda)


# ----------------------------------------
# Timing
# ----------------------------------------
microbenchmark::microbenchmark(
    f1(mu, s, g, xtx, yx, lambda, sigma, G, Gc),
    update_mu(G - 1, Gc - 1, xtx, yx, mu, s, g, sigma, lambda),
    setup = {
	X <- matrix(rnorm(n * p), nrow=n, ncol=p)
	y <- X %*% b + rnorm(n, 0, 1)
	s  <- sapply(1:p, function(i) g1(i, xtx, lambda, sigma))
	xtx <- t(X) %*% X
	yx <- t(t(y) %*% X)
    }
)



# ----------------------------------------
# Third stage optimization
# ----------------------------------------
testx  <- replicate(1000, {
    f3(rep(0, 5), mu, s, g, xtx, yx, sigma, lambda, G, Gc, 1e3) -
    update_mu_3(rep(0, 5) , xtx, yx, mu, s, g, sigma, lambda, G-1, Gc-1, 1e3)
})

mean(testx)
plot(testx)

microbenchmark::microbenchmark(
    update_mu_3(mu[G], xtx, yx, mu, s, g, sigma, lambda, G-1, Gc-1, 4e3),
    f3(mu[G], mu, s, g, xtx, yx, sigma, lambda, G, Gc, 4e3),
    setup = {
	X <- matrix(rnorm(n * p), nrow=n, ncol=p)
	y <- X %*% b + rnorm(n, 0, 1)
	s  <- sapply(1:p, function(i) g1(i, xtx, lambda, sigma))
	xtx <- t(X) %*% X
	yx <- t(t(y) %*% X)
	G <- which(groups == 1)
	Gc <- which(groups != 1)
    }
)

