
Rcpp::sourceCpp("../src/fit.cpp")


# ----------------------------------------
# Generate data
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

for (i in G) {
    print(optim(s[i], fn=
	function(si) g2(si, i, mu, s, g, xtx, yx, lambda, sigma, G, Gc),
    method="L-BFGS-B", lower=1e-3)$par)
    print(update_s2(xtx, mu, s, sigma, lambda, i - 1, G - 1)[1])
}

print(update_s(xtx, mu, s, sigma, lambda, G - 1))

g2 <- function(si, i, mu, s, g, xtx, yx, lambda, sigma, G, Gc) {
    s[i] <- si
    0.5 * sigma^-2 * xtx[i, i] * si^2 - log(si) + lambda * sqrt(sum(s[G]^2 + mu[G]^2))
}

