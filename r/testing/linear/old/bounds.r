
Rcpp::sourceCpp(code = '
#include "Rcpp.h"

// [[Rcpp::export]]
double l2(Rcpp::NumericVector a) {
    return sqrt(sum(pow(a, 2.0)));
}
')

Rcpp::sourceCpp(code = '
#include "Rcpp.h"

// [[Rcpp::export]]
double sigmoid(double x) {
    return 1.0 / (1.0 + exp(-x));
}
')


f0 <- function(m, i,  mu, s, g, xtx, yx, lambda, G, Gc) {
    mu[i] <- m
    mcn <- 500
    p <- length(G)

    mci <- mean(apply(matrix(rnorm(p * mcn, mu[G], s[G]), byrow=T, ncol=p), 1, l2))
    0.5 * (sum(xtx[G, i] * mu[G] * m)) +
	sum(xtx[Gc, i] * g[Gc] * mu[Gc] * m) -
	m * yx[i] + lambda * mci
}


f1 <- function(m, i,  mu, s, g, xtx, yx, lambda, G, Gc) {
    mu[i] <- m
    0.5 * (sum(xtx[G, i] * mu[G] * m)) +
	sum(xtx[Gc, i] * g[Gc] * mu[Gc] * m) -
	m * yx[i] + lambda * sqrt(sum(s[G]^2 + mu[G]^2))
}


f2 <- function(m, i, mu, s, g, xtx, yx, lambda, G, Gc) {
    mu[i] <- m
    0.5 * (sum(xtx[G, i] * mu[G] * m)) +
	sum(xtx[Gc, i] * g[Gc] * mu[Gc] * m) -
	m * yx[i] + lambda * (1 + sum(s[G]^2 + mu[G]^2))
}


f3 <- function(m, i, mu, s, g, xtx, yx, lambda, G, Gc) {
    mu[i] <- m
    0.5 * (sum(xtx[G, i] * mu[G] * m)) +
	sum(xtx[Gc, i] * g[Gc] * mu[Gc] * m) -
	m * yx[i] + lambda * (sum(abs(s[G]) + abs(mu[G])))
}


n <- 100
p <- 100
X <- matrix(rnorm(n * p), ncol=p)
y <- rnorm(n)

xtx <- t(X) %*% X
yx <- t(y %*% X)

mu <- rnorm(p)
s <- rgamma(p, 1, 1)
g <- rep(0.5, p)
lambda <- 0.5

G <- c(1:15)
Gc <- setdiff(1:p, G)

# ms <- seq(1.2, 1.25, by=0.001)
ms <- seq(-2, 2, by=0.001)
f0s <- sapply(ms, function(m) f0(m, 1, mu, s, g, xtx, yx, lambda, G, Gc))
f1s <- sapply(ms, function(m) f1(m, 1, mu, s, g, xtx, yx, lambda, G, Gc))
f2s <- sapply(ms, function(m) f2(m, 1, mu, s, g, xtx, yx, lambda, G, Gc))
f3s <- sapply(ms, function(m) f3(m, 1, mu, s, g, xtx, yx, lambda, G, Gc))

plot(ms, f0s, type="l", col=1, lwd=5, ylim=c(min(f1s), max(f2s)))
lines(ms, f1s, lwd=4, 	col=2, )
lines(ms, f2s, lwd=4, 	col=3, )
lines(ms, f3s, lwd=4, 	col=4, )

abline(v=ms[which.min(f0s)], col=1)
abline(v=ms[which.min(f1s)], col=2)
abline(v=ms[which.min(f2s)], col=3)
abline(v=ms[which.min(f3s)], col=4)

ms[which.min(f0s)]
ms[which.min(f1s)]
ms[which.min(f2s)]
ms[which.min(f3s)]

f2.obj <- function(m, i, mu, s, g, xtx, yx, lambda, G, Gc) {
    mu[i] <- m
    0.5 * (sum(xtx[G, i] * mu[G] * m)) +
	sum(xtx[Gc, i] * g[Gc] * mu[Gc] * m) -
	m * yx[i] + lambda * m^2
}

f2.minimizer <- function(i, mu, s, g, xtx, yx, lambda, G, Gc) {
    gni <- setdiff(G, i)
    res <- - 0.5 * sum(xtx[gni, i] * mu[gni]) - sum(xtx[Gc, i]*mu[Gc]*g[Gc]) + yx[i]
    res / (xtx[i, i] + 2 * lambda)
}

optim(0, function(x) f2.obj(x, i=1, mu=mu, s=s, g=g, xtx=xtx, yx=yx, lambda=lambda, G=G, Gc=Gc))
f2.minimizer(1, mu, s, g, xtx, yx, lambda, G, Gc) 


# sigma updates
g0 <- function(si, i, mu, s, g, xtx, yx, lambda, G, Gc) {
    s[i] <- si
    mcn <- 500
    p <- length(G)

    mci <- mean(apply(matrix(rnorm(p * mcn, mu[G], s[G]), byrow=T, ncol=p), 1, l2))
    0.5 * xtx[i, i] * si^2 - log(si) + lambda * mci
}

g1 <- function(si, i, mu, s, g, xtx, yx, lambda, G, Gc) {
    s[i] <- si
    0.5 * xtx[i, i] * si^2 - log(si) + lambda * sqrt(sum(s[G]^2 + mu[G]^2))
}

g2 <- function(si, i, mu, s, g, xtx, yx, lambda, G, Gc) {
    s[i] <- si
    0.5 * xtx[i, i] * si^2 - log(si) + lambda * (1 + sum(s[G]^2 + mu[G]^2))
}

ss <- seq(0.01, 2, by=0.001)
g0s <- sapply(ss, function(si) g0(si, 1, mu, s, g, xtx, yx, lambda, G, Gc))
g1s <- sapply(ss, function(si) g1(si, 1, mu, s, g, xtx, yx, lambda, G, Gc))
g2s <- sapply(ss, function(si) g2(si, 1, mu, s, g, xtx, yx, lambda, G, Gc))

plot(ss, g0s, type="l", col=1, lwd=5)
lines(ss, g1s, lwd=4, 	col=2, )
lines(ss, g2s, lwd=4, 	col=3, )

abline(v=ss[which.min(g0s)], col=1)
abline(v=ss[which.min(g1s)], col=1)
abline(v=ss[which.min(g2s)], col=1)

ss[which.min(g0s)]
ss[which.min(g1s)]
ss[which.min(g2s)]

g2.minimizer <- function(i, xtx, lambda) {
    1/sqrt(xtx[i, i] + 2 * lambda)
}

optim(0.05, function(si) g2(si, 1, mu, s, g, xtx, yx, lambda, G, Gc), control=list(reltol=1e-14))
g2.minimizer(1, xtx, lambda)


# gamma minimizer
h.minimizer <- function(mu, s, g, xtx, yx, lambda, G, Gc, w=0.5) {
    mk <- length(G)
    res <- log(w/(1-w)) - (
	sum(t(xtx[Gc, G] * g[Gc] * mu[Gc]) * mu[G]) +
	0.5 * sum(diag(xtx)[G] * (s[G]^2)) +
	0.5 * sum(xtx[G, G] * mu[G] %*% t(mu[G])) -
	t(yx[G]) %*% mu[G] +
	lambda * sqrt(sum(s[G]^2 + mu[G]^2)) - mk/2
    )
    return(sigmoid(res))
}

h.minimizer(mu, s, g, xtx, yx, lambda, G, Gc)



