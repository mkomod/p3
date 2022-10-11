b <- matrix(rnorm(10), ncol=1)
x <- matrix(rnorm(10 * 20), nrow=20)
l <- matrix(rgamma(20, 1, 1), ncol=1)

sum(sapply(1:20, function(i) l[i] * (x[i, ] %*% b)^2))
sum(sapply(1:20, function(i) (x[i, ] %*% b)^2))
# sum(diag(t(x) %*% x %*% b %*% t(b) * diag(as.numeric(l))))

# S_i (x_i' b)^2 = b' x' x b = tr(x' x b b')
sum(sapply(1:20, function(i) (x[i, ] %*% b)^2))
t(b) %*% t(x) %*% x %*% b
sum(diag(t(x) %*% x %*% b %*% t(b)))

x %*% t(x) %*% diag(as.numeric(l))
t(x) %*% diag(as.numeric(l)) %*% x 

sum(sapply(1:20, function(i) l[i] * (x[i, ] %*% b)^2))
t(b) %*% t(x) %*% diag(as.numeric(l)) %*%  x %*% b

t(x) %*% diag(as.numeric(l)) %*% x 
x %*% t(x) %*% diag(as.numeric(l)) %*% x 


X <- t(x) %*% diag(as.numeric(l)) %*% x 
S <- 0
for (i in 1:10) {
    for (j in 1:10) {
	S <- S + X[i, j] * b[i] * b[j]
    }
}
S

rep(1, 10)
x[1, ] %*% 
sum(x[1, ] * b) + sum(1 * b)
sum((x[1, ] + 1) * b)


# -------------------------------------------------------------------------------
# Testing Eq: jj_likelihood
# -------------------------------------------------------------------------------
y <- rbinom(20, 1, 0.5)
X <- matrix(rnorm(10 * 20), nrow=20)
b <- matrix(rnorm(10), ncol=1)
l <- matrix(rgamma(20, 1, 1), ncol=1)
a <- matrix(rgamma(20, 1, 1), ncol=1)

S <- 0
for (i in 1:20) {
    xb <- X[i, ] %*% b
    S <- S + - y[i] * xb + (xb + l[i])/2 + a[i]/2 * (xb^2 - l[i]^2)
}

A <- diag(as.numeric(a))
S1 <- 0.5 * (
    t(b) %*% t(X) %*% A %*% X %*% b -
    t(l) %*% A %*% l +
    sum(X %*% b + l)
) - sum(y * (X %*% b))

abs(S - S1) < 1e-3


# ----------------------------------------
# S_ij (X' A X)_ij E b_i b_j
# ----------------------------------------
P <- t(X) %*% A %*% X
G_k <- 1:5
G_j <- 6:10

S <- 0
for (i in G_k) {
    for (j in G_k) {
	S <- S + P[i, j] * b[i] * b[j]
    }
}
S
t(b[G_k]) %*% P[G_k, G_k] %*% b[G_k]
t(b[G_k]) %*% t(X[ , G_k]) %*% A %*% X[ , G_k] %*% b[G_k]

S <- 0
for (i in G_k) {
    for (j in G_j) {
	S <- S + P[i, j] * b[i] * b[j]
    }
}
S
t(b[G_k]) %*% t(X[ , G_k]) %*% A %*% X[ , G_j] %*% b[G_j]


# ----------------------------------------
# Testing update mu
# ----------------------------------------
y <- rbinom(20, 1, 0.5)
X <- matrix(rnorm(10 * 20), nrow=20)
b <- matrix(rnorm(10), ncol=1)
l <- matrix(rgamma(20, 1, 1), ncol=1)

m <- rnorm(p)
s <- runif(p, 0.1, 0.2)
g <- rep(0.1, p)
l <- rgamma(n, 1, 1)

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
b0 <- 5
w <- a0 / (a0 + b0)

m <- rep(0.1, p)
s <- rep(0.1, p)
g <- rep(0.5, p)

l <- opt_l(X, m, s, g)
A <- diag(as.numeric(a(l)))
XAX <- t(X) %*% A %*% X
lAl <- t(l) %*% A %*% l
L <- cbind(L, l)

0.5 * t(m_G) %*% XAX[G, G] %*% m_G +
t(m_G) %*% XAX[G, Gc] %*% (g[Gc] * m[Gc]) +
sum((0.5 - y) * X[ , G] %*% m_G) +
lambda * (sum(s[G]^2 + m_G^2))^(1/2)

# ----------------------------------------
# C++ TESTING
# ----------------------------------------
library(Rcpp)
sourceCpp("./jaakkola.cpp")

approxEq <- function(a, b, e=1e-3) all(abs(a - b) < e)


# ----------------------------------------
# Test data
# ----------------------------------------
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


# ----------------------------------------
# Test Algorithm
# ----------------------------------------
gsvb.logistic <- function(y, X, groups, m, s, g, niter=500, fit=NULL)
{
    n <- nrow(X) 
    p <- ncol(X)
    w <- 1 / p * length(unique(groups))

    l <- jaak_update_l(X, m, s, g)
    A <- diag(as.numeric(a(l)))
    XAX <- t(X) %*% A %*% X

    for (iter in 1:niter)
    {
	m.old <- m; s.old <- s; g.old <- g

	l <- jaak_update_l(X, m, s, g)
	A <- diag(as.numeric(a(l)))
	XAX <- t(X) %*% A %*% X

	for (group in unique(groups))
	{
	    G  <- which(groups == group) - 1
	    Gc <- which(groups != group) - 1
	    
	    m[G+1] <- jaak_update_mu(y, X, XAX, m, s, g, 1, G, Gc)
	    s[G+1] <- jaak_update_s( y, XAX, m, s, 1, G)
	    g[G+1] <- jaak_update_g( y, X, XAX, m, s, g, 1, w, G, Gc)
	}
	cat(iter)

	if (all(c(abs(m.old - m) < 1e-3, 
		  abs(s.old - s) < 1e-3, 
		  abs(g.old - g) < 1e-3))) 
	    break
    }
    
    return(list(m=m, s=s, g=g))
}

mu <- rnorm(p); s <- rep(0.1, p); g <- runif(p/gsize)[groups]
mu <- rnorm(p); s <- rep(0.1, p); g <- +!!b
mu <- rnorm(p); s <- rep(0.1, p); g <- rep(0.5, p)
f <- gsvb.logistic(y, X, groups, mu, s, g, niter=200)

plot(f$m * f$g)
plot(f$s)
plot(f$g[!duplicated(groups)])

sqrt(sum(((f$m * f$g) - b)^2))
sum(abs((f$m * f$g) - b))


# ----------------------------------------
# Test update mu [PASS]
# ----------------------------------------
mu <- rnorm(p); s <- rep(0.1, p); g <- runif(p/gsize)[groups]
mu <- rep(0.1, p); s <- rep(0.1, p); g <- rep(0.5, p)
G  <- which(groups == 3) - 1
Gc <- which(groups != 3) - 1

l <- opt_l(X, mu, s, g)
A <- diag(as.numeric(a(l)))
XAX <- t(X) %*% A %*% X

jaak_update_mu(y, X, XAX, mu, s, g, 1, G, Gc), 
t(t(optim(mu[G+1], 
    fn=function(mG) opt_mu(mG, y, X, XAX, mu, s, g, G+1, Gc+1, 1),
    control=list(maxit=20),
    method="BFGS")$par))


# ----------------------------------------
# Test update s [PASS]
# ----------------------------------------
jaak_update_s(y, XAX, mu, s, 1, G)

optim(s[G+1], 
    fn=function(sG) opt_s(sG, XAX, mu, G+1, 1),
    control=list(maxit=20),
    method="L-BFGS-B", lower=1e-3, upper=s[G+1][1] + 0.2)$par


# ----------------------------------------
# Test update g [PASS]
# ----------------------------------------
jaak_update_g(y, X, XAX, mu, s, g, 1, 1/200, G, Gc)
opt_g(y, X, XAX, mu, s, g, G+1, Gc+1, 1, 1/200)


# ----------------------------------------
# Test update l [PASS]
# ----------------------------------------
jaak_update_l(X, mu, s, g)
opt_l(X, mu, s, g)

