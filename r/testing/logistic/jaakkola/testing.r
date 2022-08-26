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

