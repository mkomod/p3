# ----------------------------------------
# Setup
# ----------------------------------------
n <- 100
p <- 200

y <- matrix(rnorm(n), ncol=1)
X <- matrix(rnorm(n * p), nrow=n)
b <- matrix(rnorm(p), ncol=1)

approxEq <- function(a, b, e=1e-3) abs(a - b) < e

# ----------------------------------------
# Matrix identities
# ----------------------------------------
# || Xb || = b' X' X b = tr(X'X bb')

xb <- X %*% b
bbt <- b %*% t(b)
xtx <- crossprod(X)

approxEq(t(xb) %*% xb, t(b) %*% xtx %*% b)
approxEq(t(b) %*% xtx %*% b, sum(diag(xtx %*% bbt)))


sum(sapply(1:p, function(i) {
    sapply(1:p, function(j) xtx[i, j] * b[i] * b[j])
}))

t(b) %*% xtx %*% b


# ----------------------------------------
# Vector identities
# ----------------------------------------
# <y, X b> = sum <y, Xi bi>

approxEq(crossprod(X %*% b, y), 
	 sum(sapply(1:p, function(i) crossprod(X[ , i] * b[i], y))))


approxEq(
t(b[1:5]) %*% xtx[1:5, 8:12] %*% b[8:12],
sum(sapply(1:5, function(i) {
    sapply(8:12, function(j) xtx[i, j] * b[i] * b[j])
}))
)


crossprod(crossprod(xtx, b), b)

# ----------------------------------------
# Variational family identities
# ----------------------------------------
g  <- c(1, 0.5)
mu <- c(2, 4)
s  <- c(2, 2)

b <- sapply(1:2, function(i) rnorm(1e5, mu[i], s[i]) * (runif(1e5) < g[i]))

mean(b[ , 1] * b[ , 2]) - mean(b[ , 1]) * mean(b[ , 2])
cov(b)



# ----------------------------------------
# l2 of normal vector
# ----------------------------------------
mu <- c(2, 9)
s  <- c(2, 4)
mcn <- 5e2
a <- matrix(rnorm(mcn * length(mu), mean=mu, sd=s), nrow=length(mu))
norm(a, type="F") / sqrt(mcn)

microbenchmark::microbenchmark(
    mean(apply(a, 2, l2)),
    norm(a, type="F") / sqrt(mcn)
)


# ----------------------------------------
# Constrained optimization
# ----------------------------------------
curve((exp(x) - 3) *  exp(x), -1.5, 1.5)
curve((x - 3) *  x, -5, 5)
optim(1, fn=function(x) (x - 3) * x)
optim(1, fn=function(x) (exp(x) - 3) * exp(x), 
      gr=function(x) ((exp(x) - 3) + exp(x)) * exp(x), method="CG")



# ----------------------------------------
# Minima
# ----------------------------------------
f <- function(x, S=1, ap=1, a=1e-3, n=500, b=1e-3) 
{
    ((b + S/2) * ap) / x + (n/2 + a) * log(x)
}

xs <- seq(1e-3, 1e4, by=1e-1)
plot(f(xs, 100))
min(f(xs, 100))


# ----------------------------------------
# Sampling from spike and slab
# ----------------------------------------
library(mvtnorm)

mk <- 4
sk <- 4
gk <- 0.5
n <- 1e6
samples <- rnorm(n, mk, sk) * (runif(n) < gk)

mean(1/(2 * sk^2) * (samples - mk)^2)

gk * 0.5 - gk * mk^2 * 0.5 / sk^2 + mk^2 * 0.5 / sk^2


# ----------------------------------------
# S_i S_j x_i x_j b_i b_j = b' x x' b
# ----------------------------------------
x <- rnorm(50)
b <- rnorm(50)

S <- 0
for(i in 1:50) {
    for (j in 1:50) {
	S <- S + x[i]*x[j]*b[i]*b[j]
    }
}
t(b) %*% x %*% t(x) %*% b


# ----------------------------------------
# optim
# ----------------------------------------
optim(10, fn=function(x) x^2)


# ----------------------------------------
# (x' b)^2
# ----------------------------------------
x <- rnorm(10)
b <- rnorm(10)

approxEq(sum(outer(x, x) * outer(b, b)), (sum(x * b))^2)


